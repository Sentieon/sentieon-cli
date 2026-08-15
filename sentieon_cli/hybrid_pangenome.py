"""
The Sentieon hybrid-pangenome pipeline

Combines the human pangenome with sample short-read and long-read data
for highly accurate variant calling.
"""

import argparse
import copy
import json
import pathlib
import shutil
import sys
from typing import Dict, List, Optional, Set

import packaging.version

from . import command_strings as cmds
from .archive import ar_load
from .base_pangenome import BasePangenome
from .dag import DAG
from .driver import (
    DNAModelApply,
    DNAscope,
    Driver,
    LongReadSV,
    PangenomeSV,
    PGHapUpdateAlgo,
)
from .job import Job
from .logging import get_logger
from .shard import (
    GRCH38_CONTIGS,
    determine_shards_from_fai,
    parse_fai,
    vcf_contigs,
)
from .shell_pipeline import Command, Pipeline
from .transfer import build_transfer_jobs
from .util import (
    __version__,
    check_kmc_patch,
    check_version,
    parse_rg_line,
    path_arg,
    total_memory,
    vcf_id,
)

HYBRID_PANGENOME_MIN_VERSIONS = {
    "kmc": None,
    "sentieon driver": packaging.version.Version("202503.04"),
    "vg": None,
    "bcftools": packaging.version.Version("1.22"),
    "samtools": packaging.version.Version("1.16"),
    "bedtools": None,
}

# LongReadSV settings for finding graph update regions
LONGREADSV_MIN_SV_SIZE = 20

# PangenomeSV settings
PANGENOME_SV_MIN_AF = 0.1

logger = get_logger(__name__)


class HybridPangenome(BasePangenome):
    """The Sentieon hybrid-pangenome pipeline"""

    params = copy.deepcopy(BasePangenome.params)
    params.update(
        {
            # Required arguments
            "lr_aln": {
                "nargs": "*",
                "help": (
                    "Long-read BAM or CRAM files aligned to the linear "
                    "reference genome."
                ),
                "required": True,
                "type": path_arg(exists=True, is_file=True),
            },
            "pop_vcf": {
                "flags": ["--pop_vcf"],
                "help": (
                    "A VCF containing annotations for use with DNAModelApply."
                ),
                "type": path_arg(exists=True, is_file=True),
                "required": True,
            },
            "readgroup": {
                "help": "Readgroup information for the fastq files.",
            },
            # Additional arguments
            "bed": {
                "flags": ["-b", "--bed"],
                "help": (
                    "Region BED file. Supplying this file will limit "
                    "small-variant calling to the intervals inside the BED "
                    "file."
                ),
                "type": path_arg(exists=True, is_file=True),
            },
            "pangenome_ref_name": {
                "default": "GRCh38",
                "help": "Reference name in the pangenome (GRCh38).",
            },
            "rgsm": {
                "help": (
                    "Overwrite the SM tag of the input readgroups for "
                    "compatibility"
                ),
            },
            "skip_metrics": {
                "help": "Skip metrics collection and multiQC",
                "action": "store_true",
            },
            "skip_multiqc": {
                "help": "Skip multiQC report generation",
                "action": "store_true",
            },
            # Hidden arguments
            "skip_contig_checks": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_model_apply": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_pangenome_name_checks": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_pop_vcf_id_check": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_small_variants": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
            "skip_svs": {
                "help": argparse.SUPPRESS,
                "action": "store_true",
            },
        }
    )

    positionals = BasePangenome.positionals

    def __init__(self) -> None:
        super().__init__()
        self.lr_aln: List[pathlib.Path] = []
        self.pop_vcf: Optional[pathlib.Path] = None
        self.readgroup: Optional[str] = None
        self.bed: Optional[pathlib.Path] = None
        self.pangenome_ref_name = "GRCh38"
        self.rgsm: Optional[str] = None
        self.extract_model_name = "extract.model"
        self.skip_metrics = False
        self.skip_multiqc = False
        self.skip_contig_checks = False
        self.skip_model_apply = False
        self.skip_pangenome_name_checks = False
        self.skip_pop_vcf_id_check = False
        self.skip_small_variants = False
        self.skip_svs = False

    def validate(self) -> None:
        """Validate pipeline inputs"""
        self.validate_ref()
        self.fai_data = parse_fai(pathlib.Path(str(self.reference) + ".fai"))
        self.shards = determine_shards_from_fai(
            self.fai_data, 10 * 1000 * 1000
        )
        self.pop_vcf_contigs: Dict[str, Optional[int]] = {}
        if self.pop_vcf:
            self.pop_vcf_contigs = vcf_contigs(self.pop_vcf, self.dry_run)
            self.logger.debug("VCF contigs are: %s", self.pop_vcf_contigs)

        self.validate_bundle()
        self.validate_output_vcf()

        if not self.r1_fastq or not self.readgroup:
            self.logger.error(
                "Please supply the short reads with the `--r1_fastq`, "
                "`--r2_fastq`, and `--readgroup` arguments"
            )
            sys.exit(2)
        if len(self.r1_fastq) != len(self.r2_fastq):
            self.logger.error(
                "The number of input `--r1_fastq` files does not equal the "
                "number of `--r2_fastq` files"
            )
            sys.exit(2)
        if not self.lr_aln:
            self.logger.error(
                "Please supply the long-read alignments with the `--lr_aln` "
                "argument"
            )
            sys.exit(2)

        self.validate_bwa_index()
        self.collect_readgroups()
        self.validate_readgroups()

        if not self.skip_version_check:
            for cmd, min_version in HYBRID_PANGENOME_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

            if not check_kmc_patch("kmc"):
                self.logger.error(
                    "Error: The 'kmc' executable in the PATH does not "
                    "support reading from stdin. Please ensure "
                    "you are using the patched version of KMC from "
                    "https://github.com/Sentieon/KMC/releases."
                )
                sys.exit(2)

        if self.bed is None:
            self.logger.info(
                "A BED file is recommended to avoid small-variant calling "
                "across decoy and unplaced contigs."
            )

        if not self.skip_pangenome_name_checks:
            if not str(self.gbz).endswith("grch38.gbz"):
                self.logger.error(
                    "The `--gbz` file does not have the expected suffix. "
                    "Check that you are using a GRCh38 pangenome."
                )
                sys.exit(2)

            if not str(self.hapl).endswith("grch38.hapl"):
                self.logger.error(
                    "The `--hapl` file does not have the expected suffix. "
                    "Check that you are using a GRCh38 pangenome."
                )
                sys.exit(2)

        if not self.skip_contig_checks:
            # Check the fai file contigs
            mismatch_contigs: Set[str] = set()
            for ctg, length in GRCH38_CONTIGS.items():
                d = self.fai_data.get(ctg, {})
                fai_length = d.get("length", -1)
                if length != fai_length:
                    mismatch_contigs.add(ctg)
            if mismatch_contigs:
                mismatch_contigs_s = ", ".join(mismatch_contigs)
                self.logger.error(
                    "Reference contigs with unexpected lengths: %s",
                    mismatch_contigs_s,
                )
                sys.exit(2)

            # Check the pop VCF file contigs
            if not self.dry_run:
                mismatch_contigs = set()
                for ctg, length in GRCH38_CONTIGS.items():
                    vcf_length = self.pop_vcf_contigs.get(ctg, -1)
                    if length != vcf_length:
                        mismatch_contigs.add(ctg)
                if mismatch_contigs:
                    mismatch_contigs_s = ", ".join(mismatch_contigs)
                    self.logger.error(
                        "Pop VCF contigs with unexpected lengths: %s",
                        mismatch_contigs_s,
                    )
                    sys.exit(2)

    def validate_bundle(self) -> None:
        """Validate the model bundle"""
        bundle_info_bytes = ar_load(
            str(self.model_bundle) + "/bundle_info.json"
        )
        if isinstance(bundle_info_bytes, list):
            bundle_info_bytes = b"{}"
        bundle_info = json.loads(bundle_info_bytes.decode())

        req_version_s = bundle_info.get("minScriptVersion")
        if req_version_s:
            req_version = packaging.version.Version(req_version_s)
            if req_version > packaging.version.Version(__version__):
                self.logger.error(
                    "The model bundle requires version %s or later of the "
                    "sentieon-cli.",
                    req_version,
                )
                sys.exit(2)

        bundle_pipeline = bundle_info.get("pipeline")
        if bundle_pipeline and bundle_pipeline != "Hybrid pangenome":
            self.logger.error("The model bundle is for a different pipeline.")
            sys.exit(2)

        bundle_members = set(ar_load(str(self.model_bundle)))

        # Prefer a reference-specific extract model. Fall back to the generic
        # 'extract.model' only for the default 'GRCh38' reference.
        extract_candidate = f"extract.{self.pangenome_ref_name}.model"
        if extract_candidate in bundle_members:
            self.extract_model_name = extract_candidate
        elif self.pangenome_ref_name == "GRCh38":
            self.extract_model_name = "extract.model"
        else:
            self.extract_model_name = extract_candidate

        required_members = {
            "bwa.model",
            "dnascope.model",
            "longreadsv.model",
            "minimap2.model",
            self.extract_model_name,
        }
        missing_members = required_members - bundle_members
        if missing_members:
            self.logger.error(
                "Expected model files not found in the model bundle file: %s",
                ", ".join(sorted(missing_members)),
            )
            sys.exit(2)

        bundle_vcf_id = bundle_info.get("SentieonVcfID")
        if (
            bundle_vcf_id
            and not self.skip_pop_vcf_id_check
            and not self.dry_run
        ):
            assert self.pop_vcf is not None
            pop_vcf_id = vcf_id(self.pop_vcf)
            if bundle_vcf_id != pop_vcf_id:
                self.logger.error(
                    "The ID of the `--pop_vcf` does not match the model bundle"
                )
                sys.exit(2)

    def collect_readgroups(self) -> None:
        """Collect readgroup tags from the inputs"""
        assert self.readgroup is not None
        try:
            parsed_rg = parse_rg_line(self.readgroup.replace(r"\t", "\t"))
        except ValueError as e:
            self.logger.error(
                "Invalid --readgroup value '%s': %s", self.readgroup, e
            )
            sys.exit(2)
        if not parsed_rg.get("ID"):
            self.logger.error(
                "Readgroup '%s' does not have a RGID tag",
                self.readgroup,
            )
            sys.exit(2)
        if not parsed_rg.get("SM"):
            self.logger.error(
                "Readgroup '%s' does not have a RGSM tag",
                self.readgroup,
            )
            sys.exit(2)
        self.fastq_readgroup: Dict[str, str] = parsed_rg

        self.lr_readgroups: List[List[Dict[str, str]]] = []
        for aln in self.lr_aln:
            self.lr_readgroups.append([])
            for rg_line in cmds.get_rg_lines(aln, self.dry_run):
                self.lr_readgroups[-1].append(parse_rg_line(rg_line))

    def validate_readgroups(self) -> None:
        """Confirm that all readgroups have a consistent SM tag"""
        rg_sm = self.fastq_readgroup.get("SM")
        for aln, aln_rgs in zip(self.lr_aln, self.lr_readgroups):
            for aln_rg in aln_rgs:
                if not aln_rg.get("ID"):
                    self.logger.error(
                        "Found a readgroup without an ID tag in '%s': %s",
                        aln,
                        str(aln_rg),
                    )
                    sys.exit(2)
                sm = aln_rg.get("SM")
                if not sm:
                    self.logger.error(
                        "Found a readgroup without a SM tag in '%s': %s",
                        aln,
                        str(aln_rg),
                    )
                    sys.exit(2)
                if self.dry_run or self.rgsm:
                    continue
                if sm != rg_sm:
                    self.logger.error(
                        "Input readgroup '%s' has a different RG-SM tag "
                        "from the `--readgroup` argument.\n"
                        "found='%s' expected='%s'. Please set the `--rgsm` "
                        "argument to override the SM tag in the input files",
                        str(aln_rg),
                        sm,
                        rg_sm,
                    )
                    sys.exit(2)
        self.sample_sm: str = self.rgsm if self.rgsm else str(rg_sm)

    def configure(self) -> None:
        """Configure pipeline parameters"""
        pass

    def find_unzip(self) -> str:
        """The decompression tool for fastq input"""
        unzip = "igzip"
        if not shutil.which(unzip):
            self.logger.info(
                "igzip is recommended for decompression, but is not "
                "available. Falling back to gzip."
            )
            unzip = "gzip"
        return unzip

    def build_dag(self) -> DAG:
        """Build the DAG for the hybrid-pangenome pipeline"""
        if not self.reference:
            self.logger.error("reference is required")
            sys.exit(2)
        if not self.model_bundle:
            self.logger.error("model_bundle is required")
            sys.exit(2)
        if not self.output_vcf:
            self.logger.error("output_vcf is required")
            sys.exit(2)
        if not self.pop_vcf:
            self.logger.error("pop_vcf is required")
            sys.exit(2)

        self.logger.info("Building the hybrid-pangenome DAG")
        dag = DAG()

        ref_fai = pathlib.Path(str(self.reference) + ".fai")

        # Output files
        suffix = "bam" if self.bam_format else "cram"
        out_bwa_aln = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", f"_bwa_deduped.{suffix}")
        )
        out_lift_aln = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", f"_lift_deduped.{suffix}")
        )
        sv_vcf = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_sv.vcf.gz")
        )

        # Intermediate file paths
        bwa_bam = self.tmp_dir.joinpath("sample-bwa.bam")
        lift_bam = self.tmp_dir.joinpath("sample-lift.bam")
        ext_fastq = self.tmp_dir.joinpath("sample-extract.fq.gz")
        kmer_prefix = self.tmp_dir.joinpath("sample")
        kmer_file = pathlib.Path(str(kmer_prefix) + ".kff")
        sample_pangenome = self.tmp_dir.joinpath("sample_pangenome.gbz")
        hap_raw_gfa = self.tmp_dir.joinpath("sample-hap.raw.gfa")
        pangenome_raw_gfa = self.tmp_dir.joinpath("sample-pangenome-raw.gfa")
        pangenome_gfa = self.tmp_dir.joinpath("sample-pangenome.gfa")
        pangenome_fasta = self.tmp_dir.joinpath("sample-pangenome.fa")
        longread_sv_vcf = self.tmp_dir.joinpath("sample-longread-sv.vcf.gz")
        sv_bed = self.tmp_dir.joinpath("sample-sv.bed")
        raw_vcf = self.tmp_dir.joinpath("sample-dnascope.vcf.gz")
        transfer_vcf = self.tmp_dir.joinpath("sample-dnascope_transfer.vcf.gz")

        total_mem_gb = total_memory() / (1024.0**3)

        # KMC k-mer counting across the short and long reads
        kmc_job = self.build_hybrid_kmc_job(kmer_prefix)
        dag.add_job(kmc_job)
        haplotype_dependencies: Set[Job] = {kmc_job}

        # BWA alignment and extraction
        bwa_job = self.build_alignment_job(bwa_bam, ext_fastq)
        dag.add_job(bwa_job)
        # Do not run vg-haplotypes with bwa in low-mem environments
        if total_mem_gb < 70:
            haplotype_dependencies.add(bwa_job)

        # vg haplotypes - create a sample-specific pangenome
        haplotypes_job = self.build_haplotypes_job(sample_pangenome, kmer_file)
        dag.add_job(haplotypes_job, haplotype_dependencies)

        # convert the sample pangenome to GFA
        gfa_job = self.build_gfa_job(hap_raw_gfa, sample_pangenome)
        dag.add_job(gfa_job, {haplotypes_job})

        # Graph update without the SV BED
        update_raw_job = self.build_graph_update_job(
            pangenome_raw_gfa,
            hap_raw_gfa,
            name="graph-update-raw",
        )
        dag.add_job(update_raw_job, {gfa_job})

        # Call SVs from the long reads and collect graph update regions
        longreadsv_job = self.build_longreadsv_job(longread_sv_vcf)
        dag.add_job(longreadsv_job)
        sv_bed_job = Job(
            cmds.cmd_longread_sv_bed(sv_bed, longread_sv_vcf),
            "longread-sv-bed",
            0,
            task_name="pangenome-update",
        )
        dag.add_job(sv_bed_job, {longreadsv_job})

        # Graph update with the SV BED
        update_job = self.build_graph_update_job(
            pangenome_gfa,
            pangenome_raw_gfa,
            bed=sv_bed,
            name="graph-update",
        )
        dag.add_job(update_job, {update_raw_job, sv_bed_job})

        # FASTA generation from the updated graph
        gfa2fa_job = Job(
            cmds.cmd_pgutil_gfa2fa(pangenome_fasta, ref_fai, pangenome_gfa),
            "gfa2fa",
            0,
            task_name="pangenome",
        )
        dag.add_job(gfa2fa_job, {update_job})
        faidx_job = Job(
            cmds.cmd_samtools_faidx(pangenome_fasta),
            "faidx",
            0,
            task_name="pangenome",
        )
        dag.add_job(faidx_job, {gfa2fa_job})

        # minimap2 alignment of the extracted reads with liftover
        mm2_job = self.build_minimap2_lift_job(
            lift_bam,
            ext_fastq,
            pangenome_fasta,
            pangenome_gfa,
        )
        dag.add_job(mm2_job, {bwa_job, faidx_job, update_job})

        # Deduplicate the short-read alignments. The `--lr_aln` input is
        # assumed to be deduplicated already.
        # Emit Dedup metrics for the primary (bwa) short-read alignment so
        # they land in the metrics directory scanned by MultiQC.
        dedup_metrics: Optional[pathlib.Path] = None
        if not self.skip_metrics:
            metrics_dir = pathlib.Path(
                str(self.output_vcf).replace(".vcf.gz", "_metrics")
            )
            if not self.dry_run:
                metrics_dir.mkdir(exist_ok=True)
            sample_name = self.output_vcf.name.replace(".vcf.gz", "")
            dedup_metrics = metrics_dir.joinpath(
                sample_name + ".txt.dedup_metrics.txt"
            )

        bwa_lc_job, bwa_dedup_job = self.build_dedup_job(
            out_bwa_aln, [bwa_bam], "bwa", metrics=dedup_metrics
        )
        lift_lc_job, lift_dedup_job = self.build_dedup_job(
            out_lift_aln,
            [lift_bam],
            "lift",
        )
        dag.add_job(bwa_lc_job, {bwa_job})
        dag.add_job(bwa_dedup_job, {bwa_lc_job})
        dag.add_job(lift_lc_job, {mm2_job})
        dag.add_job(lift_dedup_job, {lift_lc_job})

        # Alignment metrics from the deduplicated short reads
        if not self.skip_metrics:
            metrics_job, rehead_job = self.build_metrics_job(
                [out_bwa_aln, out_lift_aln],
            )
            dag.add_job(metrics_job, {bwa_dedup_job, lift_dedup_job})
            dag.add_job(rehead_job, {metrics_job})
            if not self.skip_multiqc:
                multiqc_job = self.multiqc()
                if multiqc_job:
                    dag.add_job(multiqc_job, {rehead_job})

        # Variant calling with the original, lifted, and long reads
        calling_bams = [out_bwa_aln, out_lift_aln] + list(self.lr_aln)
        replace_rg = self.build_replace_rg()
        calling_dependencies = {bwa_dedup_job, lift_dedup_job}

        if not self.skip_svs:
            pangenomesv_job = self.build_pangenomesv_job(
                sv_vcf,
                calling_bams,
                pangenome_gfa,
                replace_rg,
            )
            dag.add_job(pangenomesv_job, calling_dependencies | {update_job})

        if self.skip_small_variants:
            return dag

        dnascope_job = self.build_dnascope_job(
            raw_vcf,
            calling_bams,
            replace_rg,
        )
        dag.add_job(dnascope_job, calling_dependencies)

        # transfer annotations from the pop_vcf
        transfer_target = (
            transfer_vcf if not self.skip_model_apply else self.output_vcf
        )
        transfer_jobs, concat_job = build_transfer_jobs(
            transfer_target,
            self.pop_vcf,
            raw_vcf,
            self.tmp_dir,
            self.shards,
            self.pop_vcf_contigs,
            self.fai_data,
            self.dry_run,
            self.cores,
        )
        for job in transfer_jobs:
            dag.add_job(job, {dnascope_job})
        dag.add_job(concat_job, set(transfer_jobs))

        if not self.skip_model_apply:
            apply_job = self.build_dnamodelapply_job(
                transfer_vcf, self.output_vcf
            )
            dag.add_job(apply_job, {concat_job})

        return dag

    def build_hybrid_kmc_job(self, kmer_prefix: pathlib.Path) -> Job:
        """Count k-mers across the short-read fastq and long-read
        alignments"""
        assert self.reference is not None
        kmc_job = Job(
            cmds.cmd_hybrid_kmc(
                kmer_prefix,
                self.r1_fastq + self.r2_fastq,
                self.lr_aln,
                self.reference,
                self.tmp_dir,
                memory=self.kmer_memory,
                threads=self.cores,
                unzip=self.find_unzip(),
            ),
            "kmc",
            0,  # run in the background
            task_name="kmer-counting",
        )
        return kmc_job

    def build_alignment_job(
        self,
        sample_bam: pathlib.Path,
        sample_fastq: pathlib.Path,
    ) -> Job:
        """Build the bwa alignment and read extraction job"""
        assert self.reference is not None
        assert self.model_bundle is not None

        rg = copy.deepcopy(self.fastq_readgroup)
        rg["SM"] = self.sample_sm
        rg["LR"] = "0"
        bwa_job = Job(
            cmds.cmd_bwa_extract(
                sample_bam,
                sample_fastq,
                self.reference,
                self.r1_fastq,
                self.r2_fastq,
                "@RG\\t" + "\\t".join([f"{x[0]}:{x[1]}" for x in rg.items()]),
                self.model_bundle.joinpath(self.extract_model_name),
                self.model_bundle.joinpath("bwa.model"),
                self.cores,
                unzip=self.find_unzip(),
            ),
            "bwa-extract",
            self.cores,
            task_name="alignment",
        )
        return bwa_job

    def build_haplotypes_job(
        self, output_gbz: pathlib.Path, kmer_file: pathlib.Path
    ) -> Job:
        """Build vg haplotypes job"""
        assert self.hapl is not None
        assert self.gbz is not None

        haplotypes_job = Job(
            cmds.cmd_vg_haplotypes(
                output_gbz,
                kmer_file,
                self.hapl,
                self.gbz,
                threads=self.cores,
                xargs=[
                    "--include-reference",
                    "--diploid-sampling",
                    "--set-reference",
                    self.pangenome_ref_name,
                ],
            ),
            "vg-haplotypes",
            self.cores,
            task_name="pangenome",
        )
        return haplotypes_job

    def build_gfa_job(
        self, output_gfa: pathlib.Path, input_gbz: pathlib.Path
    ) -> Job:
        """Build vg convert to GFA job"""
        gfa_job = Job(
            cmds.cmd_vg_convert_gfa(
                output_gfa,
                input_gbz,
                threads=self.cores,
                reference_name=self.pangenome_ref_name,
            ),
            "vg-convert-gfa",
            0,
            task_name="pangenome",
        )
        return gfa_job

    def build_graph_update_job(
        self,
        out_gfa: pathlib.Path,
        in_gfa: pathlib.Path,
        bed: Optional[pathlib.Path] = None,
        name: str = "graph-update",
    ) -> Job:
        """Update the personalized graph using the long-read alignments"""
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=list(self.lr_aln),
        )
        driver.add_algo(
            PGHapUpdateAlgo(out_gfa, gfa_file=in_gfa, target_bed=bed)
        )
        return Job(
            Pipeline(Command(*driver.build_cmd())),
            name,
            self.cores,
            task_name="pangenome-update",
        )

    def build_longreadsv_job(self, out_vcf: pathlib.Path) -> Job:
        """Call SVs from the long reads for the graph update"""
        assert self.model_bundle is not None
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=list(self.lr_aln),
        )
        driver.add_algo(
            LongReadSV(
                out_vcf,
                model=self.model_bundle.joinpath("longreadsv.model"),
                min_sv_size=LONGREADSV_MIN_SV_SIZE,
            )
        )
        return Job(
            Pipeline(Command(*driver.build_cmd())),
            "longreadsv",
            self.cores,
            task_name="pangenome-update",
        )

    def build_minimap2_lift_job(
        self,
        out_bam: pathlib.Path,
        ext_fastq: pathlib.Path,
        pangenome_fasta: pathlib.Path,
        pangenome_gfa: pathlib.Path,
    ) -> Job:
        """Align the extracted reads to the personalized pangenome and lift
        the alignments back to the linear reference"""
        assert self.reference is not None
        assert self.model_bundle is not None

        rg = copy.deepcopy(self.fastq_readgroup)
        rg["ID"] = rg["ID"] + "-pg"
        rg["SM"] = self.sample_sm
        rg["LR"] = "2"
        mm2_job = Job(
            cmds.cmd_minimap2_lift(
                out_bam,
                pangenome_fasta,
                ext_fastq,
                pangenome_gfa,
                self.reference,
                "@RG\\t" + "\\t".join([f"{x[0]}:{x[1]}" for x in rg.items()]),
                self.model_bundle.joinpath("minimap2.model"),
                threads=self.cores,
                mm2_xargs=["--secondary=yes"],
            ),
            "mm2-lift",
            self.cores,
            task_name="pangenome-alignment",
        )
        return mm2_job

    def build_replace_rg(self) -> List[List[str]]:
        """`--replace_rg` arguments setting the LR readgroup attribute for
        the variant calling stages.

        The bwa (LR:0) and lifted (LR:2) alignments are generated by the
        pipeline with the LR attribute already in place; only the input
        long-read readgroups are rewritten (LR:1), as they cannot be
        assumed to carry the attribute.
        """
        replace_rg: List[List[str]] = [[], []]  # bwa and lifted alignments
        for aln_rgs in self.lr_readgroups:
            replace_rg.append([])
            for aln_rg in aln_rgs:
                rg_id = aln_rg.get("ID")
                sm = self.rgsm if self.rgsm else aln_rg.get("SM")
                replace_rg[-1].append(f"{rg_id}=ID:{rg_id}\\tSM:{sm}\\tLR:1")
        return replace_rg

    def build_pangenomesv_job(
        self,
        out_vcf: pathlib.Path,
        input_bams: List[pathlib.Path],
        pangenome_gfa: pathlib.Path,
        replace_rg: List[List[str]],
    ) -> Job:
        """Call SVs with the original, lifted, and long reads"""
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=input_bams,
            replace_rg=replace_rg,
        )
        driver.add_algo(
            PangenomeSV(
                out_vcf,
                gfa_file=pangenome_gfa,
                min_af=PANGENOME_SV_MIN_AF,
            )
        )
        return Job(
            Pipeline(Command(*driver.build_cmd())),
            "pangenome-sv",
            self.cores,
            task_name="sv-calling",
        )

    def build_dnascope_job(
        self,
        out_vcf: pathlib.Path,
        input_bams: List[pathlib.Path],
        replace_rg: List[List[str]],
    ) -> Job:
        """Call small variants with the original, lifted, and long reads"""
        assert self.model_bundle is not None
        pcr_indel_model = "NONE" if self.pcr_free else "CONSERVATIVE"
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=input_bams,
            interval=self.bed,
            replace_rg=replace_rg,
        )
        driver.add_algo(
            DNAscope(
                out_vcf,
                model=self.model_bundle.joinpath("dnascope.model"),
                pcr_indel_model=pcr_indel_model,
                dbsnp=self.dbsnp,
            )
        )
        return Job(
            Pipeline(Command(*driver.build_cmd())),
            "dnascope-raw",
            self.cores,
            task_name="variant-calling",
        )

    def build_dnamodelapply_job(
        self,
        in_vcf: pathlib.Path,
        out_vcf: pathlib.Path,
    ) -> Job:
        """Apply the DNAscope model"""
        assert self.model_bundle is not None
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            DNAModelApply(
                model=self.model_bundle.joinpath("dnascope.model"),
                vcf=in_vcf,
                output=out_vcf,
            )
        )
        return Job(
            Pipeline(Command(*driver.build_cmd())),
            "model-apply",
            self.cores,
            task_name="model-apply",
        )
