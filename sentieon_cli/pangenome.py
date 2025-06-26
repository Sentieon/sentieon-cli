"""
Pangenome alignment and variant calling pipeline
"""

import argparse
from enum import Enum
import itertools
import json
import multiprocessing as mp
import os
import pathlib
import shlex
import shutil
import sys
from typing import Any, Dict, List, Optional, Tuple

import packaging.version

from importlib_resources import files

from . import command_strings as cmds
from .dag import DAG
from .driver import (
    AlignmentStat,
    BaseDistributionByCycle,
    CNVscope,
    CNVModelApply,
    Dedup,
    DNAModelApply,
    DNAscope,
    Driver,
    GCBias,
    InsertSizeMetricAlgo,
    MeanQualityByCycle,
    LocusCollector,
    QualDistribution,
    Realigner,
    WgsMetricsAlgo,
)
from .job import Job
from .pipeline import BasePipeline
from .util import __version__, check_version, parse_rg_line, path_arg, tmp


PANGENOME_MIN_VERSIONS = {
    "kmc": None,
    "sentieon driver": packaging.version.Version("202503"),
    "vg": None,
    "bcftools": packaging.version.Version("1.10"),
    "samtools": packaging.version.Version("1.16"),
    "run-t1k": None,
}

MULTIQC_MIN_VERSION = {
    "multiqc": packaging.version.Version("1.18"),
}


class SampleSex(Enum):
    FEMALE = 1
    MALE = 2
    UNKNOWN = 3


class PangenomePipeline(BasePipeline):
    """The Pangenome pipeline"""

    params: Dict[str, Dict[str, Any]] = {
        # Required arguments
        "reference": {
            "flags": ["-r", "--reference"],
            "required": True,
            "help": "fasta for reference genome.",
            "type": path_arg(exists=True, is_file=True),
        },
        "gbz": {
            "help": "The pangenome graph file in GBZ format.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "hapl": {
            "help": "The haplotype file.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "snarls": {
            "help": (
                "The vg snarls file for the .gbz file. Can be created "
                "using vg."
            ),
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "xg": {
            "help": "The xg file for the .gbz file. Can be created using vg",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "model_bundle": {
            "flags": ["-m", "--model_bundle"],
            "help": "The model bundle file.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "r1_fastq": {
            "nargs": "*",
            "help": "Sample R1 fastq files.",
            "type": path_arg(exists=True, is_file=True),
        },
        "r2_fastq": {
            "nargs": "*",
            "help": "Sample R2 fastq files.",
            "type": path_arg(exists=True, is_file=True),
        },
        "readgroups": {
            "nargs": "*",
            "help": "Readgroup information for the fastq files.",
        },
        "t1k_hla": {
            "help": "The DNA HLA FASTA file for T1K.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "t1k_kir": {
            "help": "The DNA KIR FASTA file for T1K.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        "expansion_catalog": {
            "help": "An ExpansionHunter variant catalog.",
            "required": True,
            "type": path_arg(exists=True, is_file=True),
        },
        # Additional arguments
        "cores": {
            "flags": ["-t", "--cores"],
            "help": (
                "Number of threads/processes to use. Defaults to all "
                "available."
            ),
            "default": mp.cpu_count(),
        },
        "dbsnp": {
            "flags": ["-d", "--dbsnp"],
            "help": (
                "dbSNP vcf file Supplying this file will annotate variants "
                "with their dbSNP refSNP ID numbers."
            ),
            "type": path_arg(exists=True, is_file=True),
        },
        "dry_run": {
            "help": "Print the commands without running them.",
            "action": "store_true",
        },
        "kmer_memory": {
            "help": "Memory limit for KMC in GB.",
            "default": 128,
            "type": int,
        },
        "known_sites": {
            "nargs": "*",
            "type": path_arg(exists=True, is_file=True),
            "help": "Known sites for realignment.",
        },
        # Hidden arguments
        "skip_version_check": {
            "help": argparse.SUPPRESS,
            "action": "store_true",
        },
    }

    positionals: Dict[str, Dict[str, Any]] = {
        "output_vcf": {
            "help": "Output VCF file. The file name must end in .vcf.gz",
            "type": path_arg(),
        },
    }

    def __init__(self) -> None:
        super().__init__()
        self.output_vcf: Optional[pathlib.Path] = None
        self.reference: Optional[pathlib.Path] = None
        self.gbz: Optional[pathlib.Path] = None
        self.hapl: Optional[pathlib.Path] = None
        self.snarls: Optional[pathlib.Path] = None
        self.xg: Optional[pathlib.Path] = None
        self.r1_fastq: List[pathlib.Path] = []
        self.r2_fastq: List[pathlib.Path] = []
        self.readgroups: List[str] = []
        self.t1k_hla: Optional[pathlib.Path] = None
        self.t1k_kir: Optional[pathlib.Path] = None
        self.expansion_catalog: Optional[pathlib.Path] = None
        self.model_bundle: Optional[pathlib.Path] = None
        self.dbsnp: Optional[pathlib.Path] = None
        self.cores = mp.cpu_count()
        self.kmer_memory = 128
        self.known_sites: List[pathlib.Path] = []
        self.skip_version_check = False

    def main(self, args: argparse.Namespace) -> None:
        """Run the pipeline"""
        self.handle_arguments(args)
        self.setup_logging(args)
        self.validate()
        self.configure()

        tmp_dir_str = tmp()
        self.tmp_dir = pathlib.Path(tmp_dir_str)

        dag = self.build_first_dag()
        executor = self.run(dag)
        self.check_execution(dag, executor)

        self.get_sex()
        dag = self.build_second_dag()
        executor = self.run(dag)
        self.check_execution(dag, executor)

        if not self.retain_tmpdir:
            shutil.rmtree(tmp_dir_str)

    def build_dag(self) -> DAG:
        return DAG()

    def validate(self) -> None:
        """Validate pipeline inputs"""
        assert self.output_vcf
        assert self.reference

        if not self.r1_fastq:
            self.logger.error("Please supply --r1_fastq arguments")
            sys.exit(2)

        if not str(self.output_vcf).endswith(".vcf.gz"):
            self.logger.error("The output file should end with '.vcf.gz'")
            sys.exit(2)

        if len(self.r1_fastq) != len(self.readgroups):
            self.logger.error(
                "The number of readgroups does not equal the number of fastq "
                "files"
            )
            sys.exit(2)

        # Confirm the presence of the reference index file
        fai_file = str(self.reference) + ".fai"
        if not os.path.isfile(fai_file):
            self.logger.error(
                "Fasta index file %s does not exist. Please index the "
                "reference genome with 'samtools faidx'",
                fai_file,
            )
            sys.exit(2)

        # Confirm the known sites files are indexed
        for sites_file in self.known_sites:
            idx_suffix = (
                ".tbi" if str(sites_file).endswith(".vcf.gz") else ".idx"
            )
            idx_file = str(sites_file) + idx_suffix
            if not os.path.isfile(idx_file):
                self.logger.error(
                    "The index for the known sites file %s does not exist. "
                    "Please index the VCF file.",
                    str(sites_file),
                )
                sys.exit(2)

        # Validate readgroups
        rg_sample = None
        for rg in self.readgroups:
            rg_dict = parse_rg_line(rg.replace(r"\t", "\t"))
            rg_sm = rg_dict.get("SM")
            if rg_sm:
                self.logger.error(
                    "Found a readgroup without a SM tag: %s",
                    str(rg),
                )
            if rg_sample and rg_sample != rg_sm:
                self.logger.error(
                    "Inconsistent readgroup sample information found in: %s",
                    str(rg),
                )
            rg_sample = rg_sm
            if "ID" not in rg_dict:
                self.logger.error(
                    "Found a readgroup without an ID tag: %s",
                    str(rg),
                )

        if not self.skip_version_check:
            for cmd, min_version in PANGENOME_MIN_VERSIONS.items():
                if not check_version(cmd, min_version):
                    sys.exit(2)

    def configure(self) -> None:
        """Configure pipeline parameters"""
        pass

    def build_first_dag(self) -> DAG:
        """Build the first DAG for the pangenome pipeline"""
        self.logger.info("Building the first pangenome DAG")
        dag = DAG()

        # Find the sample name
        sample_name = parse_rg_line(self.readgroups[0].replace(r"\t", "\t"))[
            "SM"
        ]

        # Output files
        out_svs = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_svs.vcf.gz")
        )
        self.realigned_cram = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_pangenome-aligned.cram")
        )
        self.ploidy_json = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_ploidy.json")
        )
        out_hla = pathlib.Path(str(self.output_vcf).replace(".vcf.gz", "_hla"))
        out_kir = pathlib.Path(str(self.output_vcf).replace(".vcf.gz", "_kir"))
        out_segdup = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_segdups")
        )

        # Intermediate file paths
        kmer_prefix = self.tmp_dir.joinpath("sample.fq")
        kmer_file = pathlib.Path(str(kmer_prefix) + ".kff")
        sample_pangenome = self.tmp_dir.joinpath("sample_pangenome.gbz")
        sample_pack = self.tmp_dir.joinpath("sample_coverage.pack")
        lc_score = self.tmp_dir.joinpath("sample_lc_score.txt.gz")
        deduped_bam = self.tmp_dir.joinpath("sample_deduped.bam")
        dnascope_tmp = self.tmp_dir.joinpath("sample_dnascope_tmp.vcf.gz")

        # Create an alignment head for surjection
        alignment_headers = self.create_alignment_headers()

        # KMC k-mer counting
        kmc_job = self.build_kmc_job(kmer_prefix)
        dag.add_job(kmc_job)

        # vg haplotypes - create sample-specific pangenome
        haplotypes_job = self.build_haplotypes_job(sample_pangenome, kmer_file)
        dag.add_job(haplotypes_job, set([kmc_job]))

        # vg giraffe - align reads to pangenome
        sample_gams, giraffe_jobs = self.build_alignment_jobs(sample_pangenome)
        alignment_dependencies = {haplotypes_job}
        for job in giraffe_jobs:
            dag.add_job(job, alignment_dependencies)

        # vg pack - compute read support
        pack_job = self.build_pack_job(sample_pack, sample_gams)
        dag.add_job(pack_job, set(giraffe_jobs))

        # vg call - call SVs
        call_job = self.build_call_job(out_svs, sample_pack, sample_name)
        dag.add_job(call_job, {pack_job})

        # surject and sort - move reads back to the linear reference
        aligned_bams, surject_jobs = self.build_surject_jobs(
            sample_gams, alignment_headers
        )
        for job in surject_jobs:
            dag.add_job(job, set(giraffe_jobs))

        # duplicate marking
        lc_job, dedup_job = self.build_dedup_job(
            lc_score, deduped_bam, aligned_bams
        )
        dag.add_job(lc_job, set(surject_jobs))
        dag.add_job(dedup_job, {lc_job})

        # indel realignment
        realign_job, rehead_job = self.build_realign_job(
            self.realigned_cram, deduped_bam
        )
        dag.add_job(realign_job, {dedup_job})
        dag.add_job(rehead_job, {realign_job})

        # multiqc
        multiqc_job = self.multiqc()
        if multiqc_job:
            dag.add_job(multiqc_job, {rehead_job})

        # small variant calling
        dnascope_job, dnamodelapply_job = self.build_smallvar_job(
            dnascope_tmp, self.realigned_cram
        )
        dag.add_job(dnascope_job, {realign_job})
        dag.add_job(dnamodelapply_job, {dnascope_job})

        # estimate ploidy
        ploidy_job = self.build_ploidy_job(self.ploidy_json, deduped_bam)
        dag.add_job(ploidy_job, {dedup_job})

        # HLA and KIR calling with T1K
        hla_job, kir_job = self.build_diverse_jobs(out_hla, out_kir)
        dag.add_job(hla_job)
        dag.add_job(kir_job)

        # special-caller
        segdup_job = self.build_segdup_job(out_segdup, self.realigned_cram)
        dag.add_job(segdup_job, {realign_job})

        return dag

    def create_alignment_headers(self) -> List[pathlib.Path]:
        """Create header files for the surjected alignments"""
        aln_headers: List[pathlib.Path] = []

        contigs: List[Tuple[str, str]] = []
        fai_file = str(self.reference) + ".fai"
        if not self.dry_run:
            with open(fai_file) as fh:
                for line in fh:
                    line_split = line.rstrip().split("\t")
                    contigs.append((line_split[0], line_split[1]))

        for i, rg in enumerate(self.readgroups):
            rg_dict = parse_rg_line(rg.replace(r"\t", "\t"))

            aln_header = self.tmp_dir.joinpath(f"sample_aligned_{i}.hdr")
            if not self.dry_run:
                with open(aln_header, "w") as fh:
                    print("@HD\tVN:1.5", file=fh)
                    print(
                        "@RG\t"
                        + "\t".join([x + ":" + y for x, y in rg_dict.items()])
                    )
                    for contig, contig_l in contigs:
                        print(f"@SQ\tSN:{contig}\tLN:{contig_l}")
            aln_headers.append(aln_header)

        return aln_headers

    def build_kmc_job(self, kmer_prefix: pathlib.Path) -> Job:
        """Build KMC k-mer counting jobs"""
        # Create file list for KMC
        file_list = pathlib.Path(str(kmer_prefix) + ".paths")
        all_fastqs = []

        # Add R1 files
        all_fastqs.extend(self.r1_fastq)

        # Add R2 files if present
        if self.r2_fastq:
            all_fastqs.extend(self.r2_fastq)

        # Write file list
        if not self.dry_run:
            with open(file_list, "w") as f:
                for fq in all_fastqs:
                    f.write(f"{fq}\n")

        # Create KMC job
        kmc_job = Job(
            cmds.cmd_kmc(
                kmer_prefix,
                file_list,
                self.tmp_dir,
                memory=self.kmer_memory,
                threads=self.cores,
            ),
            "kmc",
            self.cores,
        )

        return kmc_job

    def build_haplotypes_job(
        self, output_gbz: pathlib.Path, kmer_file: pathlib.Path
    ) -> Job:
        """Build vg haplotypes job"""
        assert self.hapl
        assert self.gbz

        haplotypes_job = Job(
            cmds.cmd_vg_haplotypes(
                output_gbz,
                kmer_file,
                self.hapl,
                self.gbz,
                threads=self.cores,
            ),
            "vg-haplotypes",
            self.cores,
        )
        return haplotypes_job

    def build_alignment_jobs(
        self, sample_pangenome: pathlib.Path
    ) -> Tuple[List[pathlib.Path], List[Job]]:
        """Build vg giraffe alignment jobs"""
        # Process each pair of fastq files
        sample_gams: List[pathlib.Path] = []
        giraffe_jobs: List[Job] = []
        for i, (r1, r2, rg) in enumerate(
            itertools.zip_longest(
                self.r1_fastq, self.r2_fastq, self.readgroups
            )
        ):
            sample_gam = self.tmp_dir.joinpath(f"sample_aligned_{i}.gam")
            rg_dict = parse_rg_line(rg.replace(r"\t", "\t"))
            giraffe_job = Job(
                cmds.cmd_vg_giraffe(
                    sample_gam,
                    sample_pangenome,
                    fastq1=r1,
                    fastq2=r2,
                    readgroup=rg_dict["ID"],
                    sample=rg_dict["SM"],
                    threads=self.cores,
                ),
                "vg-giraffe",
                self.cores,
            )

            sample_gams.append(sample_gam)
            giraffe_jobs.append(giraffe_job)
        return (sample_gams, giraffe_jobs)

    def build_pack_job(
        self,
        sample_pack: pathlib.Path,
        sample_gams: List[pathlib.Path],
    ) -> Job:
        """Compute read support with vg pack"""
        assert self.gbz

        pack_job = Job(
            cmds.cmd_vg_pack(
                sample_pack,
                sample_gams,
                self.gbz,
            ),
            "vg-pack",
            self.cores,
        )
        return pack_job

    def build_call_job(
        self,
        out_svs: pathlib.Path,
        sample_pack: pathlib.Path,
        sample_name="",
    ) -> Job:
        """Call SVs with vg call"""
        assert self.gbz
        assert self.snarls

        call_job = Job(
            cmds.cmd_sv_call(
                out_svs,
                sample_pack,
                self.gbz,
                self.snarls,
                sample_name=sample_name,
            ),
            "vg-call",
            self.cores,
        )
        return call_job

    def build_surject_jobs(
        self,
        sample_gams: List[pathlib.Path],
        alignment_headers: List[pathlib.Path],
    ) -> Tuple[List[pathlib.Path], List[Job]]:
        """Surject gam files on the linear reference"""
        assert self.xg

        aligned_bams: List[pathlib.Path] = []
        surject_jobs: List[Job] = []

        for i, (sample_gam, sample_hdr) in enumerate(
            zip(sample_gams, alignment_headers)
        ):
            sample_bam = self.tmp_dir.joinpath(f"sample_aligned_{i}.bam")
            surject_job = Job(
                cmds.cmd_vg_surject(
                    sample_bam,
                    sample_hdr,
                    sample_gam,
                    self.xg,
                    threads=self.cores,
                ),
                "vg_surject",
                self.cores,
            )

            aligned_bams.append(sample_bam)
            surject_jobs.append(surject_job)
        return (aligned_bams, surject_jobs)

    def build_dedup_job(
        self,
        lc_score: pathlib.Path,
        deduped_bam: pathlib.Path,
        aligned_bams: List[pathlib.Path],
    ) -> Tuple[Job, Job]:
        """Mark duplicates"""
        assert self.output_vcf

        # Output metrics
        metric_base = self.output_vcf.name.replace(".vcf.gz", "") + ".txt"
        metrics_dir = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_metrics")
        )
        mqbc_metrics = metrics_dir.joinpath(
            metric_base + ".mean_qual_by_cycle.txt"
        )
        bdbc_metrics = metrics_dir.joinpath(
            metric_base + ".base_distribution_by_cycle.txt"
        )
        qualdist_metrics = metrics_dir.joinpath(
            metric_base + ".qual_distribution.txt"
        )
        as_metrics = metrics_dir.joinpath(metric_base + ".alignment_stat.txt")
        gc_metrics = metrics_dir.joinpath(metric_base + ".gc_bias.txt")
        gc_summary = metrics_dir.joinpath(metric_base + ".gc_bias_summary.txt")
        dedup_metrics = metrics_dir.joinpath(
            metric_base + ".dedup_metrics.txt"
        )

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=aligned_bams,
        )
        driver.add_algo(LocusCollector(lc_score))
        driver.add_algo(MeanQualityByCycle(mqbc_metrics))
        driver.add_algo(BaseDistributionByCycle(bdbc_metrics))
        driver.add_algo(QualDistribution(qualdist_metrics))
        driver.add_algo(AlignmentStat(as_metrics))
        driver.add_algo(GCBias(gc_metrics, summary=gc_summary))

        lc_job = Job(
            shlex.join(driver.build_cmd()),
            "locuscollector",
            self.cores,
        )

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=aligned_bams,
        )
        driver.add_algo(
            Dedup(
                deduped_bam,
                lc_score,
                metrics=dedup_metrics,
            )
        )
        dedup_job = Job(shlex.join(driver.build_cmd()), "dedup", self.cores)

        return (lc_job, dedup_job)

    def build_realign_job(
        self,
        realigned_cram: pathlib.Path,
        deduped_bam: pathlib.Path,
    ) -> Tuple[Job, Job]:
        assert self.output_vcf

        # Output metrics
        metric_base = self.output_vcf.name.replace(".vcf.gz", "") + ".txt"
        metrics_dir = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_metrics")
        )

        # Output metrics
        # - Run some metrics after duplicate marking
        is_metrics = metrics_dir.joinpath(metric_base + ".insert_size.txt")
        wgs_metrics = metrics_dir.joinpath(metric_base + ".wgs.txt")
        wgs_metrics_tmp = metrics_dir.joinpath(metric_base + ".wgs.txt.tmp")

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=[deduped_bam],
        )
        driver.add_algo(InsertSizeMetricAlgo(is_metrics))
        driver.add_algo(WgsMetricsAlgo(wgs_metrics, include_unpaired="true"))
        driver.add_algo(
            Realigner(
                realigned_cram,
                known_sites=self.known_sites,
            )
        )
        realigner_job = Job(
            shlex.join(driver.build_cmd()),
            "realigner",
            self.cores,
        )

        # Rehead WGS metrics so they are recognized by MultiQC
        rehead_job = Job(
            cmds.rehead_wgsmetrics(wgs_metrics, wgs_metrics_tmp),
            "Rehead metrics",
            0,
        )
        return (realigner_job, rehead_job)

    def multiqc(self) -> Optional[Job]:
        """Run MultiQC on the metrics files"""
        assert self.output_vcf
        if not self.skip_version_check:
            if not all(
                [
                    check_version(cmd, min_version)
                    for (cmd, min_version) in MULTIQC_MIN_VERSION.items()
                ]
            ):
                self.logger.warning(
                    "Skipping MultiQC. MultiQC version %s or later not found",
                    MULTIQC_MIN_VERSION["multiqc"],
                )
                return None

        metrics_dir = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_metrics")
        )
        multiqc_job = Job(
            cmds.cmd_multiqc(
                metrics_dir,
                metrics_dir,
                f"Generated by the Sentieon-CLI version {__version__}",
            ),
            "multiqc",
            0,
        )
        return multiqc_job

    def build_smallvar_job(
        self,
        dnascope_tmp: pathlib.Path,
        realigned_cram: pathlib.Path,
    ) -> Tuple[Job, Job]:
        """Small variant calling"""
        assert self.model_bundle
        assert self.output_vcf

        model = self.model_bundle.joinpath("dnascope.model")
        interval_ctgs = list(["chr" + str(x) for x in range(1, 23)])
        interval_ctgs.extend(["chrX", "chrY"])
        interval = ",".join(interval_ctgs)

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=[realigned_cram],
            interval=interval,
        )
        driver.add_algo(
            DNAscope(
                dnascope_tmp,
                dbsnp=self.dbsnp,
                pcr_indel_model="NONE",
                model=model,
            )
        )
        dnascope_job = Job(
            shlex.join(driver.build_cmd()), "dnascope", self.cores
        )

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            DNAModelApply(
                model,
                dnascope_tmp,
                self.output_vcf,
            )
        )
        dnamodelapply_job = Job(
            shlex.join(driver.build_cmd()), "model-apply", self.cores
        )

        return (dnascope_job, dnamodelapply_job)

    def build_ploidy_job(
        self,
        ploidy_json: pathlib.Path,
        deduped_bam: pathlib.Path,
    ) -> Job:
        """Estimate sample ploidy and sex"""
        estimate_ploidy = pathlib.Path(
            str(files("sentieon_cli.scripts").joinpath("estimate_ploidy.py"))
        ).resolve()
        ploidy_job = Job(
            cmds.cmd_estimate_ploidy(
                ploidy_json,
                deduped_bam,
                estimate_ploidy,
            ),
            "estimate-ploidy",
            0,
        )
        return ploidy_job

    def get_sex(self) -> None:
        """Retrieve the sample sex"""
        with open(self.ploidy_json) as fh:
            data = json.load(fh)
            sex = data["sex"]
            if sex == "female":
                self.sample_sex = SampleSex.FEMALE
            elif sex == "male":
                self.sample_sex = SampleSex.MALE
            else:
                self.sample_sex = SampleSex.UNKNOWN

    def build_diverse_jobs(
        self,
        out_hla: pathlib.Path,
        out_kir: pathlib.Path,
    ) -> Tuple[Job, Job]:
        """Genotype diverse sequences with T1K"""
        assert self.t1k_hla
        assert self.t1k_kir

        hla_job = Job(
            cmds.cmd_t1k(
                out_hla,
                self.r1_fastq,
                self.r2_fastq,
                self.t1k_hla,
                preset="hla-wgs",
                threads=self.cores,
            ),
            "T1K-HLA",
            self.cores,
        )

        kir_job = Job(
            cmds.cmd_t1k(
                out_kir,
                self.r1_fastq,
                self.r2_fastq,
                self.t1k_kir,
                preset="kir-wgs",
                threads=self.cores,
            ),
            "T1K-HLA",
            self.cores,
        )

        return (hla_job, kir_job)

    def build_segdup_job(
        self,
        out_segdup: pathlib.Path,
        realigned_cram: pathlib.Path,
    ) -> Job:
        """Call variants in difficult SegDups"""
        assert self.reference
        assert self.model_bundle

        segdup_job = Job(
            cmds.cmd_segdup_caller(
                out_segdup,
                realigned_cram,
                reference=self.reference,
                sr_bundle=self.model_bundle,
                genes="CFH,CFHR3,CYP11B1,CYP2D6,GBA,NCF1,PMS2,SMN1,STRC",
            ),
            "segdup-caller",
            self.cores,
        )
        return segdup_job

    def build_second_dag(self) -> DAG:
        """Build the second DAG for the pangenome pipeline"""
        self.logger.info("Building the second pangenome DAG")
        dag = DAG()

        # Output files
        out_cnv = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_cnv.vcf.gz")
        )
        out_expansions = pathlib.Path(
            str(self.output_vcf).replace(".vcf.gz", "_expansion")
        )

        # Intermediate file paths
        cnvscope_tmp = self.tmp_dir.joinpath("sample_cnvcope_tmp.vcf.gz")
        cnvscope_autosomes = self.tmp_dir.joinpath(
            "sample_cnvcope_autosomes.vcf.gz"
        )
        cnvscope_haploid_tmp = self.tmp_dir.joinpath(
            "sample_cnvcope_haploid_tmp.vcf.gz"
        )
        cnvscope_haploid = self.tmp_dir.joinpath(
            "sample_cnvcope_haploid.vcf.gz"
        )

        # CNV calling
        (
            cnvscope_job,
            cnvapply_job,
            haploid_cnv_job,
            haploid_apply_job,
            concat_job,
        ) = self.build_cnv_jobs(
            out_cnv,
            cnvscope_tmp,
            cnvscope_autosomes,
            cnvscope_haploid_tmp,
            cnvscope_haploid,
            self.realigned_cram,
        )
        dag.add_job(cnvscope_job)
        dag.add_job(cnvapply_job, {cnvscope_job})
        if haploid_cnv_job and haploid_apply_job and concat_job:
            dag.add_job(haploid_cnv_job)
            dag.add_job(haploid_apply_job, {haploid_cnv_job})
            dag.add_job(concat_job, {cnvapply_job, haploid_apply_job})

        # Repeat expansions
        expansion_job = self.build_expansion_job(
            out_expansions, self.realigned_cram
        )
        dag.add_job(expansion_job)

        return dag

    def build_cnv_jobs(
        self,
        out_cnv: pathlib.Path,
        cnvscope_tmp: pathlib.Path,
        cnvscope_autosomes: pathlib.Path,
        cnvscope_haploid_tmp: pathlib.Path,
        cnvscope_haploid: pathlib.Path,
        realigned_cram: pathlib.Path,
    ) -> Tuple[Job, Job, Optional[Job], Optional[Job], Optional[Job]]:
        """Call CNVs with CNVscope"""
        assert self.model_bundle

        model = self.model_bundle.joinpath("cnv.model")

        # First pass will be autosomes for males
        # or autosomes + chrX for females
        first_contigs = list(["chr" + str(x) for x in range(1, 23)])
        if self.sample_sex != SampleSex.MALE:
            first_contigs.append("chrX")

        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=[realigned_cram],
            interval=first_contigs,
        )
        driver.add_algo(
            CNVscope(
                cnvscope_tmp,
                model=model,
            )
        )
        cnvscope_job = Job(
            shlex.join(driver.build_cmd()), "cnvscope-autosomes", self.cores
        )

        # Autosome model apply
        outfile = (
            cnvscope_autosomes
            if self.sample_sex == SampleSex.MALE
            else out_cnv
        )
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            CNVModelApply(
                outfile,
                model,
                cnvscope_tmp,
            )
        )
        cnvapply_job = Job(
            shlex.join(driver.build_cmd()),
            "cnvmodelapply-autosomes",
            self.cores,
        )

        if self.sample_sex != SampleSex.MALE:
            return (cnvscope_job, cnvapply_job, None, None, None)

        # Sex chrom CNVscope
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
            input=[realigned_cram],
            interval="chrX,chrY",
        )
        driver.add_algo(
            CNVscope(
                cnvscope_haploid_tmp,
                model=model,
            )
        )
        haploid_cnv_job = Job(
            shlex.join(driver.build_cmd()), "cnvscope-haploid", self.cores
        )

        # Sex chrom CNVModelApply
        driver = Driver(
            reference=self.reference,
            thread_count=self.cores,
        )
        driver.add_algo(
            CNVModelApply(
                cnvscope_haploid,
                model,
                cnvscope_haploid_tmp,
            )
        )
        haploid_apply_job = Job(
            shlex.join(driver.build_cmd()), "cnvmodelapply-haploid", self.cores
        )

        # Concate autosome and sex chromsome VCFs
        concat_job = Job(
            cmds.bcftools_concat(
                out_cnv, [cnvscope_autosomes, cnvscope_haploid]
            ),
            "concat-calls",
            0,
        )

        return (
            cnvscope_job,
            cnvapply_job,
            haploid_cnv_job,
            haploid_apply_job,
            concat_job,
        )

    def build_expansion_job(
        self,
        out_expansions: pathlib.Path,
        realigned_cram: pathlib.Path,
    ) -> Job:
        """Identify repeat expansions"""
        assert self.reference
        assert self.expansion_catalog

        expansion_job = Job(
            cmds.cmd_expansion_hunter(
                out_expansions,
                realigned_cram,
                reference=self.reference,
                variant_catalog=self.expansion_catalog,
                sex="male" if self.sample_sex == SampleSex.MALE else "female",
                threads=self.cores,
            ),
            "expansion-hunter",
            0,
        )
        return expansion_job
