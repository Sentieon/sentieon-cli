"""
Classes for interacting with the Sentieon driver
"""

import itertools
import pathlib
from typing import Optional, List, Union

from .logging import get_logger

logger = get_logger(__name__)


class BaseAlgo:
    """A base class for Sentieon algos"""

    name = "BaseAlgo"

    def build_cmd(self) -> List[str]:
        """Build a command line for the algo"""
        cmd: List[str] = ["--algo", self.name]

        for k, v in self.__dict__.items():
            if k == "output":
                continue
            elif v is None:
                continue
            elif isinstance(v, list):
                for i in v:
                    cmd.append(f"--{k}")
                    cmd.append(str(i))
            elif isinstance(v, bool):
                if v:
                    cmd.append(f"--{k}")
            else:
                cmd.append(f"--{k}")
                cmd.append(str(v))

        if "output" in self.__dict__:
            cmd.append(str(self.__dict__["output"]))

        return cmd


class VariantPhaser(BaseAlgo):
    """algo VariantPhaser"""

    name = "VariantPhaser"

    def __init__(
        self,
        vcf: pathlib.Path,
        output: pathlib.Path,
        out_bed: Optional[pathlib.Path] = None,
        out_ext: Optional[pathlib.Path] = None,
        max_depth: int = 1000,
    ):
        self.vcf = vcf
        self.output = output
        self.out_bed = out_bed
        self.out_ext = out_ext
        self.max_depth = max_depth


class RepeatModel(BaseAlgo):
    """algo RepeatModel"""

    name = "RepeatModel"

    def __init__(
        self,
        output: pathlib.Path,
        phased: bool = False,
        min_map_qual: int = 1,
        min_group_count: int = 10000,
        read_flag_mask: Optional[str] = None,
        repeat_extension: int = 5,
        max_repeat_unit_size: int = 2,
        min_repeat_count: int = 6,
    ):
        self.output = output
        self.phased = phased
        self.min_map_qual = min_map_qual
        self.min_group_count = min_group_count
        self.read_flag_mask = read_flag_mask
        self.repeat_extension = repeat_extension
        self.max_repeat_unit_size = max_repeat_unit_size
        self.min_repeat_count = min_repeat_count


class DNAModelApply(BaseAlgo):
    """algo DNAModelApply"""

    name = "DNAModelApply"

    def __init__(
        self,
        model: pathlib.Path,
        vcf: pathlib.Path,
        output: pathlib.Path,
    ):
        self.model = model
        self.vcf = vcf
        self.output = output


class DNAscope(BaseAlgo):
    """algo DNAscope"""

    name = "DNAscope"

    def __init__(
        self,
        output: pathlib.Path,
        dbsnp: Optional[pathlib.Path] = None,
        emit_mode: Optional[str] = None,
        model: Optional[pathlib.Path] = None,
        pcr_indel_model: Optional[str] = None,
        var_type: Optional[str] = None,
        ploidy: Optional[int] = None,
    ):
        self.output = output
        self.dbsnp = dbsnp
        self.emit_mode = emit_mode
        self.model = model
        self.pcr_indel_model = pcr_indel_model
        self.var_type = var_type
        self.ploidy = ploidy


class DNAscopeHP(BaseAlgo):
    """algo DNAscopeHP"""

    name = "DNAscopeHP"

    def __init__(
        self,
        output: pathlib.Path,
        dbsnp: Optional[pathlib.Path] = None,
        model: Optional[pathlib.Path] = None,
        pcr_indel_model: Optional[pathlib.Path] = None,
    ):
        self.output = output
        self.dbsnp = dbsnp
        self.model = model
        self.pcr_indel_model = pcr_indel_model


class LongReadSV(BaseAlgo):
    """algo LongReadSV"""

    name = "LongReadSV"

    def __init__(
        self,
        output: pathlib.Path,
        model: Optional[pathlib.Path] = None,
        min_map_qual: Optional[int] = None,
        min_sv_size: Optional[int] = None,
        min_dp: Optional[int] = None,
        min_af: Optional[float] = None,
    ):
        self.output = output
        self.model = model
        self.min_map_qual = min_map_qual
        self.min_sv_size = min_sv_size
        self.min_dp = min_dp
        self.min_af = min_af


class LocusCollector(BaseAlgo):
    """algo LocusCollector"""

    name = "LocusCollector"

    def __init__(
        self,
        output: pathlib.Path,
        rna: bool = False,
        consensus: bool = False,
        umi_tag: Optional[str] = None,
        umi_ecc_dist: Optional[int] = None,
        umi_ecc_lev_dist: Optional[int] = None,
    ):
        self.output = output
        self.rna = rna
        self.consensus = consensus
        self.umi_tag = umi_tag
        self.umi_ecc_dist = umi_ecc_dist
        self.umi_ecc_lev_dist = umi_ecc_lev_dist


class Dedup(BaseAlgo):
    """algo Dedup"""

    name = "Dedup"

    def __init__(
        self,
        output: pathlib.Path,
        score_info: pathlib.Path,
        bam_compression: Optional[int] = None,
        cram_write_options: Optional[str] = None,
        metrics: Optional[pathlib.Path] = None,
        rmdup: bool = False,
    ):
        self.output = output
        self.score_info = score_info
        self.bam_compression = bam_compression
        self.cram_write_options = cram_write_options
        self.metrics = metrics
        self.rmdup = rmdup


class GVCFtyper(BaseAlgo):
    """algo GVCFtyper"""

    name = "GVCFtyper"

    def __init__(
        self,
        output: pathlib.Path,
        vcf: pathlib.Path,
    ):
        self.output = output
        self.vcf = vcf


class SVSolver(BaseAlgo):
    """algo SVSolver"""

    name = "SVSolver"

    def __init__(
        self,
        output: pathlib.Path,
        vcf: pathlib.Path,
    ):
        self.output = output
        self.vcf = vcf


class ReadWriter(BaseAlgo):
    """algo ReadWriter"""

    name = "ReadWriter"

    def __init__(
        self,
        output: pathlib.Path,
        bam_compression: Optional[int] = None,
        cram_write_options: Optional[str] = None,
        output_flag_filter: Optional[str] = None,
    ):
        self.output = output
        bam_compression = bam_compression
        cram_write_options = cram_write_options
        output_flag_filter = output_flag_filter


class InsertSizeMetricAlgo(BaseAlgo):
    """algo InsertSizeMetricAlgo"""

    name = "InsertSizeMetricAlgo"

    def __init__(self, output: pathlib.Path):
        self.output = output


class MeanQualityByCycle(BaseAlgo):
    """algo MeanQualityByCycle"""

    name = "MeanQualityByCycle"

    def __init__(self, output: pathlib.Path):
        self.output = output


class BaseDistributionByCycle(BaseAlgo):
    """algo BaseDistributionByCycle"""

    name = "BaseDistributionByCycle"

    def __init__(self, output: pathlib.Path):
        self.output = output


class QualDistribution(BaseAlgo):
    """algo QualDistribution"""

    name = "QualDistribution"

    def __init__(self, output: pathlib.Path):
        self.output = output


class GCBias(BaseAlgo):
    """algo GCBias"""

    name = "GCBias"

    def __init__(
        self,
        output: pathlib.Path,
        summary: Optional[pathlib.Path] = None,
    ):
        self.output = output
        self.summary = summary


class AlignmentStat(BaseAlgo):
    """algo AlignmentStat"""

    name = "AlignmentStat"

    def __init__(
        self,
        output: pathlib.Path,
        adapter_seq: str = "",
    ):
        self.output = output
        self.adapter_seq = adapter_seq


class CoverageMetrics(BaseAlgo):
    """algo CoverageMetrics"""

    name = "CoverageMetrics"

    def __init__(
        self,
        output: pathlib.Path,
        omit_base_output: bool = True,
    ):
        self.output = output
        self.omit_base_output = omit_base_output


class HsMetricAlgo(BaseAlgo):
    """algo HsMetricAlgo"""

    name = "HsMetricAlgo"

    def __init__(
        self,
        output: pathlib.Path,
        targets_list: pathlib.Path,
        baits_list: pathlib.Path,
    ):
        self.output = output
        self.targets_list = targets_list
        self.baits_list = baits_list


class SequenceArtifactMetricsAlgo(BaseAlgo):
    """algo SequenceArtifactMetricsAlgo"""

    name = "SequenceArtifactMetricsAlgo"

    def __init__(
        self,
        output: pathlib.Path,
        dbsnp: Optional[pathlib.Path] = None,
    ):
        self.output = output
        self.dbsnp = dbsnp


class WgsMetricsAlgo(BaseAlgo):
    """algo WgsMetricsAlgo"""

    name = "WgsMetricsAlgo"

    def __init__(
        self,
        output: pathlib.Path,
        include_unpaired: Optional[str] = None,
    ):
        self.output = output
        self.include_unpaired = include_unpaired


class HybridStage1(BaseAlgo):
    """algo HybridStage1"""

    name = "HybridStage1"

    def __init__(
        self,
        output: Union[pathlib.Path, str],
        model: Optional[pathlib.Path] = None,
        fa_file: Optional[pathlib.Path] = None,
        bed_file: Optional[pathlib.Path] = None,
        cut_indel: Optional[int] = None,
        hap_bam: Optional[pathlib.Path] = None,
        hap_bed: Optional[pathlib.Path] = None,
        cut_len: Optional[int] = None,
        split_size: Optional[int] = None,
        split_overlap: Optional[int] = None,
        hap_vcf: Optional[pathlib.Path] = None,
        kmer_missing_rate: Optional[float] = None,
        mk_avg_bq: Optional[int] = None,
        prune_factor: Optional[int] = None,
        mk_max_init_hapcnt: Optional[int] = None,
        ploidy: Optional[int] = None,
    ):
        self.output = output
        self.model = model
        self.fa_file = fa_file
        self.bed_file = bed_file
        self.cut_indel = cut_indel
        self.hap_bam = hap_bam
        self.hap_bed = hap_bed
        self.cut_len = cut_len
        self.split_size = split_size
        self.split_overlap = split_overlap
        self.hap_vcf = hap_vcf
        self.kmer_missing_rate = kmer_missing_rate
        self.mk_avg_bq = mk_avg_bq
        self.prune_factor = prune_factor
        self.mk_max_init_hapcnt = mk_max_init_hapcnt
        self.ploidy = ploidy


class HybridStage2(BaseAlgo):
    """algo HybridStage2"""

    name = "HybridStage2"

    def __init__(
        self,
        all_bed: pathlib.Path,
        out_bed: Optional[pathlib.Path] = None,
        model: Optional[pathlib.Path] = None,
        mode: Optional[str] = None,
        debug: Optional[str] = None,
        hap_bed: Optional[pathlib.Path] = None,
        unmap_bam: Optional[pathlib.Path] = None,
        alt_bam: Optional[pathlib.Path] = None,
        min_kmer_size: Optional[int] = None,
        prune_factor: Optional[int] = None,
        target_extension: Optional[int] = None,
    ):
        self.out_bed = out_bed
        self.model = model
        self.mode = mode
        self.debug = debug
        self.all_bed = all_bed
        self.hap_bed = hap_bed
        self.unmap_bam = unmap_bam
        self.alt_bam = alt_bam
        self.min_kmer_size = min_kmer_size
        self.prune_factor = prune_factor
        self.target_extension = target_extension


class HybridStage3(BaseAlgo):
    """algo HybridStage3"""

    name = "HybridStage3"

    def __init__(
        self,
        output: Union[pathlib.Path, str],
        model: Optional[pathlib.Path] = None,
        region_ext: Optional[int] = 1000,
    ):
        self.output = output
        self.model = model
        region_ext = region_ext


class CNVscope(BaseAlgo):
    """algo CNVscope"""

    name = "CNVscope"

    def __init__(
        self,
        output: Union[pathlib.Path, str],
        model: pathlib.Path,
    ):
        self.output = output
        self.model = model


class CNVModelApply(BaseAlgo):
    """algo CNVModelApply"""

    name = "CNVModelApply"

    def __init__(
        self,
        output: Union[pathlib.Path, str],
        model: pathlib.Path,
        vcf: pathlib.Path,
    ):
        self.output = output
        self.model = model
        self.vcf = vcf


class BaseDriver:
    """A base class for the Sentieon driver"""

    def __init__(
        self,
        reference: Optional[pathlib.Path] = None,
        thread_count: Optional[int] = None,
        interval: Optional[Union[pathlib.Path, str]] = None,
        interval_padding: int = 0,
        read_filter: Optional[List[str]] = None,
        replace_rg: Optional[List[List[str]]] = None,
        input: Optional[List[pathlib.Path]] = None,
        algo: Optional[List[BaseAlgo]] = None,
    ):
        if (
            isinstance(replace_rg, list)
            and isinstance(input, list)
            and len(replace_rg) > len(input)
        ):
            logger.error(
                "More readgroups '%s' than input files '%s' in --replace_rg",
                replace_rg,
                input,
            )
            raise RuntimeError("Error in command")

        self._name = "Unset"
        self.reference = reference
        self.replace_rg = replace_rg or []
        self.input = input or []
        self.thread_count = thread_count
        self.interval = interval
        self.interval_padding = interval_padding
        self.read_filter = read_filter
        self.algo = algo if algo else []

    def add_algo(self, algo: BaseAlgo):
        """Add an algo to the driver"""
        self.algo.append(algo)

    def build_cmd(self) -> List[str]:
        """Build a command line for the driver"""
        cmd: List[str] = ["sentieon", self._name]

        if self.input:
            for replace_rg, this_input in itertools.zip_longest(  # type: ignore  # noqa: E501
                self.replace_rg, self.input, fillvalue=[]
            ):
                for rg_arg in replace_rg:
                    cmd.append("--replace_rg")
                    cmd.append(rg_arg)
                cmd.append("--input")
                cmd.append(str(this_input))

        for k, v in self.__dict__.items():
            if k == "algo":
                continue
            elif k.startswith("_"):
                continue
            elif k == "input":
                continue
            elif k == "replace_rg":
                continue
            elif v is None:
                continue
            elif isinstance(v, list):
                for i in v:
                    cmd.append(f"--{k}")
                    cmd.append(str(i))
            elif isinstance(v, bool):
                if v:
                    cmd.append(f"--{k}")
            else:
                cmd.append(f"--{k}")
                cmd.append(str(v))

        for algo in self.algo:
            cmd.extend(algo.build_cmd())

        return cmd


class Driver(BaseDriver):
    """Representing the Sentieon driver"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._name = "driver"
