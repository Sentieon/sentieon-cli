"""
This module contains the functions that accept kwargs and return the
command strings.

We expect that the kwargs are created from user-input (to argh) or a config
file. Where appropriate, we use named args, but for flexibility, kwargs are
also used.

Command strings may be partial--that is, they could be mixed with other
command strings to get a full viable command.

"""

import io
import pathlib
import shlex
import typing
from typing import Any, Optional, List, Union, Dict
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
        emit_mode: str = "variant",
        model: Optional[pathlib.Path] = None,
    ):
        self.output = output
        self.dbsnp = dbsnp
        self.emit_mode = emit_mode
        self.model = model


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


class Driver:
    """Representing the Sentieon driver"""

    def __init__(
        self,
        reference: Optional[pathlib.Path] = None,
        thread_count: Optional[int] = None,
        interval: Optional[Union[pathlib.Path, str]] = None,
        read_filter: Optional[str] = None,
        input: Optional[List[pathlib.Path]] = None,
        algo: Optional[List[BaseAlgo]] = None,
    ):
        self.reference = reference
        self.input = input
        self.thread_count = thread_count
        self.interval = interval
        self.read_filter = read_filter
        self.algo = algo if algo else []

    def add_algo(self, algo: BaseAlgo):
        """Add an algo to the driver"""
        self.algo.append(algo)

    def build_cmd(self) -> List[str]:
        """Build a command line for the driver"""
        cmd: List[str] = ["sentieon", "driver"]

        for k, v in self.__dict__.items():
            if k == "algo":
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


def cmd_bedtools_subtract(
    regions_bed: typing.Optional[pathlib.Path],
    phased_bed: pathlib.Path,
    unphased_bed: pathlib.Path,
    tmp_dir: pathlib.Path,
    reference: pathlib.Path,
    dry_run: bool,
):
    if regions_bed is None:
        # set region to the full genome
        regions_bed = tmp_dir.joinpath("reference.bed")
        if not dry_run:
            with open(regions_bed, "wt", encoding="utf-8") as f:
                for line in open(name(reference) + ".fai", encoding="utf-8"):
                    toks = line.strip().split("\t")
                    f.write(f"{toks[0]}\t0\t{toks[1]}\n")
    cmd = [
        "bedtools",
        "subtract",
        "-a",
        str(regions_bed),
        "-b",
        str(phased_bed),
    ]
    return shlex.join(cmd) + ">" + shlex.quote(str(unphased_bed))


def name(path: typing.Union[str, io.TextIOWrapper, pathlib.Path]) -> str:
    """Return the name of a file that may also be a text-wrapper."""
    if isinstance(path, io.TextIOWrapper):
        return path.name
    elif isinstance(path, pathlib.Path):
        return str(path)
    return path


def cmd_pyexec_vcf_mod_haploid_patch(
    hap1_patch: str,
    hap2_patch: str,
    hap_patt: str,
    tech: str,
    phased_vcf: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> str:
    """
    merge dnascope and dnascope-hp variants

    >>> assert '--hap1 ' in cmd_pyexec_vcf_mod_haploid_patch("h1.vcf.gz",
    ... "h2.vcf.gz", "out_hap%d_%stmp.vcf.gz", "HiFi", "ph.vcf", 2,
    ... {'vcf_mod_py': 'vcf_mod.py'})
    >>> assert not "--hap1 " in cmd_pyexec_vcf_mod_haploid_patch("h1.vcf.gz",
    ... "h2.vcf.gz", "out_hap%d_%stmp.vcf.gz", "ONT", "ph.vcf", 2,
    ... {'vcf_mod_py': 'vcf_mod.py'})

    """
    assert tech.upper() in ("HIFI", "ONT")

    cmd = f"sentieon pyexec {kwargs['vcf_mod_py']} -t {cores} "
    cmd += "haploid_patch "
    cmd += f"--patch1 {hap1_patch} --patch2 {hap2_patch}"
    if tech.upper() == "HIFI":
        cmd += " --hap1 " + hap_patt % (1, "nohp_")
        cmd += " --hap2 " + hap_patt % (2, "nohp_")
    else:
        cmd += f" --phased {phased_vcf}"
    cmd += " --hap1_hp " + hap_patt % (1, "")
    cmd += " --hap2_hp " + hap_patt % (2, "")
    return cmd


def cmd_pyexec_vcf_mod_patch(
    out_vcf: str,
    vcf: str,
    vcf_hp: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> str:
    """Patch DNAscope and DNAscopeHP VCF files"""

    cmd = [
        "sentieon",
        "pyexec",
        str(kwargs["vcf_mod_py"]),
        "-t",
        str(cores),
        "patch",
        "--vcf",
        str(vcf),
        "--vcf_hp",
        str(vcf_hp),
        str(out_vcf),
    ]
    return shlex.join(cmd)


def cmd_pyexec_gvcf_combine(
    reference: pathlib.Path,
    gvcf: str,
    out_vcf: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> str:
    """Combine gVCF files"""

    cmd1 = [
        "sentieon",
        "pyexec",
        str(kwargs["gvcf_combine_py"]),
        "-t",
        str(cores),
        str(reference),
        gvcf,
        out_vcf,
        "-",
    ]
    cmd2 = [
        "sentieon",
        "util",
        "vcfconvert",
        "-",
        out_vcf.replace(".vcf.gz", ".g.vcf.gz"),
    ]
    return shlex.join(cmd1) + "|" + shlex.join(cmd2)


def cmd_pyexec_vcf_mod_merge(
    hap1_vcf: str,
    hap2_vcf: str,
    unphased_vcf: str,
    phased_vcf: str,
    phased_bed: str,
    out_vcf: str,
    cores: int,
    kwargs: Dict[str, Any],
) -> str:
    """Merge haploid VCF files"""

    cmd = [
        "sentieon",
        "pyexec",
        kwargs["vcf_mod_py"],
        "-t",
        str(cores),
        "merge",
        "--hap1",
        hap1_vcf,
        "--hap2",
        hap2_vcf,
        "--unphased",
        unphased_vcf,
        "--phased",
        phased_vcf,
        "--bed",
        phased_bed,
        out_vcf,
    ]
    return shlex.join(cmd)
