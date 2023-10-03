"""
This module contains the functions that accept **kwargs and return the
command strings.

We expect that the kwargs are created from user-input (to argh) or a config
file.

Command strings may be partial--that is, they could be mixed with other
command strings to get a full viable command.

"""
import io
import typing
from .logging import get_logger

logger = get_logger(__name__)


def cmd_model_apply(**kwargs):
    """
     sentieon driver -t "$_arg_threads" -r "$_arg_reference_fasta" \
        --algo DNAModelApply --model "$model" -v "$input_vcf" "$output_vcf"
    """
    inp_vcf = f"{kwargs['tmp_base']}/out_diploid_tmp.vcf.gz"
    out_vcf = f"{kwargs['tmp_base']}/diploid.vcf.gz"
    cmd = _cmd_sentieon_driver(**kwargs)
    cmd += f"--algo DNAmodelApply --model {kwargs['model']} "
    cmd += f"-v {inp_vcf} {out_vcf}"
    return cmd


def name(path: typing.Union[str, io.TextIOWrapper]) -> str:
    """Return the name of a file that may also be a text-wrapper."""
    if hasattr(path, "name"):
        path = path.name
    return path


def _cmd_sentieon_driver(**kwargs) -> str:
    """
    Common base for running sentieon driver.
    This is only called along with other commands.
    NOTE: sometimes we need to run without bed file in later command
    whereas an earlier command might have a bed file.
    Same for read-filter.
    """
    logger.debug(kwargs)
    if kwargs.get("bed"):
        bed = f" --interval {name(kwargs['bed'])}"
    else:
        bed = ""
    if kwargs.get("read-filter"):
        read_filter = f" --read_filter {kwargs['read-filter']}"
    else:
        read_filter = ""

    cmd = (
        f"sentieon driver -t {kwargs['cores']} -r {name(kwargs['reference'])} "
    )
    cmd += f"{bed} -i {name(kwargs['sample_input'])}{read_filter}"
    return cmd


def cmd_algo_dnascope(**kwargs):
    """
        ${_arg_gvcf:+--algo DNAscope --model "$_arg_model_bundle"/gvcf_model
        --emit_mode gvcf "$DIPLOID_GVCF"} \
        --algo DNAscope ${_arg_dbsnp:+--dbsnp "$_arg_dbsnp"} \
        --model "$_arg_model_bundle"/diploid_model "$DIPLOID_TMP_OUT"
    """
    # TODO: this is used elsewhere, should create once and pass.
    diploid_tmp_out = f"{kwargs['tmp_base']}/out_diploid_tmp.vcf.gz"
    # https://github.com/Sentieon/sentieon-scripts/blob/8d33f29e442d5a1e782445f06bc1f11e278d8f87/dnascope_LongRead/dnascope_HiFi.sh#L355-L359
    if kwargs.get("gvcf") or True:  # NOTE: this is always done in bash.
        diploid_gvcf = f"{kwargs['tmp_base']}/out_diploid.g.vcf.gz"
        gvcf = f" --algo DNAscope --model {kwargs['model_bundle']}/gvcf_model"
        gvcf += f" --emit_mode gvcf {diploid_gvcf} "
    else:
        gvcf = ""
    if kwargs.get("dbsnp"):
        dbsnp = f" --dbsnp {kwargs['dbsnp'].name} "
    else:
        dbsnp = ""
    cmd = f"{gvcf}--algo DNAscope {dbsnp} --model "
    cmd += f"{kwargs['model_bundle']}/diploid_model {diploid_tmp_out} "
    return _cmd_sentieon_driver(**kwargs) + cmd
