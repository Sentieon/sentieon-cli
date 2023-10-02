"""
This module contains the functions that accept **kwargs and return the
command strings.

We expect that the kwargs are created from user-input or a config file.

Command strings may be partial--that is, they could be mixed with other
command strings to get a full viable command.

"""


def cmd_algo_dnascope(**kwargs):
    """
    sentieon driver -t "$_arg_threads" -r "$_arg_reference_fasta" \
        ${_arg_regions:+--interval "$_arg_regions"} -i "$_arg_sample_input" \
        ${_arg_gvcf:+--algo DNAscope --model "$_arg_model_bundle"/gvcf_model
        --emit_mode gvcf "$DIPLOID_GVCF"} \
        --algo DNAscope ${_arg_dbsnp:+--dbsnp "$_arg_dbsnp"} \
        --model "$_arg_model_bundle"/diploid_model "$DIPLOID_TMP_OUT"
    """
    # TODO: this is used elsewhere, should create once and pass.
    diploid_tmp_out = f"{kwargs['tmp_base']}_diploid_tmp.vcf.gz"
    # https://github.com/Sentieon/sentieon-scripts/blob/8d33f29e442d5a1e782445f06bc1f11e278d8f87/dnascope_LongRead/dnascope_HiFi.sh#L355-L359
    if kwargs.get("gvcf"):
        diploid_gvcf = f"{kwargs['tmp_base']}_diploid.g.vcf.gz"
        gvcf = f"--algo DNAscope --model {kwargs['model-bundle']}/gvcf_model"
        gvcf += f" --emit_mode gvcf {gvcf} {diploid_gvcf}"
    else:
        gvcf = ""
    if kwargs.get("dbsnp"):
        dbsnp = f"--dbsnp {kwargs['dbsnp'].name}"
    else:
        dbsnp = ""
    cmd = f"{gvcf} --algo DNAscope {dbsnp} --model "
    cmd += f"{kwargs['model-bundle']}/diploid_model {diploid_tmp_out}"
    return cmd
