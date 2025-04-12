# DNAscope hybrid

Sentieon DNAscope hybrid is a pipeline for germline variant calling using combined short-read and long-read data from a single sample. The DNAscope hybrid pipeline is able to utilize the strengths of both short and long-read technologies to generate variant callsets that are more accurate than either short-read or long-read data alone.

The pipeline supports input data in the following formats; both short-read and long-read input are required:
* Unaligned short-read data in gzipped FASTQ format.
* Aligned short-reads in BAM or CRAM format.
* Unaliged long-read data in the uBAM or uCRAM format.
* Aligned long-read data in BAM or CRAM format.

By default, the pipeline will generate the following output files:
* Small variants (SNVs and indels) in the VCF format.
* Structural variants in the VCF format.
* Copy-number variants in the VCF format.

If unaligned reads are used as input, the pipeline will also output aligned reads in BAM or CRAM format.

The DNAscope hybrid pipeline is implemented using the Sentieon software package, which requires a valid license for use. Please contact info@sentieon.com for access to the Sentieon software and an evaluation license.

## Prerequisites

- Sentieon software package version 202503 or higher.
- [Python] version 3.8 or higher.
- [bcftools] version 1.10 or higher.
- [bedtools]
- [MultiQC] version 1.18 or higher for metrics report generation.
- [samtools] version 1.16 or higher for alignment of reads in uBAM or uCRAM format or re-alignment of previously aligned reads.
- [mosdepth] version 0.2.6 or higher for coverage metrics collection from long-read data.

The `sentieon`, `python`, `bcftools`, `bedtools`, `samtools`, `multiqc`, and `mosdepth` executables will be accessed through the user's `PATH` environment variable.

## Input data requirements

### The Reference genome

DNAscope LongRead will call variants present in the sample relative to a high quality reference genome in FASTA format. Besides the reference genome file, a samtools fasta index file (.fai) needs to be present. Short-read alignment also requires bwa index files.

We recommend aligning to a reference genome without alternate contigs. If alternate contigs are present in the genome and the pipeline is performing short-read alignment, please also supply a ".alt" file to activate [alt-aware alignment] in bwa.

## Usage

### Germline variant calling from aligned short and long-read data

A single command is run to call SNVs, indels, SVs and CNVs from aligned short and long reads:
```sh
sentieon-cli dnascope-hybrid \
  -r REFERENCE \
  --sr_aln SR_ALN [SR_ALN ...] \
  --lr_aln LR_ALN [LR_ALN ...] \
  -m MODEL_BUNDLE \
  [-d DBSNP] \
  [-b DIPLOID_BED] \
  [-t NUMBER_THREADS] \
  [--longread_tech TECH] \
  sample.vcf.gz
```

The DNAscope hybrid pipeline requires the following arguments:
- `-r REFERENCE`: the location of the reference FASTA file. A reference fasta index, ".fai" file, is also required.
- `--sr_aln`: the input short-read data in BAM or CRAM format. One or more files can be supplied by passing multiple files after the `--sr_aln` argument.
- `--lr_aln`: the input long-read data in BAM or CRAM format. One or more files can be supplied by passing multiple files after the `--lr_aln` argument.
- `-m MODEL_BUNDLE`: the location of the model bundle. Model bundle files can be found in the [sentieon-models] repository.
- `sample.vcf.gz`: the location of the output VCF file for SNVs and indels. The pipeline requires the output file end with the suffix, ".vcf.gz".

The Sentieon hybrid pipeline accepts the following optional arguments:
- `-d DBSNP`: the location of the Single Nucleotide Polymorphism database (dbSNP) used to label known variants in VCF (`.vcf`) or bgzip compressed VCF (`.vcf.gz`) format. Only one file is supported. Supplying this file will annotate variants with their dbSNP refSNP ID numbers. A VCF index file is required.
- `-b DIPLOID_BED`: interval in the reference to restrict diploid variant calling, in BED file format. Supplying this file will limit diploid variant calling to the intervals inside the BED file.
- `-t NUMBER_THREADS`: number of computing threads that will be used by the software to run parallel processes. The argument is optional; if omitted, the pipeline will use as many threads as the server has.
- `--longread_tech TECH`: the technology used to generate the long-read sequence data. The default is HIFI. With long-read nanopore sequence data, please set `--longread_tech ONT`.
- `-h`: print the command-line help and exit.
- `--dry_run`: print the pipeline commands, but do not actually execute them.

### Germline variant calling from unaligned short and long-read data

A single command is run to call SNVs, indels, SVs and CNVs from unaligned short and long reads:
```sh
sentieon-cli dnascope-hybrid \
  -r REFERENCE \
  --sr_r1_fastq SR_R1_FQ [SR_R1_FQ ...] \
  --sr_r2_fastq SR_R2_FQ [SR_R2_FQ ...] \
  --sr_readgroups SR_READGROUP [SR_READGROUP ...] \
  --lr_aln LR_ALN [LR_ALN ...] \
  --lr_align_input \
  -m MODEL_BUNDLE \
  [-d DBSNP] \
  [-b DIPLOID_BED] \
  [-t NUMBER_THREADS] \
  [--longread_tech TECH] \
  sample.vcf.gz
```

The DNAscope hybrid pipeline requires the following arguments:
- `--sr_r1_fastq`: the input R1 short-read data in gzipped FASTQ format. One or more files can be supplied by passing multiple files after the `--sr_r1_fastq` argument.
- `--sr_r2_fastq`: the input R2 short-read data in gzipped FASTQ format. One or more files can be supplied by passing multiple files after the `--sr_r2_fastq` argument.
- `--sr_readgroups`: readgroup information for each FASTQ. The pipeline will expect the same number of arguments to `--sr_r1_fastq` and `--sr_readgroups`. An example argument is, `--sr_readgroups "@RG\tID:HG002-1\tSM:HG002\tLB:HG002-LB-1\tPL:ILLUMINA"`
- `--lr_aln`: the input long-read data in uBAM or uCRAM format. One or more files can be supplied by passing multiple files after the `--lr_aln` argument.
- `--lr_align_input`: directs the pipeline to align the input long-reads.

The Sentieon hybrid pipeline accepts the following optional arguments:
- `--sr_duplicate_marking`: setting for duplicate marking. `markdup` will mark duplicate reads. `rmdup` will remove duplicate reads. `none` will skip duplicate marking. The default setting is `markdup`.
- `--lr_input_ref`: a reference fasta used for decoding the input long-read file(s). Required with long-read uCRAM or CRAM input. Can be different from the fasta used with the `-r` argument.
- `--bam_format`: use BAM format instead of CRAM for output aligned files.

## Pipeline output

### List of output files

The following files are output by the DNAscope hybrid pipeline:
- `sample.vcf.gz`: SNV and indel variant calls across the regions of the genome as defined in the `-b DIPLOID_BED` file.
- `sample.sv.vcf.gz`: structural variant calls from the Sentieon LongReadSV tool.
- `sample.cnv.vcf.gz`: copy-number variant calls from the Sentieon CNVscope tool.
- `sample_deduped.cram`: aligned, coordinate-sorted and duplicate-marked short-read data from the input FASTQ files.
- `sample_mm2_sorted_*.cram`: aligned and coordinate-sorted long-reads from the input uBAM, uCRAM, BAM, or CRAM files.
- `sample_metrics`: a directory containing QC metrics for the analyzed sample.


[Python]: https://www.python.org/
[bcftools]: http://samtools.github.io/bcftools/bcftools.html
[bedtools]: https://bedtools.readthedocs.io/en/latest/
[MultiQC]: https://multiqc.info/
[mosdepth]: https://github.com/brentp/mosdepth
[samtools]: https://www.htslib.org/
[alt-aware alignment]: https://github.com/lh3/bwa/blob/master/README-alt.md
[sentieon-models]: https://github.com/Sentieon/sentieon-models
