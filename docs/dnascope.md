# DNAscope

Sentieon DNAscope is a pipeline for alignment and germline variant calling (SNVs, SVs and indels) from short-read DNA sequence data. The DNAscope pipeline uses a combination of traditional statistical approaches and machine learning to achieve high variant calling accuracy. The DNAscope pipeline supports samples sequenced using whole-genome or targeted (hybrid-capture) enrichment library preps.

The pipeline accepts as input aligned reads in BAM or CRAM format, or un-aligned reads in FASTQ, uBAM, or uCRAM format. The pipeline will output variants in the VCF (or gVCF) formats and aligned reads in BAM or CRAM formats.

DNAscope is implemented using the Sentieon software package, which requires a valid license for use. Please contact info@sentieon.com for access to the Sentieon software and an evaluation license.

## Prerequisites

- Sentieon software package version 202308 or higher.
- [samtools] version 1.16 or higher for alignment of reads in uBAM or uCRAM format or re-alignment of previously aligned reads.
- [MultiQC] version 1.18 or higher for metrics report generation.

The `sentieon`, `samtools`, and `multiqc` executables will be accessed through the user's `PATH` environment variable.

## Input data requirements

### The Reference genome

DNAscope will call variants present in the sample relative to a high quality reference genome sequence. Besides the reference genome file, a samtools fasta index file (.fai) needs to be present. Read alignment also requires bwa index files.

We recommend aligning to a reference genome without alternate contigs. If alternate contigs are present in the genome, please also supply a ".alt" file to activate [alt-aware alignment] in bwa.

## Usage

### Alignment and variant calling from FASTQ

A single command is run to align, preprocess, and call SNVs, indels, and structural variants from FASTQ:
```sh
sentieon-cli dnascope [-h] \
  -r REFERENCE \
  --r1_fastq R1_FASTQ ... \
  --r2_fastq R2_FASTQ ... \
  --readgroups READGROUPS ... \
  -m MODEL_BUNDLE \
  [-d DBSNP] \
  [-b BED] \
  [--interval_padding INTERVAL_PADDING] \
  [-t NUMBER_THREADS] \
  [--pcr_free] \
  [-g] \
  [--duplicate_marking DUP_MARKING] \
  [--assay ASSAY] \
  [--consensus] \
  [--dry_run] \
  [--bam_format] \
  sample.vcf.gz
```

With FASTQ input, the DNAscope pipeline requires the following arguments:
- `-r REFERENCE`: the location of the reference FASTA file. A reference fasta index, ".fai" file, and bwa index files, are also required.
- `--r1_fastq R1_FASTQ`: the R1 input FASTQ. Can be used multiple times. `--r1_fastq` files without a corresponding `--r2_fastq` are assumed to be single-ended. Be aware that the pipeline performs single-sample processing, and all fastq are expected to be from the same sample.
- `--r2_fastq R2_FASTQ`: the R2 input FASTQ. Can be used multiple times.
- `--readgroups READGROUPS`: readgroup information for each FASTQ. The pipeline will expect the same number of arguments to `--r1_fastq` and `--readgroups`. An example argument is, `--readgroups "@RG\tID:HG002-1\tSM:HG002\tLB:HG002-LB-1\tPL:ILLUMINA"`
- `-m MODEL_BUNDLE`: the location of the model bundle. Model bundle files can be found in the [sentieon-models] repository.
- `sample.vcf.gz`: the location of the output VCF file for SNVs and indels. The pipeline requires the output file end with the suffix, ".vcf.gz". The file path without the suffix will be used as the basename for other output files.

The DNAscope pipeline accepts the following optional arguments:
- `-d DBSNP`: the location of the Single Nucleotide Polymorphism database (dbSNP) used to label known variants in VCF (`.vcf`) or bgzip compressed VCF (`.vcf.gz`) format. Only one file is supported. Supplying this file will annotate variants with their dbSNP refSNP ID numbers. A VCF index file is required.
- `-b BED`: interval in the reference to restrict variant calling, in BED file format. Supplying this file will limit variant calling to the intervals inside the BED file. If a BED file is not supplied, the software will process the whole genome.
- `--interval_padding INTERVAL_PADDING`: adds INTERVAL_PADDING bases padding to the edges of the input intervals. The default value is 0.
- `-t NUMBER_THREADS`: number of computing threads that will be used by the software to run parallel processes. The argument is optional; if omitted, the pipeline will use as many threads as the server has.
- `--pcr_free`: Call variants using `--pcr_indel_model NONE`, which is appropriate for libraries prepared with a PCR-free library prep. Deduplication is still performed to identify optical duplicates.
- `-g`: output variants in the gVCF format, in addition to the VCF output file. The tool will output a bgzip compressed gVCF file with a corresponding index file.
- `--duplicate_marking DUP_MARKING`: setting for duplicate marking. `markdup` will mark duplicate reads. `rmdup` will remove duplicate reads. `none` will skip duplicate marking. The default setting is `markdup`.
- `--assay ASSAY`: assay setting for metrics collection `WGS` or `WES`. The default setting is `WGS`.
- `--consensus`: generate consensus reads during duplicate marking.
- `-h`: print the command-line help and exit.
- `--dry_run`: print the pipeline commands, but do not actually execute them.
- `--bam_format`: use BAM format instead of CRAM for output aligned files.

### Alignment and variant calling from uBAM or uCRAM
A single command is run to align, preprocess, and call SNVs, indels, and structural variants from uBAM or uCRAM files:
```sh
sentieon-cli dnascope [-h] \
  -r REFERENCE \
  -i SAMPLE_INPUT ... \
  --align \
  [--input_ref INPUT_REF] \
  -m MODEL_BUNDLE \
  [-d DBSNP] \
  [-b BED] \
  [--interval_padding INTERVAL_PADDING] \
  [-t NUMBER_THREADS] \
  [--pcr_free] \
  [-g] \
  [--duplicate_marking DUP_MARKING] \
  [--assay ASSAY] \
  [--consensus] \
  [--dry_run] \
  [--bam_format] \
  sample.vcf.gz
```

With uBAM or uCRAM input, the DNAscope pipeline requires the following new arguments:
- `-i SAMPLE_INPUT`: the input sample file in uBAM or uCRAM format. One or more files can be supplied by passing multiple files after the `-i` argument.
- `--align`: directs the pipeline to align the input reads.

The DNAscope pipeline accepts the following new optional arguments:
- `--input_ref INPUT_REF`: a reference fasta used for decoding the input file(s). Required with uCRAM input. Can be different from the fasta used with the `-r` argument.

### Alignment and variant calling from sorted BAM or CRAM
A single command is run to align, preprocess, and call SNVs, indels, and structural variants from BAM or CRAM files:
```sh
sentieon-cli dnascope [-h] \
  -r REFERENCE \
  -i SAMPLE_INPUT ... \
  --collate_align \
  [--input_ref INPUT_REF] \
  -m MODEL_BUNDLE \
  [-d DBSNP] \
  [-b BED] \
  [--interval_padding INTERVAL_PADDING] \
  [-t NUMBER_THREADS] \
  [--pcr_free] \
  [-g] \
  [--duplicate_marking DUP_MARKING] \
  [--assay ASSAY] \
  [--consensus] \
  [--dry_run] \
  [--bam_format] \
  sample.vcf.gz
```

With BAM or CRAM input, the DNAscope pipeline requires the following new arguments:
- `--collate_align`: directs the pipeline to collate and then align the input reads.

### Variant calling from sorted BAM or CRAM
A single command is run to preprocess, and call SNVs, indels, and structural variants from BAM or CRAM files:
```sh
sentieon-cli dnascope [-h] \
  -r REFERENCE \
  -i SAMPLE_INPUT ... \
  -m MODEL_BUNDLE \
  [-d DBSNP] \
  [-b BED] \
  [--interval_padding INTERVAL_PADDING] \
  [-t NUMBER_THREADS] \
  [--pcr_free] \
  [-g] \
  [--duplicate_marking DUP_MARKING] \
  [--assay ASSAY] \
  [--consensus] \
  [--dry_run] \
  [--bam_format] \
  sample.vcf.gz
```

Not supplying the `--align` and `--collate_align` arguments will direct the pipeline to call variants directly from the input reads.

## Pipeline output

### List of output files

The following files are output when processing WGS FASTQ with default arguments:
- `sample.vcf.gz`: SNV and indel variant calls across the regions of the genome as defined in the `-b BED` file.
- `sample_deduped.cram` or `sample_deduped.bam`: aligned, coordinate-sorted and duplicate-marked read data from the input FASTQ files.
- `sample_svs.vcf.gz`: structural variant calls from DNAscope and SVSolver.
- `sample_metrics`: a directory containing QC metrics for the analyzed sample. 
  - `sample_metrics/coverage*`: coverage metrics for the processed sample. Only available for WGS samples.
  - `sample_metrics/{sample}.txt.alignment_stat.txt`: Metrics from the AlignmentStat algo.
  - `sample_metrics/{sample}.txt.base_distribution_by_cycle.txt`: Metrics from the BaseDistributionByCycle algo.
  - `sample_metrics/{sample}.txt.dedup_metrics.txt`: Metrics from the Dedup algo.
  - `sample_metrics/{sample}.txt.gc_bias*`: Metrics from the GCBias algo. Only available for WGS samples.
  - `sample_metrics/{sample}.txt.insert_size.txt`: Metrics from the InsertSizeMetricAlgo algo.
  - `sample_metrics/{sample}.txt.mean_qual_by_cycle.txt`: Metrics from the MeanQualityByCycle algo.
  - `sample_metrics/{sample}.txt.qual_distribution.txt`: Metrics from the QualDistribution algo.
  - `sample_metrics/{sample}.txt.wgs.txt`: Metrics from the WgsMetricsAlgo algo. Only available for WGS samples.
  - `sample_metrics/{sample}.txt.hybrid-selection.txt`: Metrics from the HsMetricAlgo algo.
  - `sample_metrics/multiqc_report.html`: collected QC metrics aggregated by MultiQC.

[samtools]: https://www.htslib.org/
[MultiQC]: https://multiqc.info/
[alt-aware alignment]: https://github.com/lh3/bwa/blob/master/README-alt.md
[sentieon-models]: https://github.com/Sentieon/sentieon-models
