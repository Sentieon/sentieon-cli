# DNAscope LongRead

Sentieon DNAscope LongRead is a pipeline for alignment and germline variant calling (SNVs, SVs, and indels) from long-read sequence data. The DNAscope LongRead pipeline is able to take advantage of longer read lengths to perform quick and accurate variant calling using specially calibrated machine learning models.

The pipeline will accept as input aligned reads in BAM or CRAM format, or un-aligned reads in FASTQ, uBAM, or uCRAM format. The pipeline will output variants in the VCF (or gVCF) formats and aligned reads in BAM or CRAM formats.

DNAscope LongRead is implemented using the Sentieon software package, which requires a valid license for use. Please contact info@sentieon.com for access to the Sentieon software and an evaluation license.

## Prerequisites

- Sentieon software package version 202308.01 or higher.
- [Python] version 3.8 or higher.
- [bcftools] version 1.10 or higher for SNV and indel calling.
- [bedtools] for SNV and indel calling.
- [samtools] version 1.16 or higher for alignment of read data in uBAM or uCRAM format or re-alignment of previously aligned reads.

The `sentieon`, `python`, `bcftools`, `bedtools`, and `samtools` executables will be accessed through the user's `PATH` environment variable.

## Input data requirements

### The Reference genome

DNAscope LongRead will call variants present in the sample relative to a high quality reference genome sequence. Besides the reference genome file, a samtools fasta index file (.fai) needs to be present.

We recommend aligning to a reference genome without alternate contigs.


## Usage

### Alignment and variant calling from FASTQ

A single command is run to align and call SNVs, indels, and structural variants from PacBio HiFi or ONT reads in the FASTQ format:
```sh
sentieon-cli dnascope-longread [-h] \
  -r REFERENCE \
  --fastq INPUT_FASTQ ... \
  --readgroups READGROUP ... \
  -m MODEL_BUNDLE \
  [-d DBSNP] \
  [-b DIPLOID_BED] \
  [-t NUMBER_THREADS] \
  [-g] \
  --tech HiFi|ONT \
  [--haploid-bed HAPLOID_BED] \
  sample.vcf.gz
```

With FASTQ input, the DNAscope LongRead pipeline requires the following arguments:
- `-r REFERENCE`: the location of the reference FASTA file. A reference fasta index, ".fai" file, is also required.
- `--fastq INPUT_FASTQ`: the input sample file in FASTQ format. One or more files can be supplied by passing multiple files after the `--fastq` argument.
- `--readgroups READGROUP`: readgroup information for the read data. The `--readgroups` argument is required if the input data is in the FASTQ format. This argument expects complete readgroup strings and these strings will be passed to `minimap2` through the `-R` argument. An example argument is `--readgroups '@RG\tID:foo\tSM:bar'`.
- `-m MODEL_BUNDLE`: the location of the model bundle. Model bundle files can be found in the [sentieon-models] repository.
- `--tech HiFi|ONT`: Sequencing technology used to generate the reads. Supported arguments are `ONT` or `HiFi`.
- `sample.vcf.gz`: the location of the output VCF file for SNVs and indels. The pipeline requires the output file end with the suffix, ".vcf.gz". The file path without the suffix will be used as the basename for other output files.

The Sentieon LongRead pipeline accepts the following optional arguments:
- `-d DBSNP`: the location of the Single Nucleotide Polymorphism database (dbSNP) used to label known variants in VCF (`.vcf`) or bgzip compressed VCF (`.vcf.gz`) format. Only one file is supported. Supplying this file will annotate variants with their dbSNP refSNP ID numbers. A VCF index file is required.
- `-b DIPLOID_BED`: interval in the reference to restrict diploid variant calling, in BED file format. Supplying this file will limit diploid variant calling to the intervals inside the BED file.
- `--haploid-bed HAPLOID_BED`: interval in the reference to restrict haploid variant calling, in BED file format. Supplying this file will perform haploid variant calling across the intervals inside the BED file.
- `-t NUMBER_THREADS`: number of computing threads that will be used by the software to run parallel processes. The argument is optional; if omitted, the pipeline will use as many threads as the server has.
- `-g`: output variants in the gVCF format, in addition to the VCF output file. The tool will output a bgzip compressed gVCF file with a corresponding index file.
- `-h`: print the command-line help and exit.
- `--dry-run`: print the pipeline commands, but do not actually execute them.

### Alignment and variant calling from uBAM, uCRAM, BAM, or CRAM

A single command is run to align and call SNVs, indels, and structural variants from PacBio HiFi or ONT reads in the uBAM, uCRAM, BAM, or CRAM formats:
```sh
sentieon-cli dnascope-longread [-h] \
  -r REFERENCE \
  -i SAMPLE_INPUT ... \
  --align \
  -m MODEL_BUNDLE \
  [-d DBSNP] \
  [-b DIPLOID_BED] \
  [-t NUMBER_THREADS] \
  [-g] \
  --tech HiFi|ONT \
  [--haploid-bed HAPLOID_BED] \
  [--input_ref INPUT_REF] \
  sample.vcf.gz
```

With uBAM, uCRAM, BAM, or CRAM input, the DNAscope LongRead pipeline requires the following new arguments:
- `-i SAMPLE_INPUT`: the input input sample file in uBAM or uCRAM format. One or more files can be supplied by passing multiple files after the `-i` argument.
- `--align`: re-align the input read data to the reference genome using Sentieon minimap2.

The DNAscope LongRead pipeline accepts the following new optional arguments:
- `--input_ref INPUT_REF`: a reference fasta used for decoding the input file(s). Required with uCRAM or CRAM input. Can be different from the fasta used with the `-r` argument.

### Variant calling from BAM or CRAM

A single command is run to call SNVs, indels, and structural variants from PacBio HiFi or ONT reads in the BAM, or CRAM formats:
```sh
sentieon-cli dnascope-longread [-h] \
  -r REFERENCE \
  -i SAMPLE_INPUT ... \
  -m MODEL_BUNDLE \
  [-d DBSNP] \
  [-b DIPLOID_BED] \
  [-t NUMBER_THREADS] \
  [-g] \
  --tech HiFi|ONT \
  [--haploid-bed HAPLOID_BED] \
  sample.vcf.gz
```

Not supplying the `--align` argument will direct the pipeline to call variants directly from the input reads.

## Pipeline output

### List of output files

The following files are output when processing FASTQ data or uBAM, uCRAM, BAM, or CRAM files with the `--align` argument:
- `sample.vcf.gz`: SNV and indel variant calls across the regions of the genome as defined in the `-b DIPLOID_BED` file.
- `sample.sv.vcf.gz`: structural variant calls from the Sentieon LongReadSV tool.
- `sample_mm2_sorted_fq_*.cram`: aligned and coordinate-sorted reads from the input FASTQ files.
- `sample_mm2_sorted_*.cram`: aligned and coordinate-sorted reads from the input uBAM, uCRAM, BAM, or CRAM files.

With the `--haploid-bed HAPLOID_BED` argument, the pipeline will create the following additional output files:
- `sample.haploid.vcf.gz`: SNV and indel variant calls across the haploid regions of the genome (as defined in the `--haploid-bed HAPLOID_BED` file).

## Other considerations

### Diploid and haploid variant calling

The default pipeline is recommended for use with samples from diploid organisms. For samples with both diploid and haploid chromosomes, the `-b DIPLOID_BED` option can be used to limit diploid variant calling to diploid chromosomes and the `--haploid-bed HAPLOID_BED` argument can be used to perform haploid variant calling across haploid chromosomes. Diploid and haploid variants will be output to separate VCF files.

Diploid and haploid BED files for the human hg38 reference genome (with male samples) can be found in the [/data](/data) folder in this repository.

### Modification

Scripts in this repository are made available under the [BSD 2-Clause license](/LICENSE).

The Python scripts in the `sentieon_cli/scripts` folder perform low-level manipulation of intermediate gVCF and VCF files generated by the pipeline. Due to the low-level data handling performed by these scripts, modification of these files by users is discouraged.

## References
**[Sentieon DNAscope LongRead – A highly Accurate, Fast, and Efficient Pipeline for Germline Variant Calling from PacBio HiFi reads]** - A preprint describing the DNAscope LongRead pipeline for calling variants from PacBio HiFi data.


[Python]: https://www.python.org/
[bcftools]: http://samtools.github.io/bcftools/bcftools.html
[bedtools]: https://bedtools.readthedocs.io/en/latest/
[samtools]: https://www.htslib.org/
[sentieon-models]: https://github.com/Sentieon/sentieon-models

[Sentieon DNAscope LongRead – A highly Accurate, Fast, and Efficient Pipeline for Germline Variant Calling from PacBio HiFi reads]: https://www.biorxiv.org/content/10.1101/2022.06.01.494452v1