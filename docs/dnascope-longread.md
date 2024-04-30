# DNAscope LongRead

Sentieon DNAscope LongRead is a pipeline for alignment and germline variant calling (SNVs, SVs, and indels) from long-read sequence data. The DNAscope LongRead pipeline is able to take advantage of longer read lengths to perform quick and accurate variant calling using specially calibrated machine learning models.

This folder provides a shell script implementing the DNAscope LongRead pipeline. The pipeline will accept as input aligned reads in BAM or CRAM format, or un-aligned data in FASTQ, uBAM, or uCRAM format. The pipeline will output variants in the VCF (or gVCF) formats and aligned reads in BAM or CRAM formats.

DNAscope LongRead is implemented using the Sentieon software package, which requires a valid license for use. Please contact info@sentieon.com for access to the Sentieon software and an evaluation license.

## Prerequisites

In order to run the pipeline, you need to use the Sentieon software package version 202308.01 or higher and [Python] version >3.8. SNV and indel calling additionally requires [bcftools] version 1.10 or higher and [bedtools]. Alignment of read data in uBAM or uCRAM or re-alignment of read data additional requires [samtools] version 1.16 or higher. The `sentieon`, `python`, `bcftools`, `bedtools`, and `samtools` executables will be accessed through the user's `PATH` environment variable.

## Input data requirements

### The Reference genome

DNAscope LongRead will call variants present in the sample relative to a high quality reference genome sequence. Besides the reference genome file, a samtools fasta index file (.fai) needs to be present. We recommend aligning to a reference genome without alternate contigs.


## Usage

A single command is run to call SNVs, indels, and structural variants from PacBio HiFi or ONT reads in the uBAM, uCRAM, BAM, or CRAM formats:
```sh
sentieon-cli dnascope-longread [-h] -r REFERENCE -i SAMPLE_INPUT \
  -m MODEL_BUNDLE [-d DBSNP] [-b DIPLOID_BED] [-t NUMBER_THREADS] [-g] \
  --tech HiFi|ONT [--haploid-bed HAPLOID_BED] [--] sample.vcf.gz
```

With uBAM, uCRAM, BAM, or CRAM input, the read data can be re-aligned to the reference genome using minimap2 using the `--align` argument.

A single command is run to call SNVs, indels, and structural variants from PacBio HiFi or ONT reads in the FASTQ format:
```sh
sentieon-cli dnascope-longread [-h] -r REFERENCE --fastq INPUT_FASTQ  \
  --readgroups READGROUP -m MODEL_BUNDLE [-d DBSNP] [-b DIPLOID_BED] \
  [-t NUMBER_THREADS] [-g] --tech HiFi|ONT [--haploid-bed HAPLOID_BED] [--] \
  sample.vcf.gz
```

The Sentieon LongRead pipeline requires the following arguments:
- `-r REFERENCE`: the location of the reference FASTA file. You should make sure that the reference is the same as the one used in the mapping stage. A reference fasta index, ".fai" file, is also required.
- `-i SAMPLE_INPUT`: the input sample file in BAM or CRAM format. One or more files can be suppied by passing multiple files after the `-i` argument.
- `--fastq INPUT_FASTQ`: the input sample file in FASTQ format. One or more files can be supplied by passing multiple files after the `--fastq` argument.
- `--readgroups READGROUP`: readgroup information for the read data. The `--readgroups` argument is required if the input data is in the fastq format. This argument expects complete readgroup strings and these strings will be passed to `minimap2` through the `-R` argument. An example argument is `--readgroups '@RG\tID:foo\tSM:bar'`.
- `-m MODEL_BUNDLE`: the location of the model bundle. Model bundle files can be found in the [sentieon-models] repository.
- `--tech HiFi|ONT`: Sequencing technology used to generate the reads. Supported arguments are `ONT` or `HiFi`. 

The Sentieon LongRead pipeline accepts the following optional arguments:
- `-d DBSNP`: the location of the Single Nucleotide Polymorphism database (dbSNP) used to label known variants in VCF (`.vcf`) or bgzip compressed VCF (`.vcf.gz`) format. Only one file is supported. Supplying this file will annotate variants with their dbSNP refSNP ID numbers. A VCF index file is required.
- `-b DIPLOID_BED`: interval in the reference to restrict diploid variant calling, in BED file format. Supplying this file will limit diploid variant calling to the intervals inside the BED file.
- `--haploid-bed HAPLOID_BED`: interval in the reference to restrict haploid variant calling, in BED file format. Supplying this file will perform haploid variant calling across the intervals inside the BED file.
- `-t NUMBER_THREADS`: number of computing threads that will be used by the software to run parallel processes. The argument is optional; if omitted, the pipeline will use as many threads as the server has.
- `-g`: output variants in the gVCF format, in addition to the VCF output file. The tool will output a bgzip compressed gVCF file with a corresponding index file.
- `--align`: re-align the input read data to the reference genome using Sentieon minimap2. This argument needs to be passed to obtain useful output from uBAM or uCRAM input.
- `-h`: print the command-line help and exit.
- `--dry-run`: print the pipeline commands, but do not actually execute them.

The Sentieon LongRead pipeline requires the following positional arguments:
- `sample.vcf.gz`: the location and filename of the variant calling output. The filename must end with the suffix, `.vcf.gz`. The tool will output a bgzip compressed VCF file with a corresponding index file.

## Pipeline output

The Sentieon LongRead pipeline will output a bgzip compressed file (.vcf.gz) containing variant calls in the standard VCFv4.2 format along with a tabix index file (.vcf.gz.tbi). If the `-g` option is used, the pipeline will also output a bgzip compressed file (.g.vcf.gz) containing variant calls in the gVCF format along with a tabix index file (.g.vcf.gz.tbi).

With the `--haploid-bed HAPLOID_BED` argument, the pipeline will create the following output files:
- `sample.vcf.gz`: SNV and indel variant calls across the diploid regions of the genome (as defined in the `-d DIPLOID_BED` file).
- `sample.haploid.vcf.gz`: SNV and indel variant calls across the haploid regions of the genome (as defined in the `--haploid-bed HAPLOID_BED` file).
- `sample.sv.vcf.gz`: structural variant calls from the Sentieon LongReadSV tool.

With the `--align` argument or fastq input, the software will also output a CRAM or BAM file for each input file.

## Other considerations

### Diploid and haploid variant calling

The default pipeline is recommended for use with samples from diploid organisms. For samples with both diploid and haploid chromosomes, the `-b DIPLOID_BED` option can be used to limit diploid variant calling to diploid chromosomes and the `--haploid-bed HAPLOID_BED` argument can be used to perform haploid variant calling across haploid chromosomes. Diploid and haploid variants will be output to separate VCF files.

Example diploid and haploid BED files can be found in the [/data](/data) folder in this repository.

### Modification

Scripts in this folder are made available under the [BSD 2-Clause license](/LICENSE). Modification of the shell scripts in this folder is encouraged for users who wish to retain intermediate files generated by the pipeline, add additional processing steps, or modify command line arguments.

The Python scripts in the `sentieon_cli/scripts`` folder perform low-level manipulation of intermediate gVCF and VCF files generated by the pipeline. Due to the low-level data handling performed by these scripts, modification of these files by users is discouraged.

## References
**[Sentieon DNAscope LongRead – A highly Accurate, Fast, and Efficient Pipeline for Germline Variant Calling from PacBio HiFi reads]** - A preprint describing the DNAscope LongRead pipeline for calling variants from PacBio HiFi data.


[Python]: https://www.python.org/
[bcftools]: http://samtools.github.io/bcftools/bcftools.html
[bedtools]: https://bedtools.readthedocs.io/en/latest/
[samtools]: https://www.htslib.org/
[sentieon-models]: https://github.com/Sentieon/sentieon-models

[Sentieon DNAscope LongRead – A highly Accurate, Fast, and Efficient Pipeline for Germline Variant Calling from PacBio HiFi reads]: https://www.biorxiv.org/content/10.1101/2022.06.01.494452v1