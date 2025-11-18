# Pangenome

Pangenome is a pipeline for alignment and variant calling from short-read DNA sequence data using pangenome graphs. The Pangenome pipeline leverages graph-based reference representations to improve alignment and variant calling accuracy, particularly in complex genomic regions with high sequence diversity.

The pipeline uses the open-source vg toolkit for pangenome alignment and structural variant calling. Sentieon tools are used for small variant calling, CNV detection, segdup variant calling, preprocessing and metrics collection. Specialized tools are used for HLA/KIR genotyping and repeat expansion analysis.

The pipeline accepts as input unaligned reads in FASTQ format and outputs variants in VCF format, aligned reads in BAM format, and comprehensive quality control metrics.

Pangenome is implemented using open-source tools and the Sentieon software package, which requires a valid license for use. Please contact info@sentieon.com for access to the Sentieon software and an evaluation license.

## Prerequisites

- Sentieon software package version 202503 or higher.
- [vg] toolkit for pangenome alignment operations.
- [KMC] version 3 or higher for k-mer counting.
- [samtools] version 1.16 or higher for alignment operations.
- [bcftools] version 1.10 or higher for VCF operations.
- [MultiQC] version 1.18 or higher for metrics report generation.
- [T1K] for HLA and KIR genotyping (optional).
- [ExpansionHunter] for repeat expansion calling (optional).
- [segdup-caller] for segmental duplication variant calling (optional).

The executables will be accessed through the user's `PATH` environment variable.

## Input data requirements

### The Reference genome

Pangenome will call variants present in the sample relative to a high quality reference genome sequence. Besides the reference genome file, a samtools fasta index file (.fai) needs to be present.

### Pangenome graph files  

The pipeline requires several pangenome graph files:
- **GBZ file**: The pangenome graph in GBZ format
- **Haplotype file**: Haplotype information for the pangenome
- **Snarls file**: VG snarls file for the GBZ file (created using `vg snarls`)
- **XG file**: VG XG index file for the GBZ file (created using `vg convert`)

### Model bundle

A Sentieon model bundle containing machine learning models for variant calling is required. Model bundle files can be found in the [sentieon-models] repository.

## Usage

### Pangenome alignment and variant calling from FASTQ

A single command is run to align reads to a pangenome graph, call small variants, structural variants, and CNVs:

```sh
sentieon-cli pangenome [-h] \
  -r REFERENCE \
  --gbz GBZ \
  --hapl HAPL \
  --snarls SNARLS \
  --xg XG \
  -m MODEL_BUNDLE \
  --r1_fastq R1_FASTQ ... \
  --r2_fastq R2_FASTQ ... \
  --readgroups READGROUPS ... \
  [-d DBSNP] \
  [-t CORES] \
  [--kmer_memory KMER_MEMORY] \
  [--expansion_catalog EXPANSION_CATALOG] \
  [--t1k_hla_seq T1K_HLA_SEQ] \
  [--t1k_hla_coord T1K_HLA_COORD] \
  [--t1k_kir_seq T1K_KIR_SEQ] \
  [--t1k_kir_coord T1K_KIR_COORD] \
  [--segdup_caller_genes SEGDUP_CALLER_GENES] \
  [--dry_run] \
  sample.vcf.gz
```

The Pangenome pipeline requires the following arguments:
- `-r REFERENCE`: the location of the reference FASTA file. A reference fasta index, ".fai" file, is also required.
- `--gbz GBZ`: the pangenome graph file in GBZ format.
- `--hapl HAPL`: the haplotype file for the pangenome.
- `--snarls SNARLS`: the vg snarls file for the GBZ file.
- `--xg XG`: the XG file for the GBZ file.
- `-m MODEL_BUNDLE`: the location of the model bundle containing DNAscope and CNVscope models.
- `--r1_fastq R1_FASTQ`: the R1 input FASTQ. Can be used with multiple files.
- `--r2_fastq R2_FASTQ`: the R2 input FASTQ. Can be used with multiple files.
- `--readgroups READGROUPS`: readgroup information for each FASTQ. An example argument is, `--readgroups "@RG\tID:HG002-1\tSM:HG002\tLB:HG002-LB-1\tPL:ILLUMINA"`.
- `sample.vcf.gz`: the location of the output VCF file for small variants. The pipeline requires the output file end with the suffix, ".vcf.gz". The file path without the suffix will be used as the basename for other output files.

The Pangenome pipeline accepts the following optional arguments:
- `-d DBSNP`: the location of the Single Nucleotide Polymorphism database (dbSNP) VCF used to label known variants. A VCF index file is required.
- `-t CORES`: number of computing threads/cores to use. Default is all available cores.
- `--kmer_memory KMER_MEMORY`: memory limit for KMC k-mer counting in GB. Default is 128 GB.
- `--expansion_catalog EXPANSION_CATALOG`: ExpansionHunter variant catalog for repeat expansion calling.
- `--t1k_hla_seq T1K_HLA_SEQ`: DNA HLA sequence FASTA file for T1K HLA calling. Requires `--t1k_hla_coord`.
- `--t1k_hla_coord T1K_HLA_COORD`: DNA HLA coordinate FASTA file for T1K HLA calling. Requires `--t1k_hla_seq`.
- `--t1k_kir_seq T1K_KIR_SEQ`: DNA KIR sequence FASTA file for T1K KIR calling. Requires `--t1k_kir_coord`.
- `--t1k_kir_coord T1K_KIR_COORD`: DNA KIR coordinate FASTA file for T1K KIR calling. Requires `--t1k_kir_seq`.
- `--segdup_caller_genes SEGDUP_CALLER_GENES`: Genes for segmental duplication calling. Example: 'CFH,CFHR3,CYP11B1,CYP2D6,GBA,NCF1,PMS2,SMN1,STRC'.
- `-h`: print the command-line help and exit.
- `--dry_run`: print the pipeline commands, but do not actually execute them.

## Pipeline workflow

The Pangenome pipeline consists of multiple stages executed across two DAGs:

### First DAG
1. **K-mer counting**: `kmc` counts k-mers from input FASTQ files
2. **Sample-specific pangenome**: `vg haplotypes` creates a sample-specific pangenome using k-mer frequencies
3. **Pangenome alignment**: `vg giraffe` aligns reads to the sample-specific pangenome
4. **Read support computation**: `vg pack` computes read support for variants
5. **Structural variant calling**: `vg call` identifies structural variants from read support
6. **Linear alignment**: `vg surject` projects pangenome alignments back to linear reference
7. **Duplicate marking**: Sentieon `Dedup` marks duplicate reads and computes metrics
8. **Indel realignment**: Sentieon `Realigner` performs indel realignment
9. **Small variant calling**: Sentieon `DNAscope` calls small variants (SNVs/indels)
10. **Ploidy estimation**: Estimation of sample sex and ploidy for downstream processing
11. **HLA/KIR genotyping**: T1K genotypes HLA and KIR loci (optional)
12. **Segmental duplication calling**: segdup-caller identifies variants in segmental duplications (optional)

### Second DAG  
1. **CNV calling**: Sentieon CNVscope calls copy number variants with sex-specific processing
2. **Repeat expansion calling**: ExpansionHunter identifies repeat expansions (optional)

## Pipeline output

### List of output files

The following files are output when processing WGS FASTQ:
- `sample.vcf.gz`: Small variant calls (SNVs and indels) from DNAscope.
- `sample_pangenome-aligned.bam`: Aligned, coordinate-sorted and duplicate-marked read data projected to the linear reference.
- `sample_svs.vcf.gz`: Structural variant calls from vg call.
- `sample_cnv.vcf.gz`: Copy number variant calls from CNVscope.
- `sample_ploidy.json`: Estimated sample sex and ploidy data.
- `sample_metrics`: A directory containing QC metrics for the analyzed sample.
  - `sample_metrics/{sample}.txt.alignment_stat.txt`: Metrics from the AlignmentStat algo.
  - `sample_metrics/{sample}.txt.base_distribution_by_cycle.txt`: Metrics from the BaseDistributionByCycle algo.
  - `sample_metrics/{sample}.txt.dedup_metrics.txt`: Metrics from the Dedup algo.
  - `sample_metrics/{sample}.txt.gc_bias*`: Metrics from the GCBias algo.
  - `sample_metrics/{sample}.txt.insert_size.txt`: Metrics from the InsertSizeMetricAlgo algo.
  - `sample_metrics/{sample}.txt.mean_qual_by_cycle.txt`: Metrics from the MeanQualityByCycle algo.
  - `sample_metrics/{sample}.txt.qual_distribution.txt`: Metrics from the QualDistribution algo.
  - `sample_metrics/{sample}.txt.wgs.txt`: Metrics from the WgsMetricsAlgo algo.
  - `sample_metrics/multiqc_report.html`: Collected QC metrics aggregated by MultiQC.

### Optional output files

When optional analyses are enabled:
- `sample_hla/`: HLA genotyping results from T1K (when HLA files are provided).
- `sample_kir/`: KIR genotyping results from T1K (when KIR files are provided).
- `sample_segdups/`: Segmental duplication variant calls (when genes are specified).
- `sample_expansion*`: Repeat expansion calls from ExpansionHunter (when a catalog is provided).

[vg]: https://github.com/vgteam/vg
[KMC]: https://github.com/refresh-bio/KMC
[samtools]: https://www.htslib.org/
[bcftools]: http://samtools.github.io/bcftools/bcftools.html
[MultiQC]: https://multiqc.info/
[T1K]: https://github.com/mourisl/T1K
[ExpansionHunter]: https://github.com/Illumina/ExpansionHunter
[segdup-caller]: https://github.com/Sentieon/segdup-caller
[sentieon-models]: https://github.com/Sentieon/sentieon-models

## Example data processing

### Obtaining the pangenome graph files

Download vg the public human pangenome from the HPRC repository

```sh
curl -LO "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz"
curl -LO "https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.hapl"
```

Generate a snarls file from the .gbz with `vg snarls`.
```sh
vg snarls hprc-v1.1-mc-grch38.gbz > hprc-v1.1-mc-grch38.snarls
```

Generate a xg file from the .gbz with `vg convert`.
```sh
vg convert -x --drop-haplotypes hprc-v1.1-mc-grch38.gbz > hprc-v1.1-mc-grch38.xg
```

### Download the hg38 reference genome

Download the reference fasta file and samtools faidx file
```sh
curl -L 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz' | \
  gzip -dc > GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
curl -LO 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai'
```

### Download an expansion catalog of clinically relevant expansions

```sh
curl -L 'https://github.com/Illumina/ExpansionHunter/raw/refs/tags/v5.0.0/variant_catalog/hg38/variant_catalog.json' > hg38_variant_catalog.json
```

### Download t1k fasta files

Download the two databases
```sh
perl t1k-build.pl -o hlaidx --download IPD-IMGT/HLA
perl t1k-build.pl -o kiridx --download IPD-KIR 
```

Download the gencode gtf file for hg38
```sh
curl -L 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz' | \
  gzip -dc > gencode.v38.annotation.gtf
```

Bulid the coordinate file
```sh
perl t1k-build.pl -o hlaidx -d hlaidx/hla.dat -g gencode.v38.annotation.gtf
perl t1k-build.pl -o kiridx -d kiridx/kir.dat -g gencode.v38.annotation.gtf
```

### Download a 30x WGS fastq

```sh
curl -LO 'https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.R1.fastq.gz'
curl -LO 'https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.R2.fastq.gz'
```

### Run the pangenome pipeline

```sh
sentieon-cli -v pangenome \
  --reference GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --gbz hprc-v1.1-mc-grch38.gbz \
  --hapl hprc-v1.1-mc-grch38.hapl \
  --snarls hprc-v1.1-mc-grch38.snarls \
  --xg hprc-v1.1-mc-grch38.xg \
  --model_bundle pangenome.bundle \
  --r1_fastq HG002.novaseq.pcr-free.30x.R1.fastq.gz \
  --r2_fastq HG002.novaseq.pcr-free.30x.R2.fastq.gz \
  --readgroups '@RG\tID:HG002-1\tSM:HG002\tPL:ILLUMINA' \
  --expansion_catalog hg38_variant_catalog.json \
  --segdup_caller_genes 'CFH,CFHR3,CYP11B1,CYP2D6,GBA,NCF1,PMS2,SMN1,STRC' \
  --t1k_hla_seq hlaidx/data_dna_seq.fa \
  --t1k_hla_coord hlaidx/data_dna_coord.fa \
  --t1k_kir_seq kiridx/data_dna_seq.fa \
  --t1k_kir_coord kiridx/data_dna_coord.fa \
  HG002_pangenome_analysis.vcf.gz
```
