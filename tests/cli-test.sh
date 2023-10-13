#!/bin/bash

set -eu o pipefail

mkdir -p /tmp/empty_bin
echo -e "#!/bin/sh\necho sentieon \$@" > /tmp/empty_bin/sentieon
echo -e "#!/bin/sh\necho bedtools \$@" > /tmp/empty_bin/bedtools
echo -e "#!/bin/sh\necho bcftools \$@" > /tmp/empty_bin/bcftools
chmod +x /tmp/empty_bin/*

touch reference.fasta
touch sample.input
touch model.file
touch dbsnp.vcf.gz
touch regions.bed
export PATH=/tmp/empty_bin/:$PATH

#wget -O tests/dnascope_HiFi.sh https://raw.githubusercontent.com/Sentieon/sentieon-scripts/master/dnascope_LongRead/dnascope_HiFi.sh
#wget -O tests/dnascope_ONT.sh https://raw.githubusercontent.com/Sentieon/sentieon-scripts/master/dnascope_LongRead/dnascope_ONT.sh
export SENTIEON_TMPDIR=/tmp/asdf

bash tests/dnascope_HiFi.sh -g -r reference.fasta -i sample.input -m model.file -d dbsnp.vcf.gz -b regions.bed -t 4 out.vcf.gz \
     | perl -pe 's/^\s*\++\s//' > hifi_cmds
sentieon-cli run-full-dnascope out.vcf.gz -r reference.fasta -i sample.input -m model.file -d dbsnp.vcf.gz  -b regions.bed --tech HiFi -t 4 2>err.hifi > my-hifi_cmds

bash tests/dnascope_ONT.sh -g -r reference.fasta -i sample.input -m model.file -d dbsnp.vcf.gz -b regions.bed -t 4 out.vcf.gz \
     | perl -pe 's/^\s*\++\s//' > ont_cmds
sentieon-cli run-full-dnascope out.vcf.gz -r reference.fasta -i sample.input -m model.file -d dbsnp.vcf.gz  -b regions.bed --tech ONT -t 4 2>err.ont > my-ont_cmds

