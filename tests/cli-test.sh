#!/bin/bash

set -eu o pipefail

mkdir -p /tmp/empty_bin
echo "#!/bin/sh" > /tmp/empty_bin/sentieon
echo "#!/bin/sh" > /tmp/empty_bin/bedtools
echo "#!/bin/sh" > /tmp/empty_bin/bedtools
chmod +x /tmp/empty_bin/*

touch reference.fasta
touch sample.input
touch model.file
touch dbsnp.vcf.gz
touch regions.bed
export PATH=/tmp/empty_bin/:$PATH

bash tests/fake_run.sh -r reference.fasta -i sample.input -m model.file -d dbsnp.vcf.gz -b regions.bed -t 4 out.vcf.gz \
    2>&1 | perl -pe 's/^\s*\++\s//' > cmds

exp=$(grep ^senti cmds | head -1)

poetry install

obs=$(sentieon_do run-algo-dnascope out.vcf.gz -r reference.fasta -i sample.input -m model.file -d dbsnp.vcf.gz  -b regions.bed -t 4 2>&1 | tail -1)

echo EXP
echo $exp
echo OBS
echo $obs

