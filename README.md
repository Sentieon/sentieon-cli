[![CI](https://github.com/brentp/sentieon-cli/actions/workflows/ci.yml/badge.svg)](https://github.com/brentp/sentieon-cli/actions/workflows/ci.yml)

# Sentieon CLI


## Setup
At present, manually add vcf_mod.py and gvcf_combine.py from the [sentieon-scripts repo](

```
wget -O sentieon_cli/gvcf_combine.py https://github.com/Sentieon/sentieon-scripts/raw/master/dnascope_LongRead/gvcf_combine.py
wget -O sentieon_cli/vcf_mod.py https://github.com/Sentieon/sentieon-scripts/raw/master/dnascope_LongRead/vcf_mod.py
```

Then run:
```
pip install poetry
poetry shell
poetry install
sentieon-cli
```
