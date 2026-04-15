[![CI](https://github.com/Sentieon/sentieon-cli/actions/workflows/ci.yml/badge.svg)](https://github.com/Sentieon/sentieon-cli/actions/workflows/ci.yml)

# Sentieon CLI

A command-line interface for the Sentieon software

## Install using pip (recommended)

Install the sentieon-cli into your python environment with `pip`:
```sh
pip install sentieon_cli
```

## Installation with Poetry

Create a new python virtual environment for the project, if needed:
```
# Create a new venv, if needed
python3 -m venv /path/to/new/virtual/environment/sentieon_cli

# Activate the venv
source /path/to/new/virtual/environment/sentieon_cli/bin/activate
```

`sentieon-cli` uses [poetry](https://pypi.org/project/poetry/) for packaging and dependency management. Initially, you will need to install poetry:
```
pip install poetry
```

Clone this repository and cd into the root directory:
```
git clone https://github.com/sentieon/sentieon-cli.git
cd sentieon-cli
```

Use poetry to install the `sentieon-cli` into the virtual environment:
```
poetry install
```

You can then run commands from the virtual environment:
```
sentieon-cli ...
```

## Global arguments
The `sentieon-cli` supports the following global arguments:
- `--verbose`: verbose logging.
- `--debug`: debugging mode for more verbose logging.

## Supported pipelines
- [**DNAscope**](https://support.sentieon.com/docs/sentieon_cli/#dnascope) - DNAscope pipeline implementation for germline SNV and indel calling from short read data.
- [**DNAscope LongRead**](https://support.sentieon.com/docs/sentieon_cli/#dnascope-longread) - DNAscope LongRead pipeline implementations for germline SNV and indel calling from long read data.
- [**DNAscope Hybrid**](https://support.sentieon.com/docs/sentieon_cli/#dnascope-hybrid) - DNAscope short-long-hybrid pipeline.
- [**DNAscope Pangenome**](https://support.sentieon.com/docs/sentieon_cli/#sentieon-pangenome) - DNAscope pangenome alignment and variant calling. Our recommended pipeline for short-read small variant calling.

## License
Unless otherwise indicated, files in this repository are licensed under a BSD 2-Clause License.
