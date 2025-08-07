[![CI](https://github.com/Sentieon/sentieon-cli/actions/workflows/ci.yml/badge.svg)](https://github.com/Sentieon/sentieon-cli/actions/workflows/ci.yml)

# Sentieon CLI

A command-line interface for the Sentieon software

## Installation from sdist (recommended)

Download the latest tar.gz file from the GitHub release page, https://github.com/sentieon/sentieon-cli/releases/ and install the package with pip:
```sh
curl -LO https://github.com/Sentieon/sentieon-cli/releases/download/v1.3.0/sentieon_cli-1.3.0.tar.gz
pip install sentieon_cli-1.3.0.tar.gz
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
- [**DNAscope**](docs/dnascope.md) - DNAscope pipeline implementation for germline SNV and indel calling from short read data.
- [**DNAscope LongRead**](docs/dnascope-longread.md) - DNAscope LongRead pipeline implementations for germline SNV and indel calling from long read data.
- [**DNAscope Hybrid**](docs/dnascope-hybrid.md) - DNAscope short-long-hybrid pipeline.

## License
Unless otherwise indicated, files in this repository are licensed under a BSD 2-Clause License.
