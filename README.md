[![CI](https://github.com/brentp/sentieon-cli/actions/workflows/ci.yml/badge.svg)](https://github.com/brentp/sentieon-cli/actions/workflows/ci.yml)

# Sentieon CLI

A command-line interface for the Sentieon software

## Setup

`sentieon-cli` uses [poetry](https://pypi.org/project/poetry/) for packaging and dependency management. Initially, you will need to install poetry into your environment:
```
pip install poetry
```

Clone this repository and cd into the root directory:
```
git clone https://github.com/sentieon/sentieon-cli.git
cd sentieon-cli
```

Use poetry to install the `sentieon-cli` into a virtualenv.
```
poetry install
```

You can then run commands through poetry.
```
poetry run sentieon-cli ...
```

Alternatively, you can use `poetry shell` to activate the environment with the `sentieon-cli`.
```
poetry shell
sentieon-cli ...
```

## Global arguments
The `sentieon-cli` supports the following global arguments:
- `--verbose`: verbose logging.
- `--debug`: debugging mode for more verbose logging.

## Supported pipelines
- [**DNAscope LongRead**](docs/dnascope-longread.md) - DNAscope LongRead pipeline implementations for germline SNV and indel calling from long read data.

## License
Unless otherwise indicated, files in this repository are licensed under a BSD 2-Clause License.