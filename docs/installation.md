# Installation

## Pre-requisites

1. [rustup](https://rustup.rs/)
2. (Optional) [samtools](http://www.htslib.org/) and
   [bedops](https://bedops.readthedocs.io/en/latest/) for executing `smc2bed-all`

```bash
# Install rustup
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# If you are using macOS, you can install samtools and bedops using homebrew
brew install samtools bedops
# Otherwise, install miniconda for samtools and bedops
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install -c bioconda samtools bedops
# And then, install
pip install argweavers+git://github.com/tcztzy/argweaver
# OR from pypi
pip install argweavers
```

## Build from source

```console
foo@bar:~$ git clone https://github.com/tcztzy/argweaver
foo@bar:~$ cd argweaver
foo@bar:~/argweaver$ cargo build --release
```

## Develop

```console
foo@bar:~/argweaver$ rye sync # Or `maturin dev`
foo@bar:~/argweaver$ pre-commit install
```

## Test

```console
foo@bar:~/argweaver$ cargo test  # test the Rust library
foo@bar:~/argweaver$ pytest      # test the Python binding
```
