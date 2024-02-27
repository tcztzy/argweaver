# Installation

## Pre-requisites

1. [rustup](https://rustup.rs/)
2. (Optional) [samtools](http://www.htslib.org/) for executing `smc2bed-all`
3. (Optional) R for plotting

```bash
# Install rustup
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# If you are using macOS, you can install samtools using homebrew
brew install samtools
# Otherwise, install miniconda for samtools
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install -c bioconda samtools
# Install R
sudo apt-get install r-base
# Or
brew install r
# Restore the environment
R -e "renv::restore()"
R -e "renv::install('./R/argweaver')"
# And then, install
pip install argweavers+git://github.com/tcztzy/argweavers
# OR from pypi in the future
pip install argweavers
```

## Notions for R package

### `ape`
If you are using macOS with Apple Silicon, you may also face the following error while install `ape`:

```
fatal error: pcre.h: No such file or directory
```

You can fix it by installing `pcre` using homebrew:

```bash
brew install pcre
```

And then export `CFLAGS`:
```bash
export CFLAGS="-I/opt/homebrew/include"
```

### `RPHAST`

Since the RPHAST package is not available on CRAN, you need to install it from
GitHub. Since the original RPHAST package is not compatible with the latest
version of R, you need to install the coregenomics/RPHAST@bugfix/f77_call branch from GitHub. [See here for the details](https://github.com/CshlSiepelLab/RPHAST/pull/6).

```R
devtools::install_github("coregenomics/RPHAST@bugfix/f77_call")
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
