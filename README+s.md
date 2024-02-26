# ARGweavers

## Pre-requisites

1. [rustup](https://rustup.rs/)
2. [CMake](https://cmake.org/) with any modern C++ compiler
3. [rye](https://rye-up.com) or [maturin](https://maturin.rs) for Python binding
4. [samtools](http://www.htslib.org/) for `smc2bed-all`

## Install

```console
foo@bar:~$ pip install argweaver+git://github.com/tcztzy/argweaver
```

After install the Python package, you don't need any additional steps to use the
binaries, all the binaries are installed to the Python environment's `bin`
directory.

This repository provides Python wrappers for the following programs:

- [arg-likelihood](./src/arg-likelihood.cpp)
- [arg-sample](./src/arg-sample.cpp)
- [arg-summarize](./src/arg-summarize.cpp)
- [smc2bed](./src/smc2bed.cpp)

And rewrite of the following programs:

- [smc2bed-all](./argweavers/scripts/smc2bed_all.py)

## Usage

```console
foo@bar:~$ smc2bed-all --help
usage: smc2bed-all [-h] [-s <startnum>] [-e <endnum>] [-i <interval>]
                   [-r <region>]
                   <baseout>

Convert all *.smc files from a single arg-sample run into a bed.gz file
which can be parsed by arg-summarize

positional arguments:
  <baseout>      base name of arg-sample output (specified with arg-sample
                 -o). Creates a file called <baseout>.bed.gz.

optional arguments:
  -h, --help     show this help message and exit
  -s <startnum>  start with this MCMC rep of arg-sample (default: 0)
  -e <endnum>    end with this MCMC rep (default: stop when no file exists)
  -i <interval>  sampling interval (will be auto-detected from output if not
                 specified)
  -r <region>    region (in format START-END, 1-based, inclusive) to pass to
                 smc2bed (default: run on all coordinates)
```

## Build

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
foo@bar:~/argweaver$ cargo test                # test the Rust library
foo@bar:~/argweaver$ rye sync --features test  # install pytest
foo@bar:~/argweaver$ pytest                    # test the Python binding
```
