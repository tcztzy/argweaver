# ARGweavers

## Pre-requisites

1. Rust toolchains (recommend [rustup](https://rustup.rs/))
2. [CMake](https://cmake.org/) with any modern C++ compiler
3. [rye](https://rye-up.com) or [maturin](https://maturin.rs) for Python binding

## Install

```console
foo@bar:~$ pip install argweaver+git://github.com/tcztzy/argweaver
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
