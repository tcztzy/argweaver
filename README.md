ARGweaver and ARGweaver-D
=========

*Sampling and manipulating genome-wide ancestral recombination graphs (ARGs).*  

The ARGweaver/ARGweaver-D software package contains programs and libraries for
sampling and manipulating ancestral recombination graphs (ARGs). An ARG
is a rich data structure for representing the ancestry of DNA
sequences undergoing coalescence and recombination. ARGweaver-D is an extension
built into the ARGweaver code that samples ARGs conditional on a demographic
model with population splits and migrations.

*ARGweaver citation:*
[Matthew D. Rasmussen, Melissa J. Hubisz, Ilan Gronau, Adam Siepel. Genome-wide
inference of ancestral recombination graphs. PLoS Genetics 2014.](https://doi.org/10.1371/journal.pgen.1004342)


*ARGweaver-D citation:*
Melissa J. Hubisz, Amy L. Williams, Adam Siepel. Mapping gene flow between ancient
hominins through demography-aware inference of the ancestral recombination graph.
Biorxiv 2019

## Download

ARGweaver can be [downloaded](http://mdrasmus.github.io/argweaver) or 
[forked](https://github.com/CshlSiepelLab/argweaver/) from GitHub.


## Documentation

See the [manual](http://compgen.cshl.edu/ARGweaver/doc/argweaver-manual.html)
for documentation on the programs and file formats associated with ARGweaver.

For ARGweaver-D, see additional documentation here:
[ARGweaver-D manual](http://compgen.cshl.edu/ARGweaver/doc/argweaver-d-manual.html)

## Requirements

The following dependencies must be installed to compile and run
ARGweaver:

- C++ compiler (e.g. [g++](http://gcc.gnu.org))
- [Python](http://python.org)


## Install

ARGweaver can be installed using any of the normal Python mechanisms.
For example, to install from [PyPI](https://pypi.python.org/pypi) you
can use [pip](https://github.com/pypa/pip) with the following command:

```
pip install argweaver
```

Alternatively, ARGweaver can be install using the `setup.py` file:

```
python setup.py install
```

Lastly, ARGweaver can be installed using the `Makefile`:

```
make
```

Once compiled, install the ARGweaver programs (default install in
`/usr`) using:

```
make install
```

By default this will install all files into `/usr`, which may require 
super user permissions.  To specify your own installation path use:

```
make install prefix=$HOME/local
```

If you use this option, make sure `$HOME/local/bin` is in your `PATH` and
`$HOME/local/lib/python2.X/site-packages` is in your `PYTHONPATH`.

ARGweaver can also run directly from the source directory.  Simply add the
`bin/` directory to your `PATH` environment variable or create symlinks to the
scripts within `bin/` to any directory on your `PATH`. Also add the
argweaver source directory to your `PYTHONPATH`. See `examples/` for details.


## Quick Start

Here is a brief example of an ARG simulation and analysis.
To generate simulated data containing a set of DNA sequences and an
ARG describing their ancestry the following command can be used:

```
arg-sim \
    -k 8 -L 100000 \
    -N 10000 -r 1.6e-8 -m 1.8e-8 \
    -o test1/test1
```

This will create an ARG with 8 sequences each 100kb in length evolving in
a population of effective size 10,000 (diploid), with recombination rate
1.6e-8 recombinations/site/generation and mutation rate 1.8e-8 
mutations/generation/site. The output will be stored in the following files:

```
test1/test1.arg   -- an ARG stored in *.arg format
test1/test1.sites -- sequences stored in *.sites format
```

To infer an ARG from the simulated sequences, the following command 
can be used:

```
arg-sample \
    -s test1/test1.sites \
    -N 10000 -r 1.6e-8 -m 1.8e-8 \
    --ntimes 20 --maxtime 200e3 -c 10 -n 100 \
    -o test1/test1.sample/out
```

This will use the sequences in `test1/test1.sites` and it assumes the
same population parameters as the simulation (i.e. `-N 10000 -r 1.6e-8
-m 1.8e-8`).  Also several sampling specific options are given 
(i.e. 20 discretized time steps, a maximum time of 200,000 generations, 
a compression of 10bp for the sequences, and 100 sampling iterations. 
After sampling the following files will be generated:

```
test1/test1.sample/out.log
test1/test1.sample/out.stats
test1/test1.sample/out.0.smc.gz
test1/test1.sample/out.10.smc.gz
test1/test1.sample/out.20.smc.gz
...
test1/test1.sample/out.100.smc.gz
```

The file `out.log` contains a log of the sampling procedure,
`out.stats` contains various ARG statistics (e.g. number of
recombinations, ARG posterior probability, etc), and `out.0.smc.gz`
through `out.100.smc.gz` contain 11 samples of an ARG in *.smc file
format.

To estimate the time to most recent common ancestor (TMRCA) across
these samples, the following command can be used:

```
arg-extract-tmrca test1/test1.sample/out.%d.smc.gz \
    > test1/test1.tmrca.txt
```

This will create a tab-delimited text file containing six columns:
chromosome, start, end, posterior mean TMRCA (generations),
lower 2.5 percentile TMRCA, and upper 97.5 percentile TMRCA. The first
four columns define a track of TMRCA across the genomic region in
BED file format.

Many other statistics can be extracted from sampled ARGs. For more details
see `examples/`.


## Development

The following Python libraries are needed for developing ARGweaver:

```
nose
pyflakes
pep8
```

These dependencies can be installed using

```sh
pip install -r requirements-dev.txt
```

The python tests can be run either with nose or make:
```sh
# Run tests with nose
nosetests test

# Run tests with make
make test
```

There are also C++ tests written using
[googletest](http://code.google.com/p/googletest/), Google's
unit-testing framework.  Googletest can either be installed system-wide or
within the ARGweaver source tree.  For convenience, googletest can be installed
in the source tree using
```sh
make gtest
```

Once installed, c++ unit tests can be run using
```sh
make ctest
```
