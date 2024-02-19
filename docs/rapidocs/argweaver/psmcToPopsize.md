# `psmcToPopsize`

Make a popsize file suitable for arg-sample from PSMC output


## Description

Make a popsize file suitable for arg-sample from PSMC output


## Usage

```r
psmcToPopsize(psmcFile, ntimes = 20, delta = 0.01, maxTime = 2e+05,
  mu = 1.4e-08, s = 100, sampleAge = 0, outfile = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`psmcFile`     |     File with PSMC output
`ntimes`     |     Number of time points in arg-sample model
`delta`     |     arg-sample delta parameter
`maxTime`     |     maximum time point used in arg-sample
`mu`     |     mutation rate in mutations per base pair per generation
`s`     |     The bin size used in the fq2psmcfa command; the default in that program is 100
`sampleAge`     |     The age of the sample, in generations (For ancient samples)
`outfile`     |     If not NULL, path to write output popsize file. This file is suitable for use with the "--popsize-file" option in arg-sample.


## Value

a data frame with time and popsize (invisibly if outfile is provided)
