# `multiPopsizeFileFromPsmc`

For multiple population model in arg-sample, create popsize file using
 one representative PSMC result file per population


## Description

For multiple population model in arg-sample, create popsize file using
 one representative PSMC result file per population


## Usage

```r
multiPopsizeFileFromPsmc(psmcFiles, ages = NULL, ntimes = 20,
  delta = 0.01, maxTime = 1e+06, mu = 1.4e-08, s = 100,
  outfile = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`psmcFiles`     |     A list with length = number of populations. Each element should be a character string giving the file path of PSMC result for a population, or can be a single numeric value for constant size population (giving diploid population size)
`ages`     |     (Used for ancient samples) A numeric vector of same length as psmcFiles. Each value represents the age in generations of the sample used in PSMC.
`ntimes`     |     Number of time points in arg-sample model
`delta`     |     arg-sample delta parameter
`maxTime`     |     maximum time point used in arg-sample
`mu`     |     mutation rate in mutations per base pair per generation
`s`     |     The bin size used in the fq2psmcfa command; the default in that program is 100
`outfile`     |     If not NULL, path to write output popsize file. This file is suitable for use with the "--popsize-file" option in arg-sample.


## Value

A data frame with columns: population, time, popsizse (invisibly if outfile is
 not NUL)


## Note

Note that this function is for multiple-population model of arg-sample, which is still
 experimental. Most users will want a single vector of population sizes, using the
 function psmcToPopsize function instead.
