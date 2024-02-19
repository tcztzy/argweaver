# `readStatsFile`

Read arg-sample stats file


## Description

This function will read a <outroot>.stats file written by arg-sample.


## Usage

```r
readStatsFile(file, burnin = NULL, thinVal = 1)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     A character string, either the name of the stats file or the directory where the stat file can be found. If a directory, there must be only one stat file inside or else an error will occur
`burnin`     |     If not NULL, a single numeric value. Any statistics from iterations < burnin will be discarded
`thinVal`     |     Remove all iters for which iteration mod thinVal != 0


## Value

A data frame containing the data from the stats file. Only
 MCMC iterations of type "resample" will be retained. In addition to
 the original columns from the stats file, several statistics may be
 added to the data frame as follows:
 "joint2" is the sum of prior2 and likelihood and is an alternative
 posterior likelihood statistic that is preferred over "joint",
 which is the sum of prior and likelihood. Whereas "prior" is the
 prior probability of the sampled ARG under the model, "prior2"
 integrates over all possible recombination events that could produce
 the same series of local trees.
 "total_recomb" is the sum of "recombs" and "invis_recombs"
