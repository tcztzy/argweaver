# `readArgSummary`

Read arg-sample summary file


## Description

This function will read a file written by arg-summarize


## Usage

```r
readArgSummary(file, burnin = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     A character string, either the name of the stats file or the directory where the stat file can be found. If a directory, there must be only one stat file inside or else an error will occur
`burnin`     |     If not NULL, a single numeric value. Any statistics from iterations < burnin will be discarded. Will be ignored if there is not a column named "MCMC_sample"


## Value

A data frame containing the arg-summarize results
