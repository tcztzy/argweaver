# `readPsmc`

Read a PSMC output file


## Description

Read a PSMC output file


## Usage

```r
readPsmc(psmcFile, mu = 1.4e-08, s = 100, sampleAge = 0)
```


## Arguments

Argument      |Description
------------- |----------------
`psmcFile`     |     file name containing PSMC output
`mu`     |     Mutation rate per base pair per generation
`s`     |     The site compression rate used in PSMC
`sampleAge`     |     (For ancient samples) the sample age, in generations


## Value

A data frame containing fields minTime, maxTime, popsize. MinTime
 and MaxTime are in generations and will be adjusted by sampleAge
