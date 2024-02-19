# `readSites`

Read an ARGweaver sites file


## Description

Read an ARGweaver sites file


## Usage

```r
readSites(file, readProbs = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     Sites file (may be gzipped. may NOT be streamed)
`readProbs`     |     If TRUE, return genotype probabilities as well (if they are in encoded in the file)


## Value

a list with two elements: one describing the region, and the
 other a data frame with a column for each haploid genome and a row
 for each variant site. A third element "probs" will be added if readProbs=TRUE
