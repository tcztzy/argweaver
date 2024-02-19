# `writeSites`

Write an ARGweaver sites file


## Description

Write an ARGweaver sites file


## Usage

```r
writeSites(sites, filename)
```


## Arguments

Argument      |Description
------------- |----------------
`sites`     |     Sites object
`file`     |     Sites file (may be gzipped. may NOT be streamed)


## Value

a list with two elements: one describing the region, and the
 other a data frame with a column for each haploid genome and a row
 for each variant site. A third element "probs" will be added if readProbs=TRUE
