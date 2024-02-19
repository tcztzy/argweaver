# `getCols`

Choose columns from summary file


## Description

Given an ARGweaver summary file,
 return columns of interest


## Usage

```r
getCols(x, inds, prefix = "", suffix = "", or = FALSE)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     A data frame read by readSummary()
`inds`     |     A vector of character strings containing individual (or haplotype names) relevant to desired columns
`prefix`     |     Choose only columns with this prefix
`suffix`     |     Choose only columns with this suffix
`or`     |     If TRUE, return columns that mention any of the individuals. If FALSE, return columns that mention all the individuals


## Value

A list of column names from x that meet the criteria


## Note

This may not work as desired if some individual names are subsets of others
