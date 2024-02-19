# `subSites`

Subset a sites object


## Description

Subset a sites object


## Usage

```r
subSites(x, keep = NULL, prune = NULL, region = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     a sites object (read in by readSites)
`keep`     |     A list of haplotypes to keep
`prune`     |     A list of haplotypes to prune
`region`     |     A numeric vector of length 2 giving minimum and maximum coordinates (1-based)


## Value

A sites object subsetted by haplotypes and/or region
