# `plotTreesAtSites`

call plotTreeFromBed at the sites in a site file; color individual nodes by allele


## Description

call plotTreeFromBed at the sites in a site file; color individual nodes by allele


## Usage

```r
plotTreesAtSites(bedFile, sites, pos = NULL, start = NULL,
  end = NULL, doSingletons = FALSE, colors = list(A = "red", C =
  "green", G = "turquoise", T = "purple", N = "gray"), ...)
```


## Arguments

Argument      |Description
------------- |----------------
`bedFile`     |     bed.gz file from smc2bed with trees
`sites`     |     Sites object, as returned by readSites(). Will plot one tree for each site.
`pos`     |     If not NULL, a vector of integers giving which sites to plot. If NULL, plot all sites
`start`     |     start position
`end`     |     end position
`doSingletons`     |     whether to plot singletons
`colors`     |     Color to use for each allele in leaf of tree
`...`     |     Additional arguments passed to plotTreesFromBed
