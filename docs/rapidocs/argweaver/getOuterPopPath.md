# `getOuterPopPath`

Get the population path of the parent branch of tree


## Description

Get the population path of the parent branch of tree


## Usage

```r
getOuterPopPath(tree)
```


## Arguments

Argument      |Description
------------- |----------------
`tree`     |     A character vector with a newick tree


## Value

The population path encoded in the newick string for outermost branch


## Note

The population path is usually zero; it is only non-zero for the multiple-population
 ARGweaver model
