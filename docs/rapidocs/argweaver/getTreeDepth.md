# `getTreeDepth`

Get maximum distance from a node to the leaf


## Description

Get maximum distance from a node to the leaf


## Usage

```r
getTreeDepth(tree, idx = -1)
```


## Arguments

Argument      |Description
------------- |----------------
`tree`     |     A tree, either of type "phylo" (read in by ape), or a character string with a Newick tree
`idx`     |     The index of the root node. Usually this will be -1 and the true root will be used. The value of idx refers to the row of the tree$edge data frame, for trees read by "ape"


## Value

The maximum distance from the indexed node (or root if idx=-1) to any leaf node.
