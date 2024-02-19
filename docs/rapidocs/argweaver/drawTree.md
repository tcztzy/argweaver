# `drawTree`

Draw a tree from a newick string; add to current plot


## Description

Draw a tree from a newick string; add to current plot


## Usage

```r
drawTree(tree, col = "black", lwd = 1, timeScale = 1,
  call.plotSpr = FALSE, cex.leafname = 2, leafCol = col,
  leafLabels = NULL, mod = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`tree`     |     A newick string representing a tree
`col`     |     Either a character string (single value) giving color of tree. Or, a list giving color to plot each leaf branch. Internal branches will be plotted as "average" color among child leafs.
`lwd`     |     Line width
`timeScale`     |     multiply all branches by this value
`call.plotSpr`     |     If TRUE and if a SPR event is encoded in the tree (with NHX tags), then draw the SPR event on the tree (currently only works if mod is not NULL)
`cex.leafname`     |     Character size of leaf names
`leafCol`     |     If given, this can be a list assigning a color to each leaf name (individual names also work if haploid names are in the form indName_1 and indName_2)
`leafLabels`     |     A list indexed by the leaf names in the tree, giving the label to display for each leaf name. If NULL, the leaf names are displayed.
`mod`     |     (Advanced; for use with multi-population version of ARGweaver) A model file read in with the function readPopModel, will draw population model underneath tree


## Seealso

[](plotTrees.md) for function which creates a new plot


## Note

The parameters col, leafCol, and leafLabels should all be lists indexed by leaf
 names of the tree. However, they can also be diploid individual names, if the leafs
 are named with the convention indName_1 and indName_2
