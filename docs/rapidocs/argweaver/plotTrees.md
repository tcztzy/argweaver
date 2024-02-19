# `plotTrees`

Plot trees from newick string


## Description

Plot trees from newick string


## Usage

```r
plotTrees(trees, prune = NULL, keepSeqs = NULL, treeInfo = NULL,
  col = "black", leafCol = col, leafLabels = NULL, timeScale = 1,
  drawSpr = FALSE, ylab = "Generations", xlab = "",
  logScale = FALSE, ylim = NULL, add = FALSE, mar = c(8, 4, 1, 1),
  regionSide = 1, regionLine = 4, mod = NULL, popwidth = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`trees`     |     A character vector, each element is a newick tree to plot
`prune`     |     A list of leafs to prune from the tree before plotting (Null means prune none)
`keepSeqs`     |     A list of leafs to keep in the tree (NULL means keep all)
`treeInfo`     |     If given, should be a character vector of same length as trees
`col`     |     Either a character string (single value) giving color of tree. Or, a list giving color to plot each leaf branch. Internal branches will be plotted as "average" color among child leafs.
`leafCol`     |     If given, this can be a list assigning a color to each leaf name (individual names also work if haploid names are in the form indName_1 and indName_2)
`leafLabels`     |     A list indexed by the leaf names in the tree, giving the label to display for each leaf name. If NULL, the leaf names are displayed.
`timeScale`     |     multiply all branches by this value
`drawSpr`     |     If TRUE, draw the SPR event that turns each tree into the next one
`ylab`     |     label for y axis
`xlab`     |     Label for x axis
`logScale`     |     If TRUE, plot in log scale
`ylim`     |     Range for y axis (default; use range of tree)
`add`     |     If TRUE, do not create a new plot
`mar`     |     Margins for plot
`regionSide`     |     Print treeInfo to the margins on this side (1=bottom, 2=left, 3=top, 4=right, anything else = don't print)
`regionLine`     |     Print region of eaach tree on this line of the margin
`mod`     |     (Advanced; for use with multi-population version of ARGweaver) A model file read in with the function readPopModel, will draw population model underneath tree
`popwidth`     |     If mod is not null, popwidth can be a numeric vector of length equal to the number of populations, giving relative width of each
`...`     |     Passed to plot function


## Note

This creates a new plot for each tree. If plotting to the screen, probably
 want to call par(ask=TRUE) first.
