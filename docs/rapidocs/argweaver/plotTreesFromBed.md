# `plotTreesFromBed`

Plot trees from bed file produced by smc2bed


## Description

Plot trees from bed file produced by smc2bed


## Usage

```r
plotTreesFromBed(file = NULL, iter = "max", chrom = NULL,
  start = -1, end = -1, prune = NULL, keepSeqs = NULL,
  col = "black", leafCol = col, leafLabels = NULL, interval = 1,
  timeScale = 1, drawSpr = FALSE, ylab = "Generations", xlab = "",
  logScale = FALSE, ylim = NULL, mar = c(8, 4, 1, 1), add = FALSE,
  regionSide = 1, regionLine = 4, regionRep = TRUE,
  treeInfo = NULL, mod = NULL, popwidth = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     A bed file produced by smc2bed
`iter`     |     The MCMC iteration to use. "max" means use the last iteration in the file. NULL or a value < 0 means to plot all iterations.
`chrom`     |     The chromosome of the region where plots are desired. See note below.
`start`     |     The start coordinate of region where plots are desired. See note below.
`end`     |     The end coordinate of region where plots are desired. See note below.
`prune`     |     A list of leafs to prune from the tree before plotting (Null means prune none)
`keepSeqs`     |     A list of leafs to keep in the tree (NULL means keep all)
`col`     |     Either a character string (single value) giving color of tree. Or, a list giving color to plot each leaf branch. Internal branches will be plotted as "average" color among child leafs.
`leafCol`     |     If given, this can be a list assigning a color to each leaf name (individual names also work if haploid names are in the form indName_1 and indName_2)
`leafLabels`     |     A list indexed by the leaf names in the tree, giving the label to display for each leaf name. If NULL, the leaf names are displayed.
`interval`     |     Only plot trees separated by this interval in base pairs
`timeScale`     |     multiply all branches by this value
`drawSpr`     |     If TRUE, draw the SPR event that turns each tree into the next one
`ylab`     |     label for y axis
`xlab`     |     label for x axis
`logScale`     |     If TRUE, plot in log scale
`ylim`     |     Range for y axis (default; use range of tree)
`mar`     |     Margins for plot
`add`     |     If TRUE, do not create a new plot
`regionSide`     |     Print region of each tree to the margins on this side (1=bottom, 2=left, 3=top, 4=right, anything else = don't print)
`regionLine`     |     Print region of each tree on this line of the margin
`regionRep`     |     If TRUE, include MCMC rep in region string
`treeInfo`     |     If given, should be a character vector of same length as trees
`mod`     |     (Advanced; for use with multi-population version of ARGweaver) A model file read in with the function readPopModel, will draw population model underneath tree
`popwidth`     |     If mod is not NULL, relative widths for each population
`...`     |     Passed to plot function


## Note

If the input file is bgzipp'ed and tabixed (which is done automatically when
 created with the script smc2bed-all), then tabix can be used to read in the file. This
 will be much more efficient if the region chrom:start-end covers a subset of the region
 covered by the entire file.

 This creates a new plot for each tree. If plotting to the screen, probably
 want to call par(ask=TRUE) first.
