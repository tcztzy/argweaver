# `plotLeafTrace`

Create a leafTrace plot


## Description

The leaftrace draws a horizontal line for every leaf across the ARG,
 such that the distance between adjacent lines reflects the distance between
 the leaves corresponding to those lines. There is an randomness to the
 leaf order and the distance between non-adjacent lines cannot be ascertained


## Usage

```r
plotLeafTrace(file = NULL, x = NULL, col = "black", xlim = NULL,
  ylim = NULL, xlab = "Coordinate", ylab = "", add = FALSE,
  subset = NULL, lwd = 1, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`file`     |     A leaf trace file (produced by the script arg-layout)
`x`     |     (alternative to file) A leaf trace data frame read in by readLeafTrace
`col`     |     Either a character string giving the color of the plot, or a list indexed by leaf or diploid individual names giving the color for each leaf/individual. Will assume that leafs are named with the convention <ind>_1 and <ind>_2 for the two diploid lineages
`xlim`     |     The range for x coordinates
`ylim`     |     The range for y coordinates
`xlab`     |     Label for the x-axis
`ylab`     |     Label for the y-axis
`add`     |     If TRUE, add plot to current plot
`subset`     |     If NULL, plot all leafs. Otherwise only plot individuals/leafs in this character vector
`lwd`     |     Line width
`...`     |     Additional objects passed to the plot function
