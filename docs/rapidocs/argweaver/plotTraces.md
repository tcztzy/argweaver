# `plotTraces`

Plot ARGweaver statistics vs iteration


## Description

Plots the statistics in stat files output by arg-sample vs iteration.
 By default plots all statistics given in the stat file (as well
 as any extra created by readStatFiles).


## Usage

```r
plotTraces(files = NULL, dir = NULL, outroot = NULL, stats = NULL,
  colors = if (length(files) == 1) {     "black" } else {
  rainbow(length(files)) }, xlim = NULL, callPar = TRUE, ncol = 2,
  addLegend = (length(files) != 1), ...)
```


## Arguments

Argument      |Description
------------- |----------------
`files`     |     A character vector of file names giving the stat files to plot
`dir`     |     (alternative to files; ignored if files is given). Plot all stat files found in this directory.
`outroot`     |     (only used when dir is provided). Only match stat files named <outroot>.stat
`stats`     |     A list of statistics to plot (i.e., "likelihood", "recombs", "prior", etc). If NULL, then plot all statistics found in the stat files
`colors`     |     A list of colors, one per file, to use in the plots
`xlim`     |     The range of iterations to plot. If NULL, use entire range
`callPar`     |     If TRUE, call the par() function to create a new graphics object with a plot for each statistic on the same page
`ncol`     |     (Only used if callPar=TRUE), arrange the plots in this many columns
`addLegend`     |     If TRUE, add an additional plot that contains a legend mapping files to colors
`...`     |     Arguments passed to plotStats()


## Value

invisibly returns a data.frame detailing the files and colors used
