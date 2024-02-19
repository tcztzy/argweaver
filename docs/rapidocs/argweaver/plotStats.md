# `plotStats`

Plot statistics trace vs iteration


## Description

Given a data frame read in by readStatsFile, plot one of the
 columns vs iteration number.


## Usage

```r
plotStats(x, stat, add = FALSE, xlim = NULL, ylim = NULL,
  xlab = "Iter", ylab = stat, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     A data.frame read in by readStatsFile. Alternatively, if x is type character, can be the name of the stats file.
`stat`     |     A character string naming the column of x to plot
`add`     |     If FALSE, create a new plot. Otherwise add to current plot
`xlim`     |     The x-axis range (only used when add==FALSE). Taken from the range of data if not given.
`ylim`     |     The y-axis range (only used when add==FALSE). Taken from the range of data if not given
`xlab`     |     The label for x-xais (only used when add=FALSE)
`ylab`     |     The label for y-axis (only used when add=FALSE)
`...`     |     Other options to pass to plot function
