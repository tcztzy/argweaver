# `multiHist`

Overlay several histograms


## Description

x should be a list of data frames with the same structure.
 Will plot a histogram of the column named by "stat"
 for each item in x


## Usage

```r
multiHist(x, stat, col = NULL, xlab = stat, breaks = 50, main = "",
  xlim = NULL, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     list of data frames. Each should have a column named stat
`stat`     |     gives the column to plot. Should be single character string.
`col`     |     A list or character vectors of colors to use for each histogram. If NULL, will use a rainbow color, but each color will be given transparency value alpha=0.5
`xlab`     |     Label for x axis
`breaks`     |     Passed on to hist function. Can be a number of breaks, or a numeric vector of break intervals.
`main`     |     Title for plot
`xlim`     |     limits for x-axis
`...`     |     Passed to "hist" function
