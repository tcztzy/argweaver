# `plotBw`

plot bigWig (a.k.a. bedGraph-style) data


## Description

plot bigWig (a.k.a. bedGraph-style) data


## Usage

```r
plotBw(x, y, y0 = NULL, y1 = NULL, col = "black", alpha = 0.2,
  add = FALSE, xlim = NULL, chrom = NULL, ylab = "Generations",
  xlab = "Coordinate", ylim = NULL, fill = FALSE, fillBottom = 0,
  chromLabel = TRUE, chromBreakLty = 3, xaxt = "s", ...)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     a data frame with at least three columns with names chrom, chromStart, chromEnd
`y`     |     a numeric vector whose length=nrow(x) giving one score per region represented in x
`y0`     |     A lower limit to y (should be numeric vector, same length as x and y)
`y1`     |     An upper limit to y (same format as y0; must be provided if y0 is provided)
`col`     |     The color of the line to draw
`alpha`     |     Transparency factor (0 = transparent, 1 = solid)
`add`     |     If TRUE, add to current plot
`xlim`     |     The x-axis limits
`chrom`     |     If given plot only elements of x for each x$chrom==chrom
`ylab`     |     label for y-axis
`xlab`     |     label for x-axis
`ylim`     |     limits for y-axis. If not given will use range(y)
`fill`     |     if TRUE, fill in area under curve from y to y-coordinate designated by fillBottom
`fillBottom`     |     if fill is TRUE, fill area
`chromLabel`     |     if TRUE and multiple chromosomes are being plotted, print chrom names on x-axis
`chromBreakLty`     |     If multiple chromosomes plotted, this line type used in vertical lines dividing them
`...`     |     Other arguments to be passed to plot functions


## Value

Invisibly returns the x-coordinates of the center of each chromosome plotted

 If y0 and y1 are given, then will shade in regions between y0 and y1 in the same
 color as the line drawn, but made transparent by factor "alpha'
 #TODO: have it take an optional low/high value, as well as optional alpha
