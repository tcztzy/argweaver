# `plotSites`

Plot a sites file


## Description

Plot a sites file


## Usage

```r
plotSites(x = NULL, file = NULL, start = NULL, end = NULL,
  anc = NULL, singletons = FALSE, allele1Color = "black",
  allele2Color = c("red", "orange", "yellow", "green", "turquoise",
  "purple"), missingColor = rgb(0, 0, 0, 0.1), stretch = 1,
  indOrder = NULL, textColor = NULL, useCoords = TRUE, xlab = NULL,
  textSize = 1, ...)
```


## Arguments

Argument      |Description
------------- |----------------
`x`     |     A sites file object as read by readSites().
`file`     |     (Alternative to x) Sites file (may be gzipped. may NOT be streamed)
`start`     |     The start coordinate of the plot. NULL implies start from beginning of region
`end`     |     The end coordinate of the plot. NULL implies go to end of region.
`anc`     |     Ancestral sample name. If NULL, allele1 is major allele and allele2 is minor allele
`singletons`     |     If TRUE, include singletons in plot
`allele1Color`     |     The color(s) for allele1 (will be recycled across sites)
`allele2Color`     |     The color(s) for allele2 (will be recycled across sites)
`missingColor`     |     The color for N allleles
`stretch`     |     A horizontal expansion factor for sites (so that sparse sites are visible in a large region)
`indOrder`     |     A character vector giving order to list haplotypes (if NULL use order used in sites object)
`textColor`     |     a list saying which color to use for which haplotype label
`useCoords`     |     if FALSE, just plot all variant sites with equal spacing
`xlab`     |     label for x-axis. Default is "Coordinate" if useCoords=TRUE or "Sites index" if useCoords=FALSE
`textSize`     |     the magnification size of haplotype labelscat
`...`     |     Other arguments passed to plot function
