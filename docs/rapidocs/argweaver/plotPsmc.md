# `plotPsmc`

Plot PSMC output file


## Description

Plot PSMC output file


## Usage

```r
plotPsmc(psmcFile, mu = 1.4e-08, s = 100, sampleAge = 0,
  add = FALSE, timeScale = 1, xlim = NULL, ylim = NULL,
  xlab = "Time", ylab = "Popsize", log = "x", ...)
```


## Arguments

Argument      |Description
------------- |----------------
`psmcFile`     |     file name containing PSMC output
`mu`     |     Mutation rate per base pair per generation
`s`     |     The site compression rate used in PSMC
`sampleAge`     |     (For ancient samples) the sample age, in generations
`add`     |     If TRUE, add to current plotting device
`timeScale`     |     Scale all times by this factor (i.e. to convert generations to years)
`xlim`     |     x-axis range. If NULL, use range in data
`ylim`     |     y-axis range. If NULL, use range in data
`xlab`     |     x-axis label
`ylab`     |     y-axis label
`log`     |     Which axis to use log scale for (either "x", "y", "xy", "n")
`...`     |     Other arguments to be passed to "segments" function (col, lty, lwd, etc)
