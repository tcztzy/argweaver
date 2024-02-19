# `getTimes`

Compute log-distributed discrete times used by arg-sample


## Description

This function returns a vector of times which are spaced
 on a log scale according to the parameter delta.


## Usage

```r
getTimes(n = 20, maxTime = 2e+05, delta = 0.01)
```


## Arguments

Argument      |Description
------------- |----------------
`n`     |     The number of time points between 0 and maxTime, inclusive
`maxTime`     |     The maximum time
`delta`     |     Spacing parameter. Times get denser towards t=0 as delta gets small


## Value

A vector of times of length n
