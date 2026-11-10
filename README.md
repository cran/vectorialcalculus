# vectorialcalculus 2.0.0

The package provides 107 public functions for university calculus. It
covers numerical differentiation and integration, geometry and vector
calculus, including the theorems of Green, Stokes and Gauss.

## Upgrading from 1.0.5

Version 2.0.0 introduces a new public interface using the vc_ prefix.
Existing scripts using the 1.0.5 function names need to be updated.
Some arguments and return values also change. See the installed
migration map with:

```r
migration <- system.file("extdata", "legacy-map.csv", package = "vectorialcalculus")
utils::read.csv(migration)
```

For example, spatial curve length now accepts a function that returns
the three coordinates together:

```r
library(vectorialcalculus)
curve <- function(t) c(cos(t), sin(t), t / 2)
velocity <- function(t) c(-sin(t), cos(t), 1 / 2)
vc_curve_length(curve, c(0, 2 * pi), velocity = velocity)$estimate
```

