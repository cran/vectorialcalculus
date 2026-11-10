# vectorialcalculus 2.0.0

## Breaking changes

* The public interface uses the vc_ prefix and replaces the 47 export names
  from 1.0.5. Existing calls require migration. The installed legacy-map.csv
  lists recovery destinations and key differences.
* Calculation and visualization have separate interfaces. Plotting is
  optional and uses plotly.

## Features

* Numerical differentiation, integration, and multivariable calculus with
  diagnostic results.
* Spatial curves, vector fields, line and surface integrals, and
  numerical comparisons for Green, Stokes and Gauss.
