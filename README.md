# circtools
  
![CI](https://github.com/ATpoint/circtools/actions/workflows/tests.yml/badge.svg)
  
A little package with functions related to the analysis of circadian / circular data in R. 
Note that it is mainly intended for my own work, but others are free to use it if you find it useful.

It includes functions to detect rhythmicity in circadian data using cosinor regression as well as
adaptations of basic data analysis functions such as intersection, distance and density calculation but modified to
respect the circular / periodic nature of the input data. For details, please see the function help at `?function_name`.

Functions:  
  
- `run_cosinor()` detects circadian rhythmicity using a cosinor model via the [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) framework.
It scales efficiently to almost any realistic number of samples and observations.  
- `circular_density()` calculates densities of a vector of circular/periodic values, such as acrophases, respecting the circular nature of the data, thereby avoiding skewed density estimates at the beginning and end of the data range  
- `circular_distance()` calculates directional circular distance between two timepoints  
- `circular_distance_reference()` calculates absolute circular distance between a query and one or many reference timepoints  
- `circular_intersect()` intersects a timepoint with a time window respecting periodicity of data  
- `make_data()` simulates circadian data based on a cosinor model with adjustable mesor, amplitude and acrophase  
- `make_design()` creates a design matrix based on the cosinor model for a provided set of timepoints  
- `rad2period()` converts acrophase values from radians hours to period hours
- `run_model_selection()` classifies circadian data from two conditions into differential categories using a model fitting approach, adopted from the `compareRhythms` package

## Installation

```r
# Install Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install()

# Install limma and circtools
BiocManager::install(c("limma", "atpoint/circtools"))
```
