# Phenotypic forest
`phenotypicForest` is an `R` package that produces two types of plots:

1. `phorest`: similar to epidemiological forest plots but used for comparing the effects of SNPs on phenotypes.  Phenotypes can be grouped by themes. Phorests can also be used to compare groups of SNPs. 
2. `polarHistogram`: plots a large number of histogram on a wheel, to save space.

See the vignette for some examples.

# Installation
You need to install `devtools` for `CRAN` first (`install.packages("devtools")`). Then run the command:
```
devtools::install_github("chrisla/phenotypicForest")
```

# Note
This is an update of an old package, to make it compliant with the current syntax of `ggplot2`. This is really old code that could use a revamp. Be my guest! I'm not planning to spend much time on it.
