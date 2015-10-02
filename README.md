# prodDiv

These scripts reproduce the main results from Fraser et al. 2015 ["Worldwide evidence of a unimodal relationship between productivity and plant species richness"](http://www.sciencemag.org/content/349/6245/302.short) (DOI: 10.1126/science.aab3916). However, to reflect the sampling design, we use generalized linear mixed effects models (GLMMs) to fit linear and quadratic Poisson regressions for the relationship between productivity and species richness.

## Required packages
* `rdryad`
* `lme4`

You can install these using: `install.packages(c("rdryad", "lme4"))`
