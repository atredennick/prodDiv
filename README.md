# prodDiv

These scripts reproduce the main results from Fraser et al. 2015 ["Worldwide evidence of a unimodal relationship between productivity and plant species richness"](http://www.sciencemag.org/content/349/6245/302.short) (DOI: 10.1126/science.aab3916). However, to reflect the sampling design, we use generalized linear mixed effects models (GLMMs) to fit linear and quadratic Poisson regressions for the relationship between productivity and species richness.

## Required packages
* `rdryad`
* `lme4`

You can install these using: `install.packages(c("rdryad", "lme4"))`

## License
The MIT License (MIT)

Copyright (c) 2015 Andrew Tredennick

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
