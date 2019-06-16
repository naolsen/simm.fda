# simm.fda
# R package simm.fda
Simultaneous inference for misaligned functional data.

Author: Niels Olsen <niels.olsen@math.ku.dk>

Version: 1.0.10


This is an R package for performing analysis on misaligned (multivariate) functional data. Please see simm-fda-short-desc.pdf and https://doi.org/10.1111/rssc.12276. Preprint http://arxiv.org/abs/1606.03295

See also the related package pavpop: http://github.com/larslau/pavpop and take a look at the Bochum arm movement data set https://github.com/larslau/Bochum_movement_data, which we analysed using simm.fda.


Since the first work was done, we have extended the model to include (univariate) functional data, where there the observations are modeled using exponential families. Currently supported are Negative Binomial and Poisson families, but others can easily be defined by the user.
Please see https://doi.org/10.1371/journal.pone.0215914. Preprint: https://arxiv.org/abs/1811.04446.


# Installation 
The recommended way is to use the devtools package, e.g. run `devtools::install_github("naolsen/simm.fda")` from the R interface.
As the package is only based on R code, no special compilers are needed. 
