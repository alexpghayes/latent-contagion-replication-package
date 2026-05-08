# Estimating peer effects in noisy, low-rank networks via network smoothing

Peer effect estimation requires precise network measurement, yet most empirical networks are noisy, rendering standard estimators inconsistent. To address measurement error in networks, we propose a method to estimate peer effects in networks whose expected adjacency matrix is low-rank. Our key result is that peer effects over a true unobserved network are asymptotically equivalent to peer effects over the expected adjacency matrix. This result reduces peer effect estimation in noisy networks to low-rank matrix estimation targeting the expected adjacency matrix. We develop our theory for weighted networks observed with additive noise, but the approach can be applied whenever there is a low-rank estimation method suited to the noise structure. We demonstrate via simulations that our approach applies to egocentric samples, aggregated relational data, and networks with missing edges, each requiring a different low-rank estimation method.

## To replicate our computational results

We use [`renv`](https://rstudio.github.io/renv/) to record package dependencies and [`targets`](https://books.ropensci.org/targets/) to coordinate our simulation study and data analysis.

To replicate our results, clone this Github repository. Once you have the repository cloned locally, re-create the project library by calling

``` r
# install.packages("renv")
renv::restore()
```

At this point, you should be ready to replicate our simulation results.

``` r
tar_make()
```

If the pipeline crashes you likely re-run `tar_make()` until it succeeds.

Results will appear as image files in the `figures/` folder.
