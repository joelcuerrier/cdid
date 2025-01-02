# cdid: The Chained Difference-in-Differences

The `cdid` package extends the popular `did` library by Brantly Callaway to improve efficiency and handle unbalanced panel data in staggered treatment designs. It implements the methodology introduced in:

**Bellego, Benatia, and Dortet-Bernadet (2024), "The Chained Difference-in-Differences", Journal of Econometrics.** <doi:10.1016/j.jeconom.2023.11.002>

## Features

-   **Enhanced Efficiency**: Offers greater precision for balanced panels (smaller standard errors)
-   **Handles Missing Data**: Offers greater precision for unbalanced panels, especially when errors are serially correlated.
-   **Reduced Bias**: Less prone to bias when missingness is tied to unobservable heterogeneity.
-   **Flexible Control Groups**: Supports control groups comprising either "nevertreated" or "notyettreated" units, like in the did library.
-   **Customizable Weighting**: Implements aggregation with two weighting matrix options:
    -   **Identity**: Best for smaller datasets with less overidentifying restrictions (many missings like in rotating panel surveys)
    -   **Two-step**: Best for larger datasets with many overidentifying restrictions (few missings).

Future developments: generalized attrition model (MAR & sequential MAR), doubly-robust estimator (so far only ipw is implemented), and improved computational efficiency. Hopefully, it will be directly integrated within the did library.

## Installation

The `cdid` package can be installed from CRAN using:

``` r
install.packages("cdid")
```

Alternatively, for the development version:

``` r
remotes::install_github("joelcuerrier/cdid", ref = "main", build_vignettes = TRUE, force = TRUE)
```

## Getting Started

### Example Usage with Simulated Data

``` r
library(did) #for comparison
library(cdid)

set.seed(123)

# Generate a balanced dataset with unit fixed-effects
# The true values of the coefficients are based on time-to-treatment. The treatment
# effect is zero before the treatment, 1.75 one period after, 1.5 two period after,
# 1.25 three period after, 1 four period after, 0.75 five period after, 0.5 six  
# period after, etc.

data0 <- fonction_simu_attrition(
  N = 500, T = 8,
  theta2_alpha_Gg = 0.5, lambda1_alpha_St = 0,
  sigma_alpha = 2, sigma_epsilon = 0.1, tprob = 0.5
)

# Ensure all observations are included for a balanced panel
data0$S <- 1

# Run the original `did` library estimation
did.results = did:: att_gt(
  yname="Y",
  tname="date",
  idname = "id",
  gname = "date_G",
  xformla = ~X,
  data = data0,
  weightsname = NULL,
  allow_unbalanced_panel = FALSE,
  panel = TRUE,
  control_group = "notyettreated",
  alp = 0.05,
  bstrap = TRUE,
  cband = TRUE,
  biters = 1000,
  clustervars = NULL,
  est_method = "ipw",
  base_period = "varying",
  print_details = FALSE,
  pl = FALSE,
  cores = 1
)

# Run `cdid` with 2-step weighting matrix
result_2step = att_gt_cdid(yname="Y", tname="date",
                         idname="id",
                         gname="date_G",
                         xformla=~X,
                         data=data0,
                         control_group="notyettreated",
                         alp=0.05,
                         bstrap=TRUE,
                         biters=1000,
                         clustervars=NULL,
                         cband=TRUE,
                         est_method="2-step",
                         base_period="varying",
                         print_details=FALSE,
                         pl=FALSE,
                         cores=1)

# Run `cdid` with identity weighting matrix
result_id = att_gt_cdid(yname="Y", tname="date",
                        idname="id",
                        gname="date_G",
                        xformla=~X,
                        data=data0,
                        control_group="notyettreated",
                        alp=0.05,
                        bstrap=TRUE,
                        biters=1000,
                        clustervars=NULL,
                        cband=TRUE,
                        est_method="Identity",
                        base_period="varying",
                        print_details=FALSE,
                        pl=FALSE,
                        cores=1)

# Print results
print(did.results)
print(result_2step)
print(result_id)
```

### Aggregating Results

After computing the group-time ATT estimates, aggregate results can be obtained using functions from did library

``` r
agg.es.did <- aggte(MP = did.results, type = 'dynamic')
agg.es.2step <- aggte(MP = result_2step, type = 'dynamic')
agg.es.id <- aggte(MP = result_id, type = 'dynamic')

# Print aggregate results
print(agg.es.did)
print(agg.es.2step)
print(agg.es.id)
```

### Working with Unbalanced Panels

The `cdid` library excels with unbalanced panels. Here's an example:

``` r
# Generate a dataset with missing observations based on sampling indicator S
data0 <- fonction_simu_attrition(
  N = 500, T = 8,
  theta2_alpha_Gg = 0.5, lambda1_alpha_St = 0,
  sigma_alpha = 2, sigma_epsilon = 0.1, tprob = 0.5
)

# Keep only non-missing (S==1)
data0 <- data0[data0$S == 1, ]

# Run estimations as before, but specify panel = FALSE for did::att_gt()
did.results = did:: att_gt(
  yname="Y",
  tname="date",
  idname = "id",
  gname = "date_G",
  xformla = ~X,
  data = data0,
  weightsname = NULL,
  allow_unbalanced_panel = FALSE,
  panel = FALSE,
  control_group = "notyettreated",
  alp = 0.05,
  bstrap = TRUE,
  cband = TRUE,
  biters = 1000,
  clustervars = NULL,
  est_method = "ipw",
  base_period = "varying",
  print_details = FALSE,
  pl = FALSE,
  cores = 1
)

#For cdid, there is no difference
result_2step = att_gt_cdid(yname="Y", tname="date",
                         idname="id",
                         gname="date_G",
                         xformla=~X,
                         data=data0,
                         control_group="notyettreated",
                         alp=0.05,
                         bstrap=TRUE,
                         biters=1000,
                         clustervars=NULL,
                         cband=TRUE,
                         est_method="2-step",
                         base_period="varying",
                         print_details=FALSE,
                         pl=FALSE,
                         cores=1)

result_id = att_gt_cdid(yname="Y", tname="date",
                        idname="id",
                        gname="date_G",
                        xformla=~X,
                        data=data0,
                        control_group="notyettreated",
                        alp=0.05,
                        bstrap=TRUE,
                        biters=1000,
                        clustervars=NULL,
                        cband=TRUE,
                        est_method="Identity",
                        base_period="varying",
                        print_details=FALSE,
                        pl=FALSE,
                        cores=1)
```

## Documentation

Complete documentation and detailed examples are available through the package's help pages:

``` r
?cdid
browseVignettes("cdid")
```

and a dedicated webpage: <https://www.davidbenatia.com/projects/cdid-library/>.

## References

Bellego, C., Benatia, D., and Dortet-Bernadet, V. (2024). *The Chained Difference-in-Differences*. Journal of Econometrics. [DOI: 10.1016/j.jeconom.2023.11.002](https://doi.org/10.1016/j.jeconom.2023.11.002)

Callaway, B., & Sant'Anna, P. H. C. (2021). *Difference-in-Differences with Multiple Time Periods*. Journal of Econometrics.

## License

This package is licensed under the GPL-2 license.
