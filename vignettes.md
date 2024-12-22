---
title: "Introduction to the Chained Diff-in-Diff (cdid)"
author: "David Benatia"
date:  "2024-12-25"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to the Chained Diff-in-Diff (cdid)}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  markdown: 
    wrap: 72
---

<!-- rmarkdown::pandoc_version() -->

# Introduction

------------------------------------------------------------------------

-   Recent advances in econometrics have focused on the identification,
    estimation, and inference of treatment effects in staggered designs,
    i.e. where treatment timing varies across units. Despite the variety
    of available approaches, none stand out as universally superior, and
    many overlook key challenges such as estimator efficiency or the
    frequent occurrence of unbalanced panel data.

    In our recent paper [(Bellégo, Benatia, Dortet-Bernadet, Journal of
    Econometrics,
    2024)](https://www.sciencedirect.com/science/article/pii/S0304407624001295),
    we address these gaps by extending the well-established method
    developed by [Callaway and Sant'Anna (Journal of Econometrics,
    2021)](https://www.sciencedirect.com/science/article/abs/pii/S0304407620303948)
    The cdid library is therefore intended to be used in connection with
    the **did** library by [Brantly Callaway and Pedro H.C. Sant'Anna]:
    [<https://bcallaway11.github.io/did/articles/did-basics.html>]. Our
    approach improves efficiency and adapts seamlessly to unbalanced
    panels, making it particularly valuable for researchers facing data
    with missing observations.

    This page introduces the **cdid** R library that implements the
    methods from our paper, showcasing their relevance for empirical
    research. **For those short on time, here are the key takeaways:**

    -   **cdid** offers greater precision than **did** in balanced panel
        data.

    -   **cdid** outperforms **did** with unbalanced panels, provided
        errors are serially correlated (e.g., in presence of unit-level
        fixed effects). For time-independent errors, **did** may be
        preferred (see the paper for more details).

    -   **cdid** is less prone to bias in cases of attrition, notably if
        missingness is related to unobservable heterogeneity for example
        due to correlations between individual fixed effects, treatment
        assignment, and attrition/sampling. For cases where missingness
        is related to observables, additional results from our paper
        have yet to be implemented in the library (see the paper).

### Features of the cdid library

**What it supports:**

-   Handles any form of missing data.

-   Allows for control units comprising only "never treated" or also
    "not-yet-treated."

-   Provides two weighting matrix options for GMM-based aggregation of
    $\Delta ATT(g,t)$ into $ATT(g,t)$: identity or 2-step. Small-sample
    simulations show that:

    -   **Identity**: Best for smaller datasets with more missing data,
        and less overidentifying restrictions.

    -   **2-step**: Best for larger, balanced datasets with many
        overidentifying restrictions.

-   Implements simple propensity scores for treatment assignment where
    predictor columns X are treated as constant across time periods. For
    time-varying covariates, users can include additional columns
    (Xt,Xt+1,etc.) in the dataset.

**Current limitations:**

-   The generalized attrition model from our paper, which allows for
    dynamic attrition processes (e.g., past outcomes influencing
    observation status), is not yet implemented.

-   The library currently supports only the IPW estimator, though
    extensions to doubly-robust estimators like in the **did** library
    are planned for the future.

# Getting started

Installing the library requires using remotes at the moment because we
haven't submitted the package on CRAN yet. The command is

``` R
remotes::install_github("joelcuerrier/cdid", ref = "main", build_vignettes = TRUE, force = TRUE)
```

``` r
#Load the relevant packages
library(did) #Callaway and Sant'Anna (2021) 
library(cdid) #Bellego, Benatia, and Dortet-Bernadet (2024)

#Generate a dataset: 500 units, 8 time periods, with unit fixed-effects (alpha)
data0=fonction_simu_attrition(N = 500,T = 8,theta2_alpha_Gg = 0.5, lambda1_alpha_St = 0, sigma_alpha = 2, sigma_epsilon = 0.1, tprob = 0.5)

#We keep all observations, so we have a balanced dataset
data0$S <- 1

#run did
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

#run cdid with 2step weighting matrix
result_2step = att_gt_cdid(yname="Y", tname="date",
                         idname="id",
                         gname="date_G",
                         xformla=~X,
                         data=data0,
                         control_group="notyettreated",
                         weightsname="S",
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

#run cdid with identity weighting matrix
result_id = att_gt_cdid(yname="Y", tname="date",
                        idname="id",
                        gname="date_G",
                        xformla=~X,
                        data=data0,
                        control_group="notyettreated",
                        weightsname="S",
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

Now that we have the ATT(g,t) for all three methods, we can compare
aggregate parameters using the following commands

``` r
# Aggregation
agg.es.did <- aggte(MP = did.results, type = 'dynamic')
agg.es.gmm <- aggte(MP = result_gmm, type = 'dynamic')
agg.es.id <- aggte(MP = result_id, type = 'dynamic')

# Print results
agg.es.did
agg.es.gmm
agg.es.id
```

Let us now do the same but with an unbalanced panel dataset

```{r}
#Generate a dataset: 500 units, 8 time periods, with unit fixed-effects (alpha)
set.seed(123)
data0=fonction_simu_attrition(N = 500,T = 8,theta2_alpha_Gg = 0.5, lambda1_alpha_St = 0, sigma_alpha = 2, sigma_epsilon = 0.1, tprob = 0.5)

#We discard observations based on sampling indicator S
data0 <- data0[data0$S==1,]

#run did
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

#run cdid with 2step weighting matrix
result_2step = att_gt_cdid(yname="Y", tname="date",
                         idname="id",
                         gname="date_G",
                         xformla=~X,
                         data=data0,
                         control_group="notyettreated",
                         weightsname="S",
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

#run cdid with identity weighting matrix
result_id = att_gt_cdid(yname="Y", tname="date",
                        idname="id",
                        gname="date_G",
                        xformla=~X,
                        data=data0,
                        control_group="notyettreated",
                        weightsname="S",
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

```{r}
# Aggregation
agg.es.did <- aggte(MP = did.results, type = 'dynamic')
agg.es.gmm <- aggte(MP = result_gmm, type = 'dynamic')
agg.es.id <- aggte(MP = result_id, type = 'dynamic')

# Print results
agg.es.did
agg.es.gmm
agg.es.id
```

```{=html}
<!-- source("R/fonction_simu_attrition.R")
source("R/gg.R")
source("R/pre_process_cdid.R")
source("gmm_compute_delta_att.R")
source("gmm_convert_delta_to_att.R")
source("gmm_convert_result.R")
source('process_attgt_gmm.R') -->
```
