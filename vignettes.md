---
title: "Introduction to cDiD with Multiple Time Periods"
author: "Joel Cuerrier"
date:  "2024-12-16"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to cDiD with Multiple Time Periods}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


<!-- rmarkdown::pandoc_version() -->


# Attribution of Vignettes to Original Author
The vignettes included in this library are based on those authored by [Brantly Callaway and Pedro H.C. Sant'Anna] for the 'did' package. We have adapted the content for use in this library to provide users with accessible and informative documentation. The original vignettes authored by [Brantly Callaway and Pedro H.C. Sant'Anna] can be found in the 'did' package documentation [https://bcallaway11.github.io/did/articles/did-basics.html].

# Introduction 
****
* This vignette provides an introductory overview of utilizing Chained-Difference-in-Differences (cDiD) designs to identify and estimate the average impact of engaging in a treatment, with a specific emphasis on leveraging functionalities from the cdid package. 

* The cdid package allows for multiple periods, variation in treatment timing, treatment effect heterogeneity, general missing data patterns, and sample selection on observables.

* The chained DiD rests on a parallel trends assumption, conditioned on being sampled. 

* The chained DiD is robust to some forms of attrition caused by unobservable heterogeneity.

* We identify three main advantages of our approach: (1) it does not require having a balanced panel subsample as it is the case with a standard DiD; (2) identifying assumption allows for sample selection on time-persistent unobservable (or latent) factors, unlike the cross-section DiD which treats the sample as repeated cross-sectional data; and (3) it may also deliver efficiency gains compared to existing methods, notably when the outcome variable is highly time-persistent.

* The cdid package is designed to provide group-time average treatment effects. Future updates to the package will expand its capabilities to include event-study type estimates, which encompass treatment effect parameters corresponding to various durations of exposure to the treatment, as well as overall treatment effect estimates.

# Examples with numerical simulations
****
We propose a simulation design adapted from the first section. Let us specify the
potential outcome as a components of variance:

$$
Y_{it}(D_i) = \alpha_i + \delta_t + \sum_{\tau=2}^{t} \beta_{\tau} D_{i\tau} + \varepsilon_{it}
$$

Where $D_{i\tau} ∈ {0, 1}$ denotes whether individual i has been treated in $\tau$ or earlier. Let us assume that $\tau ∈ {0, ..., T + 1}$ and treatments can only occur in t ≥ 2 so that G ∈ {2, ..., T + 1}. The data generating process is characterized by the following assumptions:

* The individual-specific unobservable heterogeneity is iid gaussian: $\alpha_i ∼ N(1, \alpha^2_{\alpha})$, where $\alpha^2 = 2;$
* The time-specific unobservable heterogeneity is iid gaussian: $\delta_t ∼ N(1, 1)$;
* The treatment effect is constant over time: $\beta_{\tau} = 1$;
* The error term is iid gaussian: $\varepsilon_{it} ∼ N(0, \alpha^2_{\varepsilon})$;
* The probability to receive the treatment at time g, conditional on being treated at g or in the control group, is defined as: 
$$Pr(G_{ig} = 1|X_i, \alpha_i, G_{ig} + C_i = 1) = \frac{1}{1 + \exp(\theta_0 + \theta_1 X_i + \theta_2 \alpha_i \times g)}$$.
Where $X_i ∼ N(1, 1)$ is observable for every i, unlike $\alpha_i$, and $\theta_0 = −1, \theta_1 = 0.4$ and $\theta_2 = 0$ or $\theta_2 = 0.2$. In the latter case, the treatment probability varies with treatment timing and the unobserved individual heterogeneity;
* The sampling probability in the consecutive periods t, t + 1 conditional on $\alpha_i$ is given by:
$$Pr((S_{it} + 1 = 1|\alpha_i)) = \frac{1}{1 + \exp(\lambda_0 + \lambda_1 \alpha_i \times t)}$$
With $\lambda_0 = −1$, and $\lambda_1 = 0$ or $\lambda_1 = 0.2$, so that the sampling process can also vary with time and the unobserved individual heterogeneity.

We simulate the sampled data in two steps. First, we generate a population sample
for each period t to represent individuals that are either treated at t or in the control
group. Second, we sample from this population using the specified process. 

1. Generate a population of individuals
(a) Draw $N = 2 x max $\frac{n}{E_{\alpha}[Pr(S_{itt+1})]}$ individuals per period in order to have $T + 2$ population samples of N individuals, where each individual is characterized by a vector $(\alpha_i, \delta_t, X_i, \epsilon_{it})$.
(b) Separately for each population sample g, draw a uniform random number $(\xi_i \in [0, 1])$ per individual. If $\xi_i \leq Pr(G_{it} = 1|X_i, \alpha_i, G_{it} + C_i = 1)$, then set $(G_{ig} = 1, C_i = 0)$, otherwise set $(G_{ig} = 0, C_i = 1)$.
(c) Compute $Y_{it}$ from $(\alpha_{i}, \delta_{t}, X_{i}, ε_{it}, G_{i0}, ..., _{iT+1}, C_{i})$;

2. Sample from this population
(a) Draw a uniform random number $η_{it} ∈ [0, 1]$ per individual i and period t. If $η_{it} ≤ Pr(S_{itt+1} = 1|\alpha_{i}), then set S_{itt+1} = 1$ and $S_{iττ+1} = 0$ for $τ /neq t$;

(b) Draw (without replacement) n individuals per period t from the population for which $S_{itt+1} = 1$;
(c) Compute the different estimators.
(d) Repeat steps 1(b)-2(c) 1,000 times and report the mean and standard deviation of the estimators.



``` r
# Load the cdid package
remotes::install_github("joelcuerrier/cdid", ref = "main", build_vignettes = FALSE, force = TRUE)
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo joelcuerrier/cdid@main
#> These packages have more recent versions available.
#> It is recommended to update all of them.
#> Which would you like to update?
#> 
#> 1: All                            
#> 2: CRAN packages only             
#> 3: None                           
#> 4: openssl (2.2.2 -> 2.3.0) [CRAN]
#> 
#> openssl (2.2.2 -> 2.3.0) [CRAN]
#> Installing 1 packages: openssl
#> 
#> The downloaded binary packages are in
#> 	/var/folders/18/86jbd8_d11ngf_w3c553vmyr0000gn/T//RtmpmNT903/downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────
#>      checking for file ‘/private/var/folders/18/86jbd8_d11ngf_w3c553vmyr0000gn/T/RtmpmNT903/remotes9d92690dd2a1/joelcuerrier-cdid-03bdacc/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/18/86jbd8_d11ngf_w3c553vmyr0000gn/T/RtmpmNT903/remotes9d92690dd2a1/joelcuerrier-cdid-03bdacc/DESCRIPTION’
#>   ─  preparing ‘cdid’:
#>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>      Removed empty directory ‘cdid/chained previous version/inst’
#>      Omitted ‘LazyData’ from DESCRIPTION
#>   ─  building ‘cdid_0.1.0.tar.gz’
#>      
#> 

library(cdid)
library(dplyr)
set.seed(123)

# 0. We generate data in several steps.
# Simulate data
data=fonction_simu_attrition(theta2_alpha_Gg=0, lambda1_alpha_St=0, sigma_alpha=2, sigma_epsilon=0.5, alpha_percentile=1)

# Filter indivs never observed across all years
data0 <- data[data$P_Y1_longDID == 1,]
data0 <- data0 %>%
  group_by(id) %>%
  filter(!(all(P_Y1_chaine == 0 & annee %in% 1:8))) %>%
  ungroup()

data0 <- data0[,c(1,2,3,6,7,8)]
data0 <- rename(data0,Y = Y1_chaine)
data0 <- rename(data0,date = annee)
data0 <- rename(data0,date_G = annee_G)
data0$P_Y1_chaine <- NULL

# Sort data0 by id and date
data0 <- data0 %>%
  arrange(id, date)

### Correct id's
list_id <-as.data.frame(unique(data0[,"id"]))
list_id$iden_num<-1:dim(list_id)[1] 
data0 <- merge(data0, list_id, by.x = "id", by.y = "id", all.x = TRUE)
data0$id <- data0$iden_num
data0$iden_num <- NULL

# Add a treatment var
data0$treated <- 0
data0$treated[data0$date_G!=0] <- 1

# Add a new variable S
data0 <- data0 %>%
  group_by(id) %>%                                # Group data by id
  mutate(S = ifelse(date == min(date), 1, 0)) %>% # S == 1 for lowest year, 0 otherwise
  ungroup()                                       # Ungroup data

#### CUE: data0 is the typical dataset. The script starts here. The above script is just to make sure we have the expected data form
# Data must have:
# binary treatment var (treated)
# sampling indicator for being observed in date t and t+1 (S)
# outcome variable (Y)
# predictor variable (X): several predictors should work
# id for individuals (id)
# dates (date)
# treatment dates / cohorts (date_G)

# 1. Chained DiD
chained.results =chained(
                    yname="Y",
                    tname="date",
                    idname="id",
                    gname="date_G",
                    xformla=~X, 
                    propensityformla=c("X"), 
                    data=data0,      
                    weightsname="S",  
                    bstrap=FALSE, 
                    biters=1000,
                    debT=3,
                    finT=8,
                    deb=1,
                    fin=8,
                    select='select',
                    treated='treated',
                    cband=TRUE)

# 2. Aggregation
agg.es.chained <- aggte(MP = chained.results, type = 'dynamic')
#> Warning in compute.aggte(MP = MP, type = type, balance_e = balance_e, min_e =
#> min_e, : Used bootstrap procedure to compute simultaneous confidence band

# 3. Print results
agg.es.chained
#> 
#> Call:
#> aggte(MP = chained.results, type = "dynamic")
#> 
#> Reference: Callaway, Brantly and Pedro H.C. Sant'Anna.  "Difference-in-Differences with Multiple Time Periods." Journal of Econometrics, Vol. 225, No. 2, pp. 200-230, 2021. <https://doi.org/10.1016/j.jeconom.2020.12.001>, <https://arxiv.org/abs/1803.09015> 
#> 
#> 
#> Overall summary of ATT's based on event-study/dynamic aggregation:  
#>    ATT    Std. Error     [ 95%  Conf. Int.]  
#>  1.605        0.2396     1.1354      2.0746 *
#> 
#> 
#> Dynamic Effects:
#>  Event time Estimate Std. Error [95% Pointwise  Conf. Band]  
#>          -6   0.1772     0.4876         -1.1263      1.4807  
#>          -5   0.0872     0.3613         -0.8787      1.0531  
#>          -4  -0.1860     0.2767         -0.9257      0.5537  
#>          -3  -0.2290     0.1999         -0.7632      0.3052  
#>          -2  -0.1439     0.1419         -0.5232      0.2355  
#>          -1  -0.2554     0.0892         -0.4937     -0.0171 *
#>           0   1.7951     0.0884          1.5589      2.0313 *
#>           1   1.7615     0.1351          1.4002      2.1227 *
#>           2   1.6541     0.2084          1.0971      2.2110 *
#>           3   1.4875     0.2813          0.7356      2.2395 *
#>           4   1.5304     0.4042          0.4499      2.6109 *
#>           5   1.4014     0.6096         -0.2281      3.0308  
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 
#> Control Group:  Never Treated,  Anticipation Periods:  0
#> Estimation Method:  chained
```



``` r

# 0. We generate data in several steps.
# Simulate data
data=fonction_simu_attrition(theta2_alpha_Gg=0, factor = 0.5, lambda1_alpha_St=0, sigma_alpha=2, sigma_epsilon=0.5, alpha_percentile=1)
#> Error in fonction_simu_attrition(theta2_alpha_Gg = 0, factor = 0.5, lambda1_alpha_St = 0, : unused argument (factor = 0.5)

# data0 <- data %>%
#   arrange(id, annee) %>% # Sort by `id` and `annee`
#   group_by(id) %>% # Group by `id`
#   mutate(previous_P_Y1_chaine = lag(P_Y1_chaine)) %>% # Create a lagged variable
#   filter((P_Y1_chaine== 0 & previous_P_Y1_chaine == 1)|P_Y1_chaine== 1) %>% # Apply your condition
#   ungroup() # Ungroup the data

# Generate the dummy variable "P" with a randomized fraction
data$P_Y1_longDID <- as.numeric(runif(dim(data)[1], min = 0, max = 1)>=0)

# Filter indivs never observed across all years
data0 <- data[data$P_Y1_longDID == 1,]
data0 <- data0 %>%
  group_by(id) %>%
  filter(!(all(P_Y1_longDID == 0 & annee %in% 1:8))) %>%
  ungroup()


data0 <- data0[,c(1,2,3,6,7,8)]
data0 <- rename(data0,Y = Y1_chaine)
data0 <- rename(data0,date = annee)
data0 <- rename(data0,date_G = annee_G)
data0$P_Y1_chaine <- NULL

# Sort data0 by id and date
data0 <- data0 %>%
  arrange(id, date)

### Correct id's
list_id <-as.data.frame(unique(data0[,"id"]))
list_id$iden_num<-1:dim(list_id)[1] 
data0 <- merge(data0, list_id, by.x = "id", by.y = "id", all.x = TRUE)
data0$id <- data0$iden_num
data0$iden_num <- NULL

# Add a treatment var
data0$treated <- 0
data0$treated[data0$date_G!=0] <- 1

# Add a new variable S
data0 <- data0 %>%
  group_by(id) %>%                   # Group data by id
  mutate(S = ifelse(date == min(date), 1, 1)) %>% # S == 1 for lowest year, 0 otherwise
  ungroup()                          # Ungroup data
  
#### CUE: data0 is the typical dataset. The script starts here. The above script is just to make sure we have the expected data form
# Data must have:
# binary treatment var (treated)
# sampling indicator for being observed in date t and t+1 (S)
# outcome variable (Y)
# predictor variable (X): several predictors should work
# id for individuals (id)
# dates (date)
# treatment dates / cohorts (date_G)
  

# 1. GMM DiD
gmm.results = GMM(
                    yname="Y",
                    tname="date",
                    idname="id",
                    gname="date_G",
                    xformla=~X, 
                    propensityformla=c("X"), 
                    data=data0,      
                    weightsname="S",  
                    bstrap=FALSE, 
                    biters=1000,
                    debT=3,
                    finT=8,
                    deb=1,
                    fin=8,
                    select='select',
                    treated='treated',
                    cband=TRUE)
#> List of 14
#>  $ group       : num [1:42] 3 3 3 3 3 3 3 4 4 4 ...
#>  $ t           : num [1:42] 1 3 4 5 6 7 8 1 2 4 ...
#>  $ att         : num [1:42] 0.0531 1.7305 1.5166 1.2282 1.0269 ...
#>  $ V_analytical: num [1:42, 1:42] 7.27 4.05 3.41 4.31 3.81 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:42] "V1" "V2" "V3" "V4" ...
#>   .. ..$ : chr [1:42] "V1" "V2" "V3" "V4" ...
#>  $ se          : num [1:42] 0.0437 0.0494 0.0445 0.0491 0.0448 ...
#>  $ c           : Named num 3.2
#>   ..- attr(*, "names")= chr "95%"
#>  $ inffunc     : num [1:3690, 1:42] 1.216 0.685 -0.214 -0.769 -0.132 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:42] "V1" "V2" "V3" "V4" ...
#>  $ n           : int 3690
#>  $ W           : num [1, 1] 15.6
#>  $ Wpval       : num [1, 1] 0.793
#>  $ aggte       : NULL
#>  $ alp         : num 0.05
#>  $ DIDparams   :List of 29
#>   ..$ yname                       : chr "Y"
#>   ..$ tname                       : chr "date"
#>   ..$ idname                      : chr "id"
#>   ..$ gname                       : chr "date_G"
#>   ..$ xformla                     :Class 'formula'  language ~X
#>   .. .. ..- attr(*, ".Environment")=<environment: 0x12c7aa998> 
#>   ..$ data                        :'data.frame':	29520 obs. of  10 variables:
#>   .. ..$ id           : int [1:29520] 1 1 1 1 1 1 1 1 2 2 ...
#>   .. ..$ date         : int [1:29520] 1 2 3 4 5 6 7 8 1 2 ...
#>   .. ..$ Y            : num [1:29520] 0.739 5.093 2.953 2.337 1.793 ...
#>   .. ..$ X            : num [1:29520] 1.56 1.56 1.56 1.56 1.56 ...
#>   .. ..$ date_G       : num [1:29520] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..$ treated      : num [1:29520] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..$ S            : num [1:29520] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..$ control_group: num [1:29520] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..$ .w           : num [1:29520] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..$ .rowid       : int [1:29520] 1 1 1 1 1 1 1 1 2 2 ...
#>   ..$ control_group               : chr "notyettreated"
#>   ..$ anticipation                : num 0
#>   ..$ weightsname                 : chr "S"
#>   ..$ alp                         : num 0.05
#>   ..$ bstrap                      : logi TRUE
#>   ..$ biters                      : num 1000
#>   ..$ clustervars                 : NULL
#>   ..$ cband                       : logi TRUE
#>   ..$ print_details               : logi FALSE
#>   ..$ pl                          : logi FALSE
#>   ..$ cores                       : num 1
#>   ..$ est_method                  : chr "chained"
#>   ..$ base_period                 : chr "varying"
#>   ..$ panel                       : logi FALSE
#>   ..$ true_repeated_cross_sections: logi FALSE
#>   ..$ n                           : int 3690
#>   ..$ nG                          : int 6
#>   ..$ nT                          : int 8
#>   ..$ tlist                       : int [1:8] 1 2 3 4 5 6 7 8
#>   ..$ glist                       : num [1:6] 3 4 5 6 7 8
#>   ..$ call                        : language GMM(yname = "Y", tname = "date", idname = "id", gname = "date_G", xformla = ~X,      propensityformla = c("X"), d| __truncated__ ...
#>   ..$ delta.att.influ             :List of 4
#>   .. ..$ :'data.frame':	168 obs. of  6 variables:
#>   .. .. ..$ attgt_id: int [1:168] 1 2 3 4 5 6 7 8 9 10 ...
#>   .. .. ..$ nobsG   : num [1:168] 249 249 249 249 249 249 249 249 249 249 ...
#>   .. .. ..$ nobsC   : num [1:168] 3441 3441 3184 2921 2654 ...
#>   .. .. ..$ date_G  : num [1:168] 3 3 3 3 3 3 3 3 3 3 ...
#>   .. .. ..$ date    : int [1:168] 2 3 4 5 6 7 8 3 4 5 ...
#>   .. .. ..$ Y       : num [1:168, 1] -0.0535 1.6762 1.461 1.1699 0.9715 ...
#>   .. ..$ : num [1:3690, 1:2, 1:168] 1 2 3 4 5 6 7 8 9 10 ...
#>   .. ..$ :'data.frame':	3690 obs. of  2 variables:
#>   .. .. ..$ date_G: num [1:3690] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ id    : int [1:3690] 1 2 3 4 5 6 7 8 9 10 ...
#>   .. ..$ :'data.frame':	168 obs. of  48 variables:
#>   .. .. ..$ ATT(3,1): num [1:168] -1 -1 -1 -1 -1 -1 -1 0 0 0 ...
#>   .. .. ..$ ATT(3,2): num [1:168] 1 0 0 0 0 0 0 -1 -1 -1 ...
#>   .. .. ..$ ATT(3,3): num [1:168] 0 1 0 0 0 0 0 1 0 0 ...
#>   .. .. ..$ ATT(3,4): num [1:168] 0 0 1 0 0 0 0 0 1 0 ...
#>   .. .. ..$ ATT(3,5): num [1:168] 0 0 0 1 0 0 0 0 0 1 ...
#>   .. .. ..$ ATT(3,6): num [1:168] 0 0 0 0 1 0 0 0 0 0 ...
#>   .. .. ..$ ATT(3,7): num [1:168] 0 0 0 0 0 1 0 0 0 0 ...
#>   .. .. ..$ ATT(3,8): num [1:168] 0 0 0 0 0 0 1 0 0 0 ...
#>   .. .. ..$ ATT(4,1): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(4,2): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(4,3): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(4,4): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(4,5): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(4,6): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(4,7): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(4,8): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(5,1): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(5,2): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(5,3): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(5,4): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(5,5): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(5,6): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(5,7): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(5,8): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(6,1): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(6,2): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(6,3): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(6,4): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(6,5): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(6,6): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(6,7): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(6,8): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(7,1): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(7,2): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(7,3): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(7,4): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(7,5): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(7,6): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(7,7): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(7,8): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(8,1): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(8,2): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(8,3): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(8,4): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(8,5): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(8,6): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(8,7): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. .. ..$ ATT(8,8): num [1:168] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ att.influ                   :List of 2
#>   .. ..$ :'data.frame':	42 obs. of  3 variables:
#>   .. .. ..$ Y     : num [1:42] 0.0531 1.7305 1.5166 1.2282 1.0269 ...
#>   .. .. ..$ date_G: num [1:42] 3 3 3 3 3 3 3 4 4 4 ...
#>   .. .. ..$ date  : num [1:42] 1 3 4 5 6 7 8 1 2 4 ...
#>   .. ..$ : num [1:3690, 1, 1:42] 1.216 0.685 -0.214 -0.769 -0.132 ...
#>   ..- attr(*, "class")= chr "DIDparams"
#>  $ call        : language GMM(yname = "Y", tname = "date", idname = "id", gname = "date_G", xformla = ~X,      propensityformla = c("X"), d| __truncated__ ...
#>  - attr(*, "class")= chr "MP"
#> NULL

# 2. Aggregation
agg.es.gmm <- aggte(MP = gmm.results, type = 'dynamic')

agg.es.gmm
#> 
#> Call:
#> aggte(MP = gmm.results, type = "dynamic")
#> 
#> Reference: Callaway, Brantly and Pedro H.C. Sant'Anna.  "Difference-in-Differences with Multiple Time Periods." Journal of Econometrics, Vol. 225, No. 2, pp. 200-230, 2021. <https://doi.org/10.1016/j.jeconom.2020.12.001>, <https://arxiv.org/abs/1803.09015> 
#> 
#> 
#> Overall summary of ATT's based on event-study/dynamic aggregation:  
#>     ATT    Std. Error     [ 95%  Conf. Int.]  
#>  1.1312         0.018     1.0959      1.1665 *
#> 
#> 
#> Dynamic Effects:
#>  Event time Estimate Std. Error [95% Simult.  Conf. Band]  
#>          -7  -0.0045     0.0034       -0.0144      0.0053  
#>          -6  -0.0012     0.0023       -0.0078      0.0055  
#>          -5  -0.0027     0.0024       -0.0097      0.0043  
#>          -4  -0.0016     0.0022       -0.0079      0.0047  
#>          -3  -0.0004     0.0023       -0.0069      0.0062  
#>          -2   0.0107     0.0072       -0.0098      0.0312  
#>           0   1.7247     0.0155        1.6803      1.7690 *
#>           1   1.4885     0.0173        1.4389      1.5380 *
#>           2   1.2658     0.0200        1.2085      1.3232 *
#>           3   1.0011     0.0254        0.9285      1.0738 *
#>           4   0.7781     0.0307        0.6903      0.8660 *
#>           5   0.5290     0.0502        0.3853      0.6727 *
#> ---
#> Signif. codes: `*' confidence band does not cover 0
#> 
#> Control Group:  Not Yet Treated,  Anticipation Periods:  0
#> Estimation Method:  chained
```





<!-- source("R/fonction_simu_attrition.R")
source("R/gg.R")
source("R/pre_process_cdid.R")
source("gmm_compute_delta_att.R")
source("gmm_convert_delta_to_att.R")
source("gmm_convert_result.R")
source('process_attgt_gmm.R') -->
