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
#> ── R CMD build ─────────────────────────────────────────────────────────
#>      checking for file ‘/private/var/folders/18/86jbd8_d11ngf_w3c553vmyr0000gn/T/RtmpbeZJQ6/remotes901a5bdab7a9/joelcuerrier-cdid-6a64b83/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/18/86jbd8_d11ngf_w3c553vmyr0000gn/T/RtmpbeZJQ6/remotes901a5bdab7a9/joelcuerrier-cdid-6a64b83/DESCRIPTION’
#>   ─  preparing ‘cdid’:
#>    checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>      Removed empty directory ‘cdid/chained previous version/inst’
#>   ─  building ‘cdid_0.0.0.9000.tar.gz’
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
#> Error in DIDparams(yname = yname, tname = tname, idname = idname, gname = gname, : could not find function "DIDparams"

# 2. Aggregation
agg.es.chained <- aggte(MP = chained.results, type = 'dynamic')
#> Error in aggte(MP = chained.results, type = "dynamic"): could not find function "aggte"

# 3. Print results
agg.es.chained
#> Error: object 'agg.es.chained' not found
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
#> Error in DIDparams(yname = yname, tname = tname, idname = idname, gname = gname, : could not find function "DIDparams"

# 2. Aggregation
agg.es.gmm <- aggte(MP = gmm.results, type = 'dynamic')
#> Error in aggte(MP = gmm.results, type = "dynamic"): could not find function "aggte"

agg.es.gmm
#> Error: object 'agg.es.gmm' not found
```





<!-- source("R/fonction_simu_attrition.R")
source("R/gg.R")
source("R/pre_process_cdid.R")
source("gmm_compute_delta_att.R")
source("gmm_convert_delta_to_att.R")
source("gmm_convert_result.R")
source('process_attgt_gmm.R') -->
