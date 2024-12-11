# Fast Bayesian Kernel Machine Regression Package

In this document, we illustrate the main features of the `fbkmr` R package through examples. Additional information on the statistical methodology and computational details are provided in the accompanying documentation and research articles.

## How to Cite

The package applies methods introduced in the paper:

**Sonabend, A., Zhang, J., Schwartz, J., Coull, B.A. and Lu, J., 2024. Scalable Gaussian process regression via median posterior inference for estimating multi-pollutant mixture health effects. arXiv preprint arXiv:2411.10858.**

## Brief Overview of Fast Bayesian Kernel Machine Regression Method

Kernel machine regression (KMR), also known as Gaussian process regression, is a popular tool in the machine learning literature. For large datasets, we implement data sketching techniques to make computation feasible. The general modeling framework is:

$$
 g(\mu_i) = h(z_{i1}, \ldots, z_{iM}) + \beta x_i, \quad i = 1, \ldots, n
$$

where $g$ is a monotonic link function, $\mu_i = E(Y_i)$, $h$ is a flexible function of the predictor variables $z_{i1}, \ldots, z_{iM}$, and $x$ is a vector of covariates assumed to have a linear relationship with the outcome ($\beta$ is the corresponding vector of coefficients).

The R package `bmkr` can be used to fit the BKMR model. However, when the data size is large, the computation can be slow.

Our **fast Bayesian kernel machine regression** (FBKMR) method addresses scalability issues for massive dataset by employing a divide-and-conquer strategy. This approach enables efficient computation for large datasets while maintaining robust inference capabilities.

The proposed framework involves splitting the data into smaller subsets, computing posterior distributions independently within these subsets using a rescaled Gaussian kernel, and then aggregating these results via a generalized geometric median. This methodology ensures the posterior distribution reflects the full data while mitigating computational complexity. 

## Installation

This package requires the Gurobi optimization solver. Please follow the instructions below to install Gurobi in R.


1. **Install the Gurobi Optimizer**: Download and install the appropriate version of the Gurobi Optimizer for your operating system from the [Gurobi website](https://www.gurobi.com/downloads/).

2. **Obtain a Gurobi License**: Gurobi requires a valid license. Academic users can obtain a free license by registering on the [Gurobi website](https://www.gurobi.com/academia/academic-program-and-licenses/). After registration, follow the provided [instructions](https://www.gurobi.com/features/academic-named-user-license/) to set up your license. Whence you run grbgetkey using the argument provided on the Academic License Detail page (ex: grbgetkey ae36ac20-16e6-acd2-f242-4da6e765fa0a). The grbgetkey program will prompt you to store the license file on your machine. Open terminal and run the following command to set up the license (change the path to the actual path to your license file):

    ```bash
    export GRB_LICENSE_FILE="/path/to/gurobi.lic"
    ```


3. **Locate the Gurobi R Package File**: The Gurobi R package (`gurobi`) is included with the Gurobi installation. You can find the package file in the `R` subdirectory of your Gurobi installation directory. For example, if you installed Gurobi 12.0.0, the default installation paths could be:

   - **Windows**: `C:\gurobi1200\win64\R\gurobi_12.0-0.zip`
   - **macOS**: `/Library/gurobi1200/macos_universal2/R/gurobi_12.0-0_R_4.4.1.tgz`
   - **Linux**: `/opt/gurobi1200/linux64/R/gurobi_12.0-0_R_4.4.1.tar.gz`

The actual path could be different depending on your Gurobi version and installation directory.
   

4. **Install the 'slam' Package in R**: The Gurobi R package depends on the `slam` package. Install it from CRAN by running the following command in your R console:

   ```R
   install.packages('slam')
   ```

5. **Install the Gurobi R Package**: Use the `install.packages` function to install the Gurobi R package from the local file. Replace `<path_to_gurobi_package>` with the actual path to the Gurobi R package file located in step 3.

   ```R
   install.packages('<path_to_gurobi_package>', repos = NULL)
   ```

   For example,  on macOS:

   ```R
   install.packages('/Library/gurobi1200/macos_universal2/R/gurobi_12.0-0_R_4.4.1.tgz', repos = NULL)
   ```

   On Linux:

   ```R
   install.packages('/opt/gurobi1200/linux64/R/gurobi_12.0-0_R_4.4.1.tar.gz', repos = NULL)
   ```

   

6. **Verify the Installation**: Load the Gurobi package in R to ensure it was installed correctly:

   ```R
   library(gurobi)
   result <- gurobi(list(A = matrix(c(1, 1), nrow = 1), obj = c(1, 1), modelsense = 'max', rhs = 1, sense = '<'))
   ```

   If the package loads without errors, the installation was successful.

## Example of Using `fbkmr` Package
First, load the R package.

```r
library(fbkmr)
```

To illustrate the main features of the R package `fbkmr`, let's first generate some data. We have built in a few functions directly into the R package for this purpose.

```r
set.seed(111)
dat <- SimData(n = 1000, M = 4) # Using larger sample size for sketching example
y <- dat$y
Z <- dat$Z
X <- dat$X
```

Rember when analyzing real datasets, we need to standardize the Z (exposure variables). To prevent potential sigularity problems, it would be great also standardize the covariate X, as long as it make sense. Too much categorical covariates may cause the program convergence problematicly, thus, it would be better using fewer categorical covariates or just combine categories in a reasonable way.

When analyzing real datasets, X here should be a model matrix, which can be generated in a similar way doing linear regression:
```r
data <- data.frame(y = y, x = X)
x=model.matrix(y ~ x,data=data)
```
Let’s view the true exposure-response function used to generate the data
```r
z1 <- seq(min(dat$Z[, 1]), max(dat$Z[, 1]), length = 20)
z2 <- seq(min(dat$Z[, 2]), max(dat$Z[, 2]), length = 20)
hgrid.true <- outer(z1, z2, function(x,y) apply(cbind(x,y), 1, dat$HFun))

res <- persp(z1, z2, hgrid.true, theta = 30, phi = 20, expand = 0.5, 
             col = "lightblue", xlab = "", ylab = "", zlab = "")
```
![True exposure-response function](https://github.com/junwei-lu/fbkmr/blob/main/figs/surf.png)

### Fit FBKMR

To fit the FBKMR model, we use the `skmbayes` function which implements data sketching for scalability. 
```r
set.seed(111)
fitkm <- skmbayes(y = y, Z = Z, X = X, iter = 1000, n_subset=5, file_path = NULL,save_loc = FALSE, linprog=FALSE, n.cores = 4)
```
Key parameters include:
- `n_subset`: Number of data splits. You can choose the number of splits based on the cpu cores you have. If you work on a single core, the program will run sequentially over the splits. We recommend using 5 splits for a single core.
- `sketch.type`: Type of sketching ("gaussian" or "sub_sampling")
- `iter`: Number of MCMC iterations
- `varsel`: Whether to conduct variable selection
- `iter`: Number of iterations of the MCMC sampler
- `y`: Vector of outcomes 
- `Z`: Matrix of exposures (each column is an exposure variable)
- `X`: Matrix of covariates (each column is a covariate)
- `n_subset`: Number of sub datasets for the algorithm to run (n_subset=1 for original BKMR)
- `file_path`: Location to save intermediate results for future use
- `save_loc`: Whether to save intermediate data (TRUE/FALSE)
- `linprog`: Whether to find overall coefficient estimates for covariates X (TRUE/FALSE)
- `n.cores`: Number of cores to use for parallel computing (defaults to half of total cores if not specified)

When `linprog=TRUE`, the algorithm returns an overall bkmrfit object similar to the `bkmr` package. When `FALSE`, it returns a list of n_subset bkmrfit objects. Since linear programming can be time-consuming with large n_subset, we recommend setting `linprog=FALSE` if you are not interested in exploring covariate effects.

For choosing the suitable number of subsets, we recommend that for each sub dataset, the sample size is around 1000. However, the number of subsets should not greater than square root of total sample size.

### Summarize model output

**Plot the predictor-response function**

One cross section of interest is the univariate relationship between each zm and the outcome, where all of the other exposures are fixed to a particular percentile. This can be done using the function PredictorResponseUnivar. The argument specifying the quantile at which to fix the other exposures is given by q.fixed (the default value is q.fixed = 0.5).

```r
pred.resp.univar <- PredictorResponseUnivar(fit = fitkm, n_subset = 5, q.fixed = 0.5, n.cores = 4)
library(ggplot2)
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
    geom_smooth(stat = "identity") + 
    facet_wrap(~ variable) +
  ylab("h(z)")
```
![Predictor-response function](https://github.com/junwei-lu/fbkmr/blob/main/figs/01_hz_plot.png)

Building upon the previous example, we can similarly visualze the bivarate exposure-response function for two predictors, where all of the other predictors are fixed at a particular percentile.

```r
pred.resp.bivar <- PredictorResponseBivar(fit = fitkm, q.fixed = 0.5,n_subset = 5, n.cores = 4)
ggplot(pred.resp.bivar, aes(z1, z2, fill = est)) + 
    geom_raster() + 
    facet_grid(variable2 ~ variable1) +
  scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  xlab("expos1") +
  ylab("expos2") +
  ggtitle("h(expos1, expos2)")
```
![Bivariate predictor-response function](https://github.com/junwei-lu/fbkmr/blob/main/figs/02_expose.png)

To investigate the predictor-response function of a single predictor in `Z` for the second predictor in `Z` fixed at various quantiles, we can use the `PredictorResponseBivarLevels` function, which takes as input the bivariate exposure-response function outputted from the previous command, where the argument `qs` specifies a sequence of quantiles at which to fix the second predictor.

```r
pred.resp.bivar.levels <- PredictorResponseBivarLevels(
  pred.resp.df = pred.resp.bivar, Z=Z,qs = c(0.1, 0.5, 0.9))
ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
    geom_smooth(aes(col = quantile), stat = "identity") + 
    facet_grid(variable2 ~ variable1) +
    ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")
```
![Predictor-response function](https://github.com/junwei-lu/fbkmr/blob/main/figs/03_est.png)
**Plot the predictor-response function**

In addition to visually inspecting the estimated predictor-response function $h$, one may also wish to calculate a range of summary statistics that highlight specific features of the (potentially) high-dimensional surface. One potential summary measure of interest is to compute the overall effect of the predictors, by comparing the value of $h$ when all predictors are at a particular percentile as compared to when all of them are at their 50th percentile. The function `OverallRiskSummaries_wasp` allows one to specify a sequence of values of quantiles using the argument `qs` and the fixed quantile (the default is the 50th percentile) using the argument `q.fixed`.

```r
risks.overall <- OverallRiskSummaries_wasp(fit = fitkm,  
                                      qs = seq(0.25, 0.75, by = 0.05), 
                                      q.fixed = 0.1,
                                      n_subset = 5,
                                      n.cores = 4)
```
We can also plot the overall risk summaries; here we use the `ggplot` package.

```r
ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
    geom_pointrange()
```
![Overall risk summaries](https://github.com/junwei-lu/fbkmr/blob/main/figs/04_risk.png)
Another summary of h that may be of interest would be to summarize the contribution of an individual predictor to the response. For example, we may wish to compare risk when a single predictor in h is at the 75th percentile as compared to when that predictor is at its 25th percentile, where we fix all the remaining predictors to a particular percentile. We refer to this as the single-predictor health risks, and these can be computed using the function `SingVarRiskSummaries_wasp`. The two different quantiles at which to compare the risk are specified using the `qs.diff` argument, and a sequence of values at which to fix the remaining pollutants can be specified using the `q.fixed` argument.

```r
risks.singvar <- SingVarRiskSummaries_wasp(fit = fitkm, y = y, Z = Z, X = X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      n_subset = 5,
                                      n.cores = 4)
ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
    geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
```
![Single-predictor health risks](https://github.com/junwei-lu/fbkmr/blob/main/figs/05_qfix.png)



