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

The Gaussian kernel is used for representing $h$:

$$ 
K(z, z') = \exp\left\{-\sum_{m=1}^{M} r_m (z_m - z'_m)^2\right\} 
$$

Here, $z$ and $z'$ represent vectors of predictors for two different individuals, and $r_m \geq 0$ denotes the tuning parameter that controls the smoothness of $h$ as a function of the exposure $z_m$.
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
library(bkmr)
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

### Fit FBKMR

To fit the FBKMR model, we use the `skmbayes` function which implements data sketching for scalability. Key parameters include:
- `n_subset`: Number of data splits. You can choose the number of splits based on the cpu cores you have. If you work on a single core, the program will run sequentially over the splits. We recommend using 5 splits for a single core.
- `sketch.type`: Type of sketching ("gaussian" or "sub_sampling")
- `iter`: Number of MCMC iterations
- `varsel`: Whether to conduct variable selection

```r
set.seed(111)
fitkm <- skmbayes(y = y, Z = Z, X = X, 
                  n_subset = 5,
                  sketch.type = "gaussian", 
                  iter = 1000, 
                  varsel = TRUE)
# You can compare the fit with the standard bkmr package
# Warning: this will take a long time to run
fitkm_bkmr <- kmbayes(y = y, Z = Z, X = X, 
                  n_subset = 5,
                  sketch.type = "gaussian", 
                  iter = 1000, 
                  varsel = TRUE)
```

