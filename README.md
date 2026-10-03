# rardl

> **This repository is superseded and no longer maintained.**
> At the request of the CRAN team, this package was merged into the CRAN package
> [ardlverse](https://cran.r-project.org/package=ardlverse). The function `rardl()` is maintained there,
> with corrections that are not in this repository. The code here is an older version
> and should not be used for new work.
>
> ```r
> install.packages("ardlverse")
> ```

**Rolling-Window and Recursive ARDL Cointegration Analysis** for R

## Overview

`rardl` provides rolling-window and recursive implementations of the ARDL
bounds cointegration test for detecting structural instability in
cointegrating relationships over time.

## Analysis Types

| `type`        | Description                                        |
|---------------|----------------------------------------------------|
| `"rolling"`   | Rolling-window ARDL bounds test                    |
| `"recursive"` | Recursive (expanding) ARDL bounds test             |
| `"radf"`      | Recursive ADF unit root test                       |
| `"rgranger"`  | Recursive Granger causality test                   |
| `"simulate"`  | Monte Carlo simulation of critical values          |

## Installation

```r
# Old version, GitHub only (see the notice at the top of this page)
# install.packages("remotes")
remotes::install_github("muhammedalkhalaf/rardl")
```

## Usage

```r
library(rardl)

set.seed(42)
n  <- 120
x1 <- cumsum(rnorm(n))
y  <- 0.5 * x1 + rnorm(n, sd = 0.5)
df <- data.frame(y = y, x1 = x1)

# Recursive ARDL
res_rec <- rardl(df, type = "recursive", depvar = "y",
                 indepvars = "x1", init_obs = 40)
print(res_rec)

# Rolling ARDL
res_rol <- rardl(df, type = "rolling", depvar = "y",
                 indepvars = "x1", window_size = 40)
print(res_rol)
```

## References

Shahbaz, M., Khan, A. and Mubarak, M. S. (2023). Energy Economics, 116, 106642.

Khan, A., Shahbaz, M. and Napari, A. (2023). Energy Economics, 117, 106683.

## License

GPL-3
