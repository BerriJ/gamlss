# %% Install from local source
rm(list = ls())
devtools::load_all()
library(tidyverse)
# devtools::build()
# devtools::install_local(force = TRUE)

# Install from repository
# remotes::install_github("BerriJ/gamlss", ref = "dev")

# Load packages
# library(gamlss)
# %%

# %%
read_csv("data/X.csv", col_names = FALSE) %>%
    as.matrix() -> X

read_csv("data/Y.csv", col_names = FALSE) %>%
    as.matrix() -> Y

# X <- X[1:30, , drop = FALSE]
# Y <- Y[1:30, , drop = FALSE]

devtools::load_all()
mod <- gamlss(
    formula = Y ~ X, #- 1,
    sigma.formula = ~X, #- 1,
    nu.formula = ~X, #- 1,
    family = "TF",
    trace = FALSE
)
head(fitted(mod, what = "nu"))
# %%

# %%
gamlss.dist::dTF(3, mu = 0, sigma = 1, nu = 8)

# %%
