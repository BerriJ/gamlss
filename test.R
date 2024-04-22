# %% Install from local source
rm(list = ls())
devtools::load_all()
library(tidyverse)
# devtools::build()
# devtools::install_local(force = TRUE)

# Install from repository
# remotes::install_github("BerriJ/gamlss", ref = "dev")

# Load packages
library(gamlss)

sigma <- matrix(nrow = sum(data$tag == "train"), ncol = 24)

for (i in 0:23) {
    data <- readr::read_csv(paste0("data/regdata_", i, ".csv"),
        show_col_types = FALSE
    ) %>%
        drop_na()

    sum(data$tag == "train")

    X <- data %>%
        dplyr::filter(tag == "train") %>%
        select(-tag, -Price, -Date)

    Y <- data %>%
        dplyr::filter(tag == "train") %>%
        pull(Price)


    mod <- gamlss(
        formula = Y ~ as.matrix(X),
        sigma.formula = ~ as.matrix(X),
        # data = mtcars[1:4, ],
        # family = "TF",
        trace = FALSE
    )

    sigma[, i + 1] <- fitted(mod, what = "sigma")

    # cat("Hour", i, ":", fitted(mod, what = "nu")[1], "\n")
}

ts.plot(sigma, col = rainbow(24), ylim = c(0, 50))
# %%

# %% Some functions that are used in the gamlss function
source("functions.R")
# %%

# %% equivalent example below
# gamlss <- function(
formula <- formula(Y ~ as.matrix(X))
sigma.formula <- formula(~ as.matrix(X))
nu.formula <- ~1
tau.formula <- ~1
family <- NO()
data <- mtcars[1:4, ]

# Irelevant for us?
# weights <- NULL # for weighted likelihood analysis
# # (not the same as in GLM's)
contrasts <- NULL # one type of contrasts for all  parameters
start.from <- NULL # starting from previous gamlss object
mu.start <- NULL # starting from given values
sigma.start <- NULL

mu.fix <- FALSE # whether the parameter is fixed
sigma.fix <- FALSE
nu.fix <- FALSE
tau.fix <- FALSE
control <- gamlss.control(trace = TRUE)
i.control <- glim.control(glm.trace = TRUE) # the inner circle control (GLIM)

rqres <- gamlss:::rqres

## -----------------------------------------------------------------------------
## =============================================================================
## here is where the proper gamlss function starts
## =============================================================================
## -----------------------------------------------------------------------------
##       Save call for future reference
gamlsscall <- match.call(gamlss, call("gamlss",
    formula = mpg ~ wt + hp,
    family = "NO",
    data = mtcars[1:4, ],
    dataDist = "NO",
    trace = TRUE
)) #   the function call

##       Evaluate the model frame
mnames <- c("", "formula", "data", "weights") #  "subset"  "na.action"
cnames <- names(gamlsscall) # get the names of the arguments of the call
cnames <- cnames[match(mnames, cnames, 0)] # keep only the ones that match with mnames
mcall <- gamlsscall[cnames] # get in mcall all the relevant information but remember
# that the first elenent will be NULL
mcall[[1]] <- as.name("model.frame") # replace NULL with model.frame
##        Specials for smoothing
mcall$formula <- if (missing(data)) {
    terms(formula, specials = .gamlss.sm.list)
} else {
    terms(formula, specials = .gamlss.sm.list, data = data)
}
mu.frame <- eval(mcall, sys.parent()) # evalute the data.frame at the model.frame

## -----------------------------------------------------------------------------------------
## This part deals with the family
family <- as.gamlss.family(family) # bring first the gamlss family
G.dev.expr <- body(family$G.dev.inc) # MS Thursday, April 11, 2002 at 10:34
#  nopar <- family$nopar # the number of parameters for the family
## -----------------------------------------------------------------------------------------
## Now extract the model components using model.extra and model.matrix
## This part deals with the response variable
Y <- model.extract(mu.frame, "response") # extracting the y variable from the formula
if (is.null(dim(Y))) { # if y not matrix
    N <- length(Y)
} else {
    N <- dim(Y)[1]
} # calculate the dimension for y
# .gamlss.bi.list <-  if (exists("gamlss.bi.list",envir=.GlobalEnv))
#                         get("gamlss.bi.list", envir=.GlobalEnv) else .gamlss.bi.list
## extracting now the y and the binomial denominator in case we use BI or BB
if (any(family$family %in% .gamlss.bi.list)) {
    if (NCOL(Y) == 1) {
        y <- if (is.factor(Y)) Y != levels(Y)[1] else Y
        bd <- rep(1, N)
        if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
    } else if (NCOL(Y) == 2) {
        if (any(abs(Y - round(Y)) > 0.001)) {
            warning("non-integer counts in a binomial GAMLSS!")
        }
        bd <- Y[, 1] + Y[, 2]
        y <- Y[, 1]
        if (any(y < 0 | y > bd)) stop("y values must be 0 <= y <= N") # MS Monday, October 17, 2005
    } else {
        stop(paste(
            "For the binomial family, Y must be",
            "a vector of 0 and 1's or a 2 column", "matrix where col 1 is no. successes",
            "and col 2 is no. failures"
        ))
    }
} else if (any(family$family %in% .gamlss.multin.list)) {
    y <- if (is.factor(Y)) {
        unclass(Y)
    } else {
        Y
    }
} else if (is.Surv(Y)) {
    ## checking that the family is censored
    if (length(grep("censored", family$family[[2]])) == 0) {
        stop(paste("the family in not a censored distribution, use cens()"))
    }
    ## checking compatability of Surv object and censored distribution
    if (length(grep(attr(Y, "type"), family$family[[2]])) == 0) {
        stop(paste("the Surv object and the censored distribution are not of the same type"))
    }
    y <- Y
    # if (NCOL(Y) == 2)
    #   {
    #      #.event <- Y[,2]
    #      #    y  <- Y[,1]
    #      y <- Y
    #   }
    # else if (NCOL(Y) == 2)
    # stop("interval censored data are not implemented in gamlss yet")
} else {
    y <- Y
}
## -----------------------------------------------------------------------------------------
## checking the permissible y values
if (!family$y.valid(y)) { # MS Thursday, June 20, 2002 at 16:30
    stop("response variable out of range")
}
## -----------------------------------------------------------------------------------------
## this part is used if start.from is used as argument
## ------------start.from fitted model--------
if (!is.null(start.from)) {
    mu.start <- NULL
    sigma.start <- NULL
    nu.start <- NULL
    tau.start <- NULL
    ##               location model
    if ("mu" %in% start.from$parameters) {
        mu.start <- start.from$mu.fv
    }
    ##               scale-dispersion submodel
    if ("sigma" %in% start.from$parameters) {
        sigma.start <- start.from$sigma.fv
    }
    ##               nu submodel
    if ("nu" %in% start.from$parameters) {
        nu.start <- start.from$nu.fv
    }
    ##               tau submodel
    if ("tau" %in% start.from$parameters) {
        tau.start <- start.from$tau.fv
    }
}
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------
## extract the weights
w <- model.extract(mu.frame, weights) # weights for the likelihood
if (is.null(w)) {
    w <- rep(1, N)
} else if (any(w < 0)) stop("negative weights not allowed") #
#   else if (!all(trunc(w)==w)) warning("weights should be integer values \n",
#         " indicating number of observations with identical values \n") #
## =========================================================================================
## -----------------------------------------------------------------------------------------
##  Set up location-mean submodel:
##             mu.X   design matrix
##        mu.offset   offset in linear predictor
##         mu.start   starting values for mu (optional)
## -----------------------------------------------------------------------------------------
mu.fit <- list() # MS Thursday, January 23, 2003 at 14:46
mu.formula <- formula # ms Wednesday, December 29, 2004
mu.terms <- attr(mu.frame, "terms") #   it peeks up the terms attribute
mu.smoothers <- get.smoothers(mu.terms)
mu.a <- attributes(mu.terms) #  from the model.frame
mu.X <- model.matrix(mu.terms, mu.frame, contrasts) # the mean model matrix
mu.offset <- model.extract(mu.frame, offset) # the mean-location offset
if (is.null(mu.offset)) mu.offset <- rep(0, N)
mu.object <- get.object("mu")
formals(mu.object$dldp, envir = new.env()) <- alist(mu = fv) # this is to get the right GLIM arguments
formals(mu.object$d2ldp2, envir = new.env()) <- alist(mu = fv) #
formals(mu.object$G.di, envir = new.env()) <- alist(mu = fv) #
formals(mu.object$valid, envir = new.env()) <- alist(mu = fv) #
## initial values for mu
if (!is.null(mu.start)) {
    mu <- if (length(mu.start) > 1) mu.start else rep(mu.start, N)
} else {
    (eval(family$mu.initial))
} # MS: Friday, March 29, 2002 at 11:27


## ---------------------------------------------------------------------------------------
##  Set up dispersion-scale submodel:
##           sigma.X   design matrix
##           sigma.offset   offset in linear predictor
##       sigma.start   starting values for sigma (optional)
## ---------------------------------------------------------------------------------------

if ("sigma" %in% names(family$parameters)) {
    orig.Envir <- attr(mcall$formula, ".Environment") # DS fix for Willem Thursday, March 18, 2010
    sigma.fit <- list() # MS Thursday, January 23, 2003 at 14:48
    form.sigma <- other.formula(form = sigma.formula)

    sigma.terms <- terms(form.sigma, specials = .gamlss.sm.list, data = data)

    mcall$formula <- sigma.terms
    attr(mcall$formula, ".Environment") <- orig.Envir # DS fix for Willem Thursday, March 18, 2010
    sigma.frame <- eval(mcall, sys.parent())
    sigma.terms <- attr(sigma.frame, "terms")
    sigma.smoothers <- get.smoothers(sigma.terms)
    sigma.a <- attributes(sigma.terms)
    sigma.X <- model.matrix(sigma.terms, sigma.frame, contrasts)
    sigma.offset <- model.extract(sigma.frame, offset)
    if (is.null(sigma.offset)) sigma.offset <- rep(0, N)
    sigma.object <- get.object("sigma")
    formals(sigma.object$dldp, envir = new.env()) <- alist(sigma = fv) #
    formals(sigma.object$d2ldp2, envir = new.env()) <- alist(sigma = fv) #
    formals(sigma.object$G.di, envir = new.env()) <- alist(sigma = fv) #
    formals(sigma.object$valid, envir = new.env()) <- alist(sigma = fv) #
    formals(family$d2ldmdd, envir = new.env()) <- alist(sigma = sigma) #  ?? I do not think is needed
    ## initial values for sigma
    if (!is.null(sigma.start)) {
        sigma <- if (length(sigma.start) > 1) sigma.start else rep(sigma.start, N)
    } else {
        eval(family$sigma.initial)
    }
}
## -----------------------------------------------------------------------------------------
##  Set up for the 3rd parameter submodel:
##            nu.X   design matrix
##       nu.offset   offset in linear predictor
##        nu.start   starting values for nu (optional)
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------
##  Set up for the 4rd parameter submodel:
##            tau.X   design matrix
##       tau.offset   offset in linear predictor
##        tau.start   starting values for tau (optional)
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------
## =========================================================================================
##  Checking whether proper algorithm  (RS, CG or mixed)
## =========================================================================================
## -----------------------------------------------------------------------------------------
name.method <- substitute(RS())
name.method <- deparse(name.method[1])
list.methods <- c("RS()", "CG()", "mixed()")
i.method <- pmatch(name.method, list.methods, nomatch = 0)
if (!i.method) stop("Method must be RS(), CG() or mixed()")
## -----------------------------------------------------------------------------------------
## fitting the model
fiter <- 0
# conv <- eval(substitute(RS())) # This is where the model is fitted! # TODO

# RS <- function(
n.cyc <- control$n.cyc
no.warn <- TRUE
# ) {

## -end of GLIM.fit------------------------------------------------------------

## ---start RS-----------------------------------------------------------------
## getting the contol papameters
c.crit <- control$c.crit
# n.cyc <- control$n.cyc
trace <- control$trace
autostep <- control$autostep
mu.step <- control$mu.step
sigma.step <- control$sigma.step
nu.step <- control$nu.step
tau.step <- control$tau.step
gd.tol <- control$gd.tol
iter <- control$iter
conv <- FALSE
## initial Gloval deviance
G.dev.incr <- eval(G.dev.expr)
G.dev <- sum(w * G.dev.incr)
G.dev.old <- G.dev + 1
## ----------------------------------------------------------------------------
## the outer iteration starts here
while (abs(G.dev.old - G.dev) > c.crit && iter < n.cyc) {
    # the mean submodel
    if ("mu" %in% names(family$parameters)) {
        if (family$parameter$mu == TRUE & mu.fix == FALSE) {
            # mu.old <- mu
            mu.fit <<- glim.fit(
                f = mu.object, X = mu.X, y = y, w = w,
                fv = mu, os = mu.offset, step = mu.step,
                control = i.control, gd.tol = gd.tol,
                auto = autostep
            )
            mu <<- mu.fit$fv
            mu.object$smooth <- mu.fit$smooth
        }
    }
    # the scale-dispersion submodel
    if ("sigma" %in% names(family$parameters)) {
        print("fitting sigma")
        if (family$parameter$sigma == TRUE & sigma.fix == FALSE) {
            # sigma.old <- sigma


            save(mu, family, sigma.object, sigma.X, y, w, sigma, sigma.offset, sigma.step, i.control, gd.tol, autostep, file = "sigma_glim_fit.RData")

            # Be aware: this function needs mu and family although not stated
            sigma.fit <<- glim.fit(
                f = sigma.object, X = sigma.X, y = y,
                w = w, fv = sigma, os = sigma.offset,
                step = sigma.step, control = i.control,
                gd.tol = gd.tol, auto = autostep
            )
            sigma <<- sigma.fit$fv
            sigma.object$smooth <- sigma.fit$smooth
        }
    }
    # # the nu submodel
    # if ("nu" %in% names(family$parameters)) {
    #     if (family$parameter$nu == TRUE & nu.fix == FALSE) {
    #         # nu.old <- nu
    #         nu.fit <<- glim.fit(
    #             f = nu.object, X = nu.X, y = y,
    #             w = w, fv = nu, os = nu.offset,
    #             step = nu.step, control = i.control,
    #             gd.tol = gd.tol, auto = autostep
    #         )
    #         nu <<- nu.fit$fv
    #         nu.object$smooth <- nu.fit$smooth
    #     }
    # }
    # # the tau submodel
    # if ("tau" %in% names(family$parameters)) {
    #     if (family$parameter$tau == TRUE & tau.fix == FALSE) {
    #         # tau.old <- tau
    #         tau.fit <<- glim.fit(
    #             f = tau.object, X = tau.X, y = y,
    #             w = w, fv = tau, os = tau.offset,
    #             step = tau.step, control = i.control,
    #             gd.tol = gd.tol, auto = autostep
    #         )
    #         tau <<- tau.fit$fv
    #         tau.object$smooth <- tau.fit$smooth
    #     }
    # }
    #   the overall Global Deviance
    G.dev.old <- G.dev
    G.dev.incr <- eval(G.dev.expr)
    G.dev <- sum(w * G.dev.incr)
    iter <- iter + 1
    fiter <<- iter
    if (trace) {
        cat("GAMLSS-RS iteration ", iter, ": Global Deviance = ",
            format(round(G.dev, 4)), " \n",
            sep = ""
        )
    }
    if (G.dev > (G.dev.old + gd.tol) && iter > 1) {
        stop(paste(
            "The global deviance is increasing", "\n",
            "Try different steps for the parameters or the model maybe inappropriate"
        ))
    }
}
if (abs(G.dev.old - G.dev) < c.crit) { # MS Wednesday, June 11, 2003 at 11:58
    # taken out (abs((G.dev-G.dev.old)/(0.1+abs(G.dev.old)))<c.crit&&iter<=n.cyc)
    conv <- TRUE
} else {
    FALSE
}
if (!conv && no.warn) warning("Algorithm RS has not yet converged")
conv
# }

method <- substitute(RS())
# %%

## -----------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------
##  Getting the GAMLSS object out
## ----------------------------------------------------------------------------------------
## first the general output
## calculate the Global deviance again
G.dev.incr <- eval(G.dev.expr)
G.dev <- sum(w * G.dev.incr)
out <- list(
    family = family$family, parameters = names(family$parameters),
    call = gamlsscall, y = y, control = control, weights = w,
    G.deviance = G.dev, N = N, rqres = family$rqres, iter = fiter,
    type = family$type, method = method, contrasts = contrasts
)
# , na.action=na.act
out$converged <- conv
out$residuals <- eval(family$rqres)
noObs <- if (all(trunc(w) == w)) sum(w) else N
out$noObs <- noObs
## binomial denominator
if (any(family$family %in% .gamlss.bi.list)) out$bd <- bd
## -----------------------------------------------------------------------------------------
saveParam <- control$save
##  Output for mean model: ----------------------------------------------------------------
if ("mu" %in% names(family$parameters)) {
    out <- c(out, mu = parameterOut(what = "mu", save = saveParam))
} else {
    out$mu.df <- 0
}
## define now the degrees of freedom for the fit and residuals
out$df.fit <- out$mu.df
out$pen <- out$mu.pen
out$df.residual <- noObs - out$mu.df
## Output for dispersion model: ----------------------------------------------------------
if ("sigma" %in% names(family$parameters)) {
    out <- c(out, sigma = parameterOut(what = "sigma", save = saveParam))
    out$df.fit <- out$mu.df + out$sigma.df
    out$pen <- out$mu.pen + out$sigma.pen
    out$df.residual <- noObs - out$mu.df - out$sigma.df
}
## =======================================================================================
out$P.deviance <- out$G.deviance + out$pen # ms Thursday, May 13, 2004
out$aic <- G.dev + 2 * out$df.fit
out$sbc <- G.dev + log(noObs) * out$df.fit
# MS Thursday, April 22, 2004 at 11:26
#  if ((ls(1,pattern="fiter")=="fiter")) rm(fiter, envir = as.environment(1))
# MS Thursday, January 8, 2004 at 17:52
class(out) <- c("gamlss", "gam", "glm", "lm")
out
# } The end of gamlss function
# %%
