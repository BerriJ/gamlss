rm(list = ls())
load("sigma_glim_fit.RData")

sigma <- rep(sd(y), length(y))

f <- sigma.object
X <- sigma.X
y <- y
w <- w
fv <- sigma
os <- sigma.offset
step <- sigma.step # This is one!
control <- i.control
gd.tol <- gd.tol
auto <- autostep

# Below it the glm.fit function
cc <- control$cc # convergence criterion-tolerance
cyc <- control$cyc # max. no. of cycles
trace <- control$glm.trace # whether to print
bf.cyc <- control$bf.cyc
bf.tol <- control$bf.tol
bf.trace <- control$bf.trace

# JB: Init from here on
# dr <- 1 / pmax(exp(eta), .Machine$double.eps)
itn <- 0
lp <- eta <- log(fv) # f$linkfun(fv)
di <- -2 * gamlss.dist::dNO(y, mu, fv, log = TRUE) # f$G.di(fv) # deviance increment
dv <- sum(w * di) # the global deviance
olddv <- dv + 1 # the old global deviance
# Moved that calculation into the assignment of wv directly
# dldp <- ((y - mu)^2 - fv^2) / (fv^3) # f$dldp(fv) # u score
# d2ldp2 <- -(2 / (fv^2)) # f$d2ldp2(fv) # second derivative of log-Likelihood
# d2ldp2 <- pmin(-(2 / (fv^2)), -1e-15) # added 26-10-07
# This corresponds to: https://github.com/BerriJ/gamlss/blob/main/R/gamlss-5.R#L197
wt <- -(pmin(-(2 / (fv^2)), -1e-15) / (pmax(exp(eta), .Machine$double.eps)^-2)) #  -(d2l/dp2)/(1/(dmu/deta))^2=- (d2l/dp2)(dmu/eta)^2
# we need to stop the weights to go to Infty

wt <- pmax(pmin(wt, 1e+10), 1e-10)

# wv <- (eta-os)+step*dldp/(dr*wt)
wv <- (eta - os) + (((y - mu)^2 - fv^2) / (fv^3)) / (wt / pmax(exp(eta), .Machine$double.eps)) # eta
if (family$type == "Mixed") wv <- ifelse(is.nan(wv), 0, wv) ## TEST
iterw <- FALSE
who <- f$who
smooth.frame <- f$smooth.frame
s <- f$smooth
## starting the recycling
while (abs(olddv - dv) > cc && itn < cyc) # MS Wednesday, June 26, 2002
{
    itn <- itn + 1 # the glim inner iteration number
    lpold <- lp
    # sold <- s
    if (length(who) > 0) {
        # fit <- additive.fit(
        #     x = X, y = wv, w = wt * w, s = s, who = who, smooth.frame, maxit = bf.cyc,
        #     tol = bf.tol, trace = bf.trace
        # )
        # # lp <- fit$fitted.values
        # lp <- if (itn == 1) fit$fitted.values else step * fit$fitted.values + (1 - step) * lpold
        # #  s <- fit$smooth # test Wednesday, January 8, 2003 at 14:37
        # s <- if (itn == 1) fit$smooth else step * fit$smooth + (1 - step) * sold
    } else {
        # JB: This is where the main fitting is done
        # print(as.numeric(X))
        print(as.numeric(wv))
        fit <- lm.wfit(X, wv, rep(2, length(wv)), method = "qr")
        # JB: orig: fit <- lm.wfit(X, wv, wt * w, method = "qr")

        lp <- if (itn == 1) {
            fit$fitted.values
        } else {
            # Beware: the default step value is 1!
            step * fit$fitted.values + (1 - step) * lpold
        }
    }
    ## method 1
    # JB: That offset (os) is zero so we can ignore that
    eta <- lp + os # fixed Wednesday, September 4, 2002 at 09:45 DS
    ## own link
    fv <- pmax(exp(eta), .Machine$double.eps) # JB: orig: f$linkinv(eta)
    ## own dist
    di <- -2 * gamlss.dist::dNO(y, mu, fv, log = TRUE) # JB: orig: f$G.di(fv)
    olddv <- dv
    dv <- sum(w * di)
    ## new for automatic steps MS BR Friday, April 15, 2005 at 18:50
    # JB: This will not be executed at the moment
    # if (dv > olddv && itn >= 2 && auto == TRUE) {
    #     for (i in 1:5) # MS Thursday, September 22, 2005
    #     {
    #         lp <- (lp + lpold) / 2
    #         eta <- lp + os
    #         fv <- f$linkinv(eta)
    #         di <- f$G.di(fv)
    #         dv <- sum(w * di)
    #         #  cat("try",i,"\n")
    #         if (length(who) > 0) s <- (s + sold) / 2
    #         if ((olddv - dv) > cc) break # MS Thursday, September 22, 2005
    #     }
    # }
    if ((dv > olddv + gd.tol) && itn >= 2 && iterw == FALSE) {
        warning(
            "The deviance has increased in an inner iteration for ",
            names(formals(f$valid)), "\n", "Increase gd.tol and if persist, try different steps", "\n", "or model maybe inappropriate"
        ) #
        iterw <- TRUE
    }
    # dr <- 1 / pmax(exp(eta), .Machine$double.eps)
    ## method 2
    # dldp <- ((y - mu)^2 - fv^2) / (fv^3) # JB: orig f$dldp(fv) #
    # d2ldp2 <- -(2 / (fv^2)) # JB: orig: f$d2ldp2(fv)
    # d2ldp2 <- pmin(-(2 / (fv^2)), -1e-15) # added 26-10-07


    wt <- -(pmin(-(2 / (fv^2)), -1e-15) / (pmax(exp(eta), .Machine$double.eps)^-2))
    wt <- pmax(pmin(wt, 1e+10), 1e-10)
    wv <- (eta - os) + (((y - mu)^2 - fv^2) / (fv^3)) / (wt / pmax(exp(eta), .Machine$double.eps))


    # if (family$type == "Mixed") wv <- ifelse(is.nan(wv), 0, wv) ## TEST
    #   olddv <- dv
    #      dv <- sum(w*di)
    print(itn)
} # end of while
pen <- 0 # DS Thursday, November 21, 2002 at 23:04
if (length(who) > 0) {
    pen <- sum(eta * wt * (wv - eta))
}

# ms Saturday, December 4, 2004
# c(fit, list(fv = fv, wv = wv, wt = wt, eta = eta, os = os, pen = pen))

# Thats the final variance estimate that we are getting when calling gamlss
print(fv)
