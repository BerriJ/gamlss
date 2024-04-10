load("sigma_glim_fit.RData") # Loads y, mu, sigma.X

os <- 0 # Offset
fv <- rep(sd(y), length(y))

eta <- log(fv) # Link function

wt <- -(pmin(-(2 / (fv^2)), -1e-15) /
    (pmax(exp(eta), .Machine$double.eps)^-2))
wt <- pmax(pmin(wt, 1e+10), 1e-10)
wv <- (eta - os) + (((y - mu)^2 - fv^2) / (fv^3)) /
    (wt / pmax(exp(eta), .Machine$double.eps))

di <- -2 * gamlss.dist::dNO(y, mu, fv, log = TRUE)
dv <- sum(w * di)

# Init olddv so that the while loop runs at least once
olddv <- dv + 1
itr <- 0

while (abs(olddv - dv) > 0.001) # MS Wednesday, June 26, 2002
{
    itr <- itr + 1
    fit <- lm.wfit(sigma.X, wv, rep(2, length(wv)), method = "qr")

    eta <- fit$fitted.values
    fv <- pmax(exp(eta), .Machine$double.eps) # Inverse link function
    di <- -2 * gamlss.dist::dNO(y, mu, fv, log = TRUE) #
    olddv <- dv
    dv <- sum(w * di)

    wt <- -(pmin(-(2 / (fv^2)), -1e-15) /
        (pmax(exp(eta), .Machine$double.eps)^-2))
    wt <- pmax(pmin(wt, 1e+10), 1e-10)
    wv <- (eta - os) + (((y - mu)^2 - fv^2) / (fv^3)) /
        (wt / pmax(exp(eta), .Machine$double.eps))

    print(itr)
}

fv
