foo <- gamlss.dist::TF()

y <- c(1.5, 2)
theta <- matrix(1:6, ncol = 3, byrow = TRUE)

mu <- theta[, 1]
sigma <- theta[, 2]
nu <- theta[, 3]

print(foo$dldm(y, mu, sigma, nu))
print(foo$dldd(y, mu, sigma, nu))
print(foo$dldv(y, mu, sigma, nu))

print(foo$d2ldm2(sigma, nu))
print(foo$d2ldd2(sigma, nu))
print(foo$d2ldv2(y, mu, sigma, nu))

print(foo$d2ldmdd(y))
print(foo$d2ldmdv(y))
print(foo$d2ldddv(sigma, nu))
