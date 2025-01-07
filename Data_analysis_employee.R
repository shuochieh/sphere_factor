# Analysis by sectors: Product of 7 49-dimensional spheres

by_cat = "sectors"
source("./Prelim_employment.R")
source("./sphere_util.R")

x = sqrt(x)
n = dim(x)[1]
d = dim(x)[2]
p = dim(x)[3]

# Estimate mu
mu_hat = array(NA, dim = dim(x[1,,]))
log_x = array(NA, dim = c(n, p * d))
for (i in 1:d) {
  mu_hat[i,] = mean_on_sphere(x[,i,])
  log_x[,((i - 1) * p + 1):(i * p)] = Log_sphere(x[,i,], mu_hat[i,])
}

r = 20
model = LYB_fm(log_x, r = r, h = 5)
V = model$V
Fac = model$f_hat

# Fraction of variance unexplained:
FVU = rep(0, r)
for (k in 1:r) {
  z_hat = Fac[,1:k] %*% t(V[,1:k])
  x_hat = array(NA, dim = c(n, d, p))
  TV = 0
  for (i in 1:d) {
    x_hat[,i,] = Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,])
    FVU[k] = FVU[k] + sum(geod_sphere(x[,i,], x_hat[,i,])^2)
    TV = TV + sum(geod_sphere(x[,i,], mu_hat[i,])^2)
  }
  FVU[k] = FVU[k] / TV
}

# Scree plot
plot(FVU, type = "b", 
     xlab = "Num. of factors", ylab = "Fraction of variance unexplained",
     main = paste0("Scree plot (by ", by_cat, ")"))