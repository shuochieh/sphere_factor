by_cat = "sectors"

# by_cat = sectors: Analysis by sectors, viewing data as product of 7 49-dimensional spheres
#                   (d = 7, p = 50)
# by_cat = states: Analysis by states, viewing data as product of 50 6-dimensional spheres
#                   (d = 50, p = 7)


source("./Prelim_employment.R")
source("./sphere_util.R")

d = dim(x)[2]
p = dim(x)[3]
n = 360 # max: 418

if (n < 418) {
  x_test = x[-c(1:n),,]
  x_test = sqrt(x_test)
  log_x_test = array(NA, dim = c(418 - n, p * d))
}

x = x[1:n,,]
x = sqrt(x)


# Estimate mu
mu_hat = array(NA, dim = dim(x[1,,]))
log_x = array(NA, dim = c(n, p * d))
for (i in 1:d) {
  mu_hat[i,] = mean_on_sphere(x[,i,])
  log_x[,((i - 1) * p + 1):(i * p)] = Log_sphere(x[,i,], mu_hat[i,])
  log_x_test[,((i - 1) * p + 1):(i * p)] = Log_sphere(x_test[,i,], mu_hat[i,])
}

r = 15
model = LYB_fm(log_x, r = r, h = 5)
V = model$V
Fac = model$f_hat

# Fraction of variance unexplained:
FVU = rep(0, r)
for (k in 1:r) {
  z_hat = Fac[,1:k] %*% t(V[,1:k])
  z_hat = t(t(z_hat) + model$mean)
  x_hat = array(NA, dim = c(n, d, p))
  TV = 0
  for (i in 1:d) {
    x_hat[,i,] = Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,])
    FVU[k] = FVU[k] + sum(geod_sphere(x[,i,], x_hat[,i,])^2)
    TV = TV + sum(geod_sphere(x[,i,], mu_hat[i,])^2)
  }
  FVU[k] = FVU[k] / TV
}
FVE = 1 - FVU

# Prediction analysis
pred_err = rep(0, r)
for (k in 1:r) {
  z_hat = log_x_test %*% V[,1:k] %*% t(V[,1:k])
  z_hat = t(t(z_hat) + model$mean)
  x_hat = array(NA, dim = c(dim(z_hat)[1], d, p))
  TV = 0
  for (i in 1:d) {
    x_hat[,i,] = Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,])
    pred_err[k] = pred_err[k] + sum(geod_sphere(x_test[,i,], x_hat[,i,])^2)
  }
}
pred_err = sqrt(pred_err / dim(z_hat)[1])

# Compare joint with separate modeling and linear factor models
r_min_p = min(r, p - 1)
FVU_s = rep(0, r_min_p)
TV = 0
for (i in 1:d) {
  model_s = LYB_fm(log_x[,((i - 1) * p + 1):(i * p)], r = r_min_p, h = 5)
  V_s = model_s$V
  Fac_s = model_s$f_hat
  TV = TV + sum(geod_sphere(x[,i,], mu_hat[i,])^2)
  for (k in 1:r_min_p) {
    z_hat = Fac_s[,1:k] %*% t(V_s[,1:k])
    z_hat = t(t(z_hat) + model_s$mean)
    x_hat = Exp_sphere(z_hat, mu_hat[i,])
    FVU_s[k] = FVU_s[k] + sum(geod_sphere(x_hat, x[,i,])^2)
  }
}
FVU_s = FVU_s / TV
FVE_s = 1 - FVU_s

pred_err_s = rep(0, r_min_p)
for (i in 1:d) {
  model_s = LYB_fm(log_x[,((i - 1) * p + 1):(i * p)], r = r_min_p, h = 5)
  V_s = model_s$V
  for (k in 1:r_min_p) {
    z_hat = log_x_test[,((i - 1) * p + 1):(i * p)] %*% V_s[,1:k] %*% t(V_s[,1:k])
    z_hat = t(t(z_hat) + model_s$mean)
    x_hat = Exp_sphere(z_hat, mu_hat[i,])
    pred_err_s[k] = pred_err_s[k] + sum(geod_sphere(x_hat, x_test[,i,])^2)
  }
}
pred_err_s = sqrt(pred_err_s / dim(z_hat)[1])



FVU_l = rep(0, r)
temp = matrix(NA, nrow = n, ncol = d * p)
for (t in 1:n) {
  temp[t,] = c(x[t,,])
}
model_l = LYB_fm(temp, r = r, h = 5)
V_l = model_l$V
Fac_l = model_l$f_hat
for (k in 1:r) {
  TV_l = sum(t(t(temp) - model_l$mean)^2)
  z_hat = Fac_l[,1:k] %*% t(V_l[,1:k])
  x_hat = t(t(z_hat) + model_l$mean)
  FVU_l[k] = sum((x_hat - temp)^2) / TV_l
}
FVE_l = 1 - FVU_l

cat("FVE (Joint):", round(FVE, 3), "\n")
cat("FVE (Separate):", round(FVE_s, 3), "\n")
cat("FVE (Linear*):", round(FVE_l, 3), "\n")









# Scree plot
par(mfrow = c(1,1))
plot(FVU, type = "b", 
     xlab = "Num. of factors", ylab = "Fraction of variance unexplained",
     main = paste0("Scree plot (by ", by_cat, ")"))

# Plot factors
par(mfrow = c(2, 2))
for (i in 1:12) {
  plot(x = seq(from = as.Date("1990-01-01"), to = as.Date("2024-10-01"), by = "month")[1:n],
       y = Fac[,i],
       type = "l",
       xlab = "Time",
       ylab = "",
       main = paste0("Factor ", i))
}

# Scatterplot of states by leading two loadings
coordinates = array(NA, dim = c(50, 7, 2))
if (by_cat == "sectors") {
  for (i in 1:50) {
    state_idx = i + 50 * c(0:6)
    coordinates[i,,] = V[state_idx, 2:3]
  }
} else if (by_cat == "states") {
  for (i in 1:50) {
    state_idx = ((i - 1) * 7 + 1):(i * 7)
    coordinates[i,,] = V[state_idx, 2:3]
  }
}

par(mfrow = c(2,2))
for (i in 1:7) {
  plot(x = coordinates[,i,1], y = coordinates[,i,2], xlab = "", ylab = "",
       main = sectors[i])
}
















