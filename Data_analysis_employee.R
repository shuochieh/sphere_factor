by_cat = "states"
# by_cat = sectors: Analysis by sectors, viewing data as product of 7 49-dimensional spheres
#                   (d = 7, p = 50)
# by_cat = states: Analysis by states, viewing data as product of 50 6-dimensional spheres
#                   (d = 50, p = 7)

source("./Prelim_employment.R")
source("./sphere_util.R")

smooth = TRUE # if moving average of x is used
if (smooth) {
  x = smooth_x
  dates = seq(from = as.Date("1991-01-01"), to = as.Date("2024-10-01"), by = "month")
} else {
  dates = seq(from = as.Date("1990-01-01"), to = as.Date("2024-10-01"), by = "month")
}

d = dim(x)[2]
p = dim(x)[3]
n = 348 # max: 418; if smooth_x is used then max = 406

if (n < ifelse(smooth, 406, 418)) {
  x_test = x[-c(1:n),,]
  x_test = sqrt(x_test)
  log_x_test = array(NA, dim = c(ifelse(smooth, 406, 418) - n, p * d))
}

x = x[1:n,,]
x = sqrt(x)


# Estimate mu
mu_hat = array(NA, dim = dim(x[1,,]))
log_x = array(NA, dim = c(n, p * d))
dis = array(NA, dim = c(n, d))
for (i in 1:d) {
  mu_hat[i,] = mean_on_sphere(x[,i,])
  log_x[,((i - 1) * p + 1):(i * p)] = Log_sphere(x[,i,], mu_hat[i,])
  log_x_test[,((i - 1) * p + 1):(i * p)] = Log_sphere(x_test[,i,], mu_hat[i,])
  
  dis[,i] = geod_sphere(mu_hat[i,], x[,i,])
}

par(mfrow = c(1, 1))
plot(x = dates[1:n],
     y = dis[,1], type = "l", lwd = 0.5, xlab = "", ylab = "", ylim = range(dis),
     main = paste0("Geodesic distance to the Frechet means (by ", by_cat, ")"))
for (i in 2:(d - 1)) {
  lines(x = dates[1:n],
        y = dis[,i], lwd = 0.5, col = i)
}

r = 15
model = LYB_fm(log_x, r = r, h = 5)
V = model$V
Fac = model$f_hat

# Fraction of variance unexplained:
FVU = rep(0, r)
for (k in 1:r) {
  TV = 0
  z_hat = predict_fm(V[,1:k], model$mean, log_x)
  for (i in 1:d) {
    FVU[k] = FVU[k] + sum(geod_sphere(x[,i,], 
                                      Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,]))^2)
    TV = TV + sum(geod_sphere(x[,i,], mu_hat[i,])^2)
  }
  FVU[k] = FVU[k] / TV
}
FVE = 1 - FVU

# Prediction analysis
pred_err = rep(0, r)
pred_err_E = rep(0, r)
for (k in 1:r) {
  z_hat = predict_fm(V[,1:k], model$mean, log_x_test)
  for (i in 1:d) {
    x_hat = Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,])
    pred_err[k] = pred_err[k] + sum(geod_sphere(x_test[,i,], x_hat)^2)
    pred_err_E[k] = pred_err_E[k] + sum((x_test[,i,] - x_hat)^2)
  }
}
pred_err = sqrt(pred_err / dim(z_hat)[1])
pred_err_E = sqrt(pred_err_E / dim(z_hat)[1])


# Compare joint with separate modeling and linear factor models
r_min_p = min(r, p - 1)
FVU_s = rep(0, r_min_p)
TV = 0
for (i in 1:d) {
  model_s = LYB_fm(log_x[,((i - 1) * p + 1):(i * p)], r = r_min_p, h = 5)
  TV = TV + sum(geod_sphere(x[,i,], mu_hat[i,])^2)
  for (k in 1:r_min_p) {
    z_hat = predict_fm(model_s$V[,1:k], model_s$mean, log_x[,((i - 1) * p + 1):(i * p)])
    x_hat = Exp_sphere(z_hat, mu_hat[i,])
    FVU_s[k] = FVU_s[k] + sum(geod_sphere(x_hat, x[,i,])^2)
  }
}
FVU_s = FVU_s / TV
FVE_s = 1 - FVU_s

pred_err_s = rep(0, r_min_p)
pred_err_sE = rep(0, r_min_p)
for (i in 1:d) {
  model_s = LYB_fm(log_x[,((i - 1) * p + 1):(i * p)], r = r_min_p, h = 5)
  V_s = model_s$V
  for (k in 1:r_min_p) {
    z_hat = predict_fm(V_s[,1:k], model_s$mean, log_x_test[,((i - 1) * p + 1):(i * p)])
    x_hat = Exp_sphere(z_hat, mu_hat[i,])
    pred_err_s[k] = pred_err_s[k] + sum(geod_sphere(x_hat, x_test[,i,])^2)
    pred_err_sE[k] = pred_err_sE[k] + sum((x_hat - x_test[,i,])^2)
  }
}
pred_err_s = sqrt(pred_err_s / dim(z_hat)[1])
pred_err_sE = sqrt(pred_err_sE / dim(z_hat)[1])


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
  x_hat = predict_fm(V_l[,1:k], model_l$mean, temp)
  FVU_l[k] = sum((x_hat - temp)^2) / TV_l
}
FVE_l = 1 - FVU_l

pred_err_l = rep(0, r)
temp_test = matrix(NA, nrow = dim(x_test)[1], ncol = d * p)
for (t in 1:dim(x_test)[1]) {
  temp_test[t,] = c(x_test[t,,])
}
for (k in 1:r) {
  x_hat = predict_fm(V_l[,1:k], model_l$mean, temp_test)
  pred_err_l[k] = pred_err_l[k] + sum((x_hat - temp_test)^2)
}
pred_err_l = sqrt(pred_err_l / dim(z_hat)[1])


cat("FVE (Joint):", round(FVE, 3), "\n")
cat("FVE (Separate):", round(FVE_s, 3), "\n")
cat("FVE (Linear*):", round(FVE_l, 3), "\n")

cat("Prediction error (Joint*):", round(pred_err_E, 3), "\n")
cat("Prediction error (Separate*):", round(pred_err_sE, 3), "\n")
cat("Prediction error (Linear*):", round(pred_err_l, 3), "\n")

# Scree plot
par(mfrow = c(1,1))
plot(FVU, type = "b", 
     xlab = "Num. of factors", ylab = "Fraction of variance unexplained",
     main = paste0("Scree plot (by ", by_cat, ")"),
     lwd = 1.5, ylim = c(0, 0.25))
lines(FVU_s, col = 2, lty = 2, lwd = 1.5, type = "b")
lines(FVU_l, col = 4, lty = 2, lwd = 1.5, type = "b")
legend("topright", 
       legend = c("Joint", "Separate", "Linear*"),
       col = c(1, 2, 4), lty = c(1, 2, 2), lwd = 1.5, 
       pch = c(1, 1, 1), bty = "n")

par(mfrow = c(1,1))
plot(pred_err_E, type = "b", 
     xlab = "Num. of factors", 
     ylab = "RMSE",
     main = paste0("Euclidean prediction errors (by ", by_cat, ")"),
     lwd = 1.5, ylim = c(0, 0.25))
lines(pred_err_sE, col = 2, lty = 2, lwd = 1.5, type = "b")
lines(pred_err_l, col = 4, lty = 2, lwd = 1.5, type = "b")
legend("topright", 
       legend = c("Joint", "Separate", "Linear"),
       col = c(1, 2, 4), lty = c(1, 2, 2), lwd = 1.5, 
       pch = c(1, 1, 1), bty = "n")

# Plot factors
par(mfrow = c(3, 2))
for (i in 1:12) {
  plot(x = dates[1:n],
       y = Fac[,i],
       type = "l",
       xlab = "Time",
       ylab = "",
       main = paste0("Factor ", i, " (by ", by_cat, ")"))
}

# Scatterplot of states by leading two loadings
coordinates = array(NA, dim = c(50, 7, 2))
order_eigen = c(2, 3)

if (by_cat == "sectors") {
  for (i in 1:50) {
    state_idx = i + 50 * c(0:6)
    coordinates[i,,] = V[state_idx, order_eigen]
  }
} else if (by_cat == "states") {
  for (i in 1:50) {
    state_idx = ((i - 1) * 7 + 1):(i * 7)
    coordinates[i,,] = V[state_idx, order_eigen]
  }
}

par(mfrow = c(2,2))
for (i in 1:7) {
  plot(x = coordinates[which(blue_codes == 1),i,1], 
       y = coordinates[which(blue_codes == 1),i,2], 
       main = sectors[i], col = "blue",
       xlim = range(coordinates[,i,1]),
       ylim = range(coordinates[,i,2]),
       xlab = paste0(order_eigen[1], "-th loading"),
       ylab = paste0(order_eigen[2], "-th loading"))
  points(x = coordinates[which(blue_codes == 0),i,1], 
         y = coordinates[which(blue_codes == 0),i,2],
         col = "red")
}
















