by_cat = "states"
# by_cat = sectors: Analysis by sectors, viewing data as product of 7 49-dimensional spheres
#                   (d = 7, p = 50)
# by_cat = states: Analysis by states, viewing data as product of 50 6-dimensional spheres
#                   (d = 50, p = 7)
reduced_states = FALSE 

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
n = 346 # max: 418; if smooth_x is used then max = 406

if (n < ifelse(smooth, 406, 418)) {
  x_test = x[-c(1:n),,]
  log_x_test = array(NA, dim = c(ifelse(smooth, 406, 418) - n, p * d))
}
x = x[1:n,,]
x = sqrt(x)
x_test = sqrt(x_test)


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
FVU_Euclidean = rep(0, r)
for (k in 1:r) {
  TV = 0
  TV_Euclidean = 0
  
  z_hat = predict_fm(V[,1:k], model$mean, log_x)
  for (i in 1:d) {
    FVU[k] = FVU[k] + sum(geod_sphere(x[,i,], 
                                      Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,]))^2)
    TV = TV + sum(geod_sphere(x[,i,], mu_hat[i,])^2)
    
    FVU_Euclidean[k] = FVU_Euclidean[k] + 
      sum(abs(x[,i,] - Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,]))) 
    TV_Euclidean = TV_Euclidean + sum(abs(t(x[,i,]) - mu_hat[i,]))
  }
  FVU[k] = FVU[k] / TV
  FVU_Euclidean[k] = FVU_Euclidean[k] / TV_Euclidean
}
FVE = 1 - FVU
FVE_Euclidean = 1 - FVU_Euclidean

pred_err = rep(0, r)
pred_err_Euclidean = rep(0, r)
for (k in 1:r) {
  z_hat = predict_fm(V[,1:k], model$mean, log_x_test)
  for (i in 1:d) {
    x_hat = Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,])
    pred_err[k] = pred_err[k] + sum(geod_sphere(x_test[,i,], x_hat)^2)
    pred_err_Euclidean[k] = pred_err_Euclidean[k] + sum(abs(x_test[,i,] - x_hat))
  }
}
pred_err = sqrt(pred_err / dim(x_test)[1])
pred_err_Euclidean = pred_err_Euclidean / dim(x_test)[1] # sqrt(pred_err_Euclidean / dim(x_test)[1])

# Compare: linear factor model applied to square-root transform data
FVU_linear = rep(0, r)
FVU_linear_Euclidean = rep(0, r)

temp_dta = matrix(NA, nrow = n, ncol = d * p)
for (t in 1:n) {
  for (j in 1:d) {
    temp_dta[t,((j - 1) * p + 1):(j * p)] = x[t,j,]
  }
}

model_linear = LYB_fm(temp_dta, r = r, h = 5)
V_linear = model_linear$V
Fac_linear = model_linear$f_hat
for (k in 1:r) {
  x_hat = predict_fm(V_linear[,1:k], model_linear$mean, temp_dta)
  FVU_linear_Euclidean[k] = sum(abs(x_hat - temp_dta)) / TV_Euclidean #^2) / TV_Euclidean
  
  x_hat_proj_sphere = array(NA, dim = c(n, d, p))
  for (j in 1:d) {
    x_hat_proj_sphere[,j,] = x_hat[,((j - 1) * p + 1):(j * p)] / apply(x_hat[,((j - 1) * p + 1):(j * p)], 1, norm, "2")
    FVU_linear[k] = FVU_linear[k] + sum(geod_sphere(x[,j,],
                                                    x_hat_proj_sphere[,j,])^2)
  }
  FVU_linear[k] = FVU_linear[k] / TV
}
FVE_linear = 1 - FVU_linear
FVE_linear_Euclidean = 1 - FVU_linear_Euclidean

pred_err_linear = rep(0, r)
pred_err_linear_Euclidean = rep(0, r)

temp_test = matrix(NA, nrow = dim(x_test)[1], ncol = d * p)
for (t in 1:dim(x_test)[1]) {
  for (j in 1:d) {
    temp_test[t,((j - 1) * p + 1):(j * p)] = x_test[t,j,]
  }
}

for (k in 1:r) {
  x_hat = predict_fm(V_linear[,1:k], model_linear$mean, temp_test)
  pred_err_linear_Euclidean[k] = pred_err_linear_Euclidean[k] + sum(abs(x_hat - temp_test))

  for (j in 1:d) {
    temp = x_hat[,((j - 1) * p + 1):(j * p)] / apply(x_hat[,((j - 1) * p + 1):(j * p)], 1, norm, "2")
    pred_err_linear[k] = pred_err_linear[k] + sum(geod_sphere(x_test[,j,],
                                                              temp)^2)
  }
}
pred_err_linear_Euclidean = pred_err_linear_Euclidean / dim(x_test)[1] # sqrt(pred_err_linear_Euclidean / dim(x_test)[1])
pred_err_linear = sqrt(pred_err_linear / dim(x_test)[1])


cat("FVE (Joint):", round(FVE, 3), "\n")
cat("FVE (Linear*):", round(FVE_linear, 3), "\n")

cat("FVE (Joint):", round(FVE_Euclidean, 3), "\n")
cat("FVE (Linear*):", round(FVE_linear_Euclidean, 3), "\n")

cat("Prediction error (Joint*):", round(pred_err_Euclidean, 3), "\n")
cat("Prediction error (Linear*):", round(pred_err_linear_Euclidean, 3), "\n")

cat("Prediction error (Joint*):", round(pred_err, 3), "\n")
cat("Prediction error (Linear*):", round(pred_err_linear, 3), "\n")

# Scree plot
par(mfrow = c(1,1))
plot(FVU_Euclidean, type = "b", 
     xlab = "Num. of factors", ylab = "Fraction of variance unexplained (l1 norm)",
     main = paste0("Scree plot"),
     lwd = 1.5, ylim = c(0, 0.5))
lines(FVU_linear_Euclidean, col = 4, lty = 2, lwd = 1.5, type = "b")
legend("topright", 
       legend = c("Multi-Sphere factor model", "Linear factor model"),
       col = c(1, 4), lty = c(1, 2), lwd = 1.5, 
       pch = c(1, 1), bty = "n")

par(mfrow = c(1,1))
plot(pred_err_Euclidean, type = "b", 
     xlab = "Num. of factors", 
     ylab = "l1-distance",
     main = paste0("Prediction errors"),
     lwd = 1.5) 
lines(pred_err_linear_Euclidean, col = 4, lty = 2, lwd = 1.5, type = "b")
legend("topright", 
       legend = c("Multi-Sphere factor model", "Linear factor model"),
       col = c(1, 4), lty = c(1, 2), lwd = 1.5, 
       pch = c(1, 1), bty = "n")

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
coordinates = array(NA, dim = c(45, 3, 2))
order_eigen = c(1, 2)

if (by_cat == "sectors") {
  for (i in 1:45) {
    state_idx = i + 45 * c(0:2)
    coordinates[i,,] = V[state_idx, order_eigen]
  }
} else if (by_cat == "states") {
  for (i in 1:45) {
    state_idx = ((i - 1) * 3 + 1):(i * 3)
    coordinates[i,,] = V[state_idx, order_eigen]
  }
}

par(mfrow = c(1, 3))
for (i in 1:3) {
  plot(x = coordinates[which(states %in% NE1),i,1], 
       y = coordinates[which(states %in% NE1),i,2], 
       main = sectors[i], col = 1, pch = 1,
       xlim = range(coordinates[,i,1]),
       ylim = range(coordinates[,i,2]),
       xlab = paste0(order_eigen[1], "st loading"),
       ylab = paste0(order_eigen[2], "nd loading"),
       cex = 2)
  points(x = coordinates[which(states %in% South1),i,1],
         y = coordinates[which(states %in% South1),i,2],
         col = "gray80", pch = 8, cex = 2)
  points(x = coordinates[which(states %in% South2),i,1],
         y = coordinates[which(states %in% South2),i,2],
         col = "gray80", pch = 7, cex = 2)
  points(x = coordinates[which(states %in% South3),i,1],
         y = coordinates[which(states %in% South3),i,2],
         col = "gray80", pch = 9, cex = 2)
  points(x = coordinates[which(states %in% NE2),i,1],
         y = coordinates[which(states %in% NE2),i,2],
         col = 1, pch = 16, cex = 2)
  points(x = coordinates[which(states %in% MW1),i,1],
         y = coordinates[which(states %in% MW1),i,2],
         col = 2, pch = 2, cex = 2)
  points(x = coordinates[which(states %in% MW2),i,1],
         y = coordinates[which(states %in% MW2),i,2],
         col = 2, pch = 17, cex = 2)
  points(x = coordinates[which(states %in% West1),i,1],
         y = coordinates[which(states %in% West1),i,2],
         col = 4, pch = 18, cex = 2)
  points(x = coordinates[which(states %in% West2),i,1],
         y = coordinates[which(states %in% West2),i,2],
         col = 4, pch = 5, cex = 2)
  if (i == 1) {
    legend("bottomright", 
           legend = c("New England", "Middle Atlantic",
                      "East North Central", "West North Central",
                      "South Atlantic", "East South Central", "West South Central",
                      "Mountain", "Pacific"),
           pch = c(1, 16,
                   2, 17, 
                   8, 7, 9,
                   18, 5),
           col = c(1, 1,
                   2, 2,
                   "gray80", "gray80", "gray80",
                   4, 4),
           cex = 1.5)
  }
}

# Time plots
par(mfrow = c(3, 3))
temp = which(states %in% c(NE1, NE2))

plot(y = x[,temp[1],1]^2, type = "l", main = states[temp[1]], 
     ylim = c(0, 0.9),
     x = seq(from = as.Date("1991-01-01"), by = "month", length.out = 406),
     xlab = "",
     ylab = "")
lines(x = seq(from = as.Date("1991-01-01"), by = "month", length.out = 406),
      y = x[,temp[1],2]^2, col = 2, lty = 4, lwd = 1.5)
lines(x = seq(from = as.Date("1991-01-01"), by = "month", length.out = 406),
      y = x[,temp[1],3]^2, col = 3, lty = 2, lwd = 1.5)
legend("top",
       legend = c("NRMN", "CONS", "MFG"),
       col = c(1, 2, 3),
       lty = c(1, 4, 2), 
       lwd = c(1, 1.5, 1.5),
       horiz = T)

for (j in 2:length(temp)) {
  plot(y = x[,temp[j],1]^2, type = "l", main = states[temp[j]], 
       ylim = range(x[,temp[j],]^2),
       x = seq(from = as.Date("1991-01-01"), by = "month", length.out = 406),
       xlab = "",
       ylab = "")
  lines(x = seq(from = as.Date("1991-01-01"), by = "month", length.out = 406),
        y = x[,temp[j],2]^2, col = 2, lty = 4, lwd = 1.5)
  lines(x = seq(from = as.Date("1991-01-01"), by = "month", length.out = 406),
        y = x[,temp[j],3]^2, col = 3, lty = 2, lwd = 1.5)
}
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


















