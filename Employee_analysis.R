### Analysis of employment data

x = array(NA, dim = c(418, 7, 50))

for (t in 1:418) {
  x[t,,] = t(res[t,,])
}

x = sqrt(x) # square root transform

# Estimate mu
mu_hat = matrix(NA, nrow = 7, ncol = 50)
trans_x = matrix(NA, nrow = 418, ncol = 7 * 50)
for (j in 1:7) {
  mu_hat[j,] = mean_on_sphere(x[,j,])
  trans_x[,((j - 1) * 50 + 1):(j * 50)] = Log_sphere(x[,j,], mu_hat[j,])
}

colnames(mu_hat) = states
for (j in 1:length(industries)) {
  barplot(mu_hat[j,], cex.names = 0.9, las = 3, ylim = c(0, ifelse(j == 2, 0.41, 0.35)), 
          main = paste("Estimated mean composition of", industries[j]))
}
