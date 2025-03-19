library(tidyr)
library(ggplot2)
library(reshape2)
source("./main_func.R")

plot_time_series <- function(x, time = NULL, series_cols = NULL, 
                             title = "Time Series Plot", 
                             x_label = "Time", 
                             y_label = "Value",
                             l_size = 1.2, 
                             p_size = 1.4,
                             legend.pos = "top",
                             legend_title = "series",
                             legend_rows = 1,
                             y_lim = NULL) {
  df_long <- pivot_longer(x, cols = all_of(series_cols), 
                          names_to = "Series", values_to = "Value")
  
  # Create a ggplot with multiple time series
  p = ggplot(df_long, aes(x = .data[[time]], y = Value, color = Series, shape = Series)) +
    geom_line(size = l_size) +  
    geom_point(size = p_size) +  
    scale_color_brewer(palette = "Dark2") + 
    scale_shape_manual(values = seq(16, 18, length.out = length(series_cols))) +  
    labs(title = title, x = x_label, y = y_label, color = legend_title, shape = legend_title) +
    theme(
      legend.position = legend.pos, 
      axis.text.x = element_text(angle = 0, size = 12),  
      axis.text.y = element_text(size = 12),  
      axis.title.x = element_text(size = 14),  
      axis.title.y = element_text(size = 14),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    ) +
    guides(color = guide_legend(nrow = legend_rows), shape = guide_legend(nrow = legend_rows))
  if (!is.null(y_lim)) {
    
    p <- p + scale_y_continuous(limits = y_lim) +
      scale_x_continuous(breaks = seq(1, 10, by = 2))  
  }
  
  return (p)
}

plot_assist = function (x, label.y, title = "", legend.if = FALSE) {
  Labels = c("Riemannian factor model", "Linear factor model", "Oracle prediction error")
  colors = c("red", "blue", "green")
  percentiles = array(NA, dim = c(3, 10, 2))
  for (i in 1:3) {
    percentiles[i,,] = cbind(apply(x[,i,], 2, quantile, 0.05), apply(x[,i,], 2, quantile, 0.95))
  }
  means = colMeans(x, dims = 1)
  
  plot(NA, xlim = range(1:10), ylim = c(0, max(percentiles)), 
       xlab = "Number of factors", ylab = label.y, main = title)
  for (i in 1:3) {
    polygon(c(1:10, 10:1),
            c(percentiles[i,,1], rev(percentiles[i,,2])),
            col = adjustcolor(colors[i], alpha.f = 0.2), border = NA)
    lines(1:10, means[i,], col = colors[i], lwd = 2, type = "b",
          pch = i - 1)
  }
  if (legend.if) {
    legend("topright", legend = Labels, col = colors, lwd = 2, pch = c(0, 1, 2), 
           bty = "n")
  }
}

dta_gen = function (n, type, manifold.type, fac_sd = 1) {
  # type: specifications
  # manifold.type: "BWS" or "Sphere"
  
  if (type == 1) {
    # BWS space p by p SPD matrices
    p = 10
    r = 5
    inj_noise = 0.0
    fac_noise = 2.0
    alpha = 0.8 # AR coefficient for factor process
  }
  
  if (type == 2) {
    # BWS space p by p SPD matrices
    p = 10
    r = 5
    inj_noise = 0.5
    fac_noise = 1.0
    alpha = 0.8 # AR coefficient for factor process
  }
  
  
  if (manifold.type == "BWS") {
    mu_P = matrix(rnorm(500 * p), ncol = p)
    mu = 3 * (t(mu_P) %*% mu_P / 500) 

    Factors = array(0, dim = c(n + 100, r))
    for (t in 2:(n + 100)) {
      Factors[t,] = alpha * Factors[t - 1,] + rnorm(r, sd = fac_sd)
    }
    Factors = Factors[-c(1:100),]
    
    A = array(rnorm(p * (p + 1) * r / 2), dim = c(p * (p + 1) / 2, r))
    A = t(t(A) / apply(A, 2, norm, "2"))

    Z = array(NA, dim = c(n, p * (p + 1) / 2))
    Z_nless = array(NA, dim = c(n, p * (p + 1) / 2))
    for (t in 1:n) {
      Z_nless[t,] = A %*% c(Factors[t,])
      Z[t,] = A %*% c(Factors[t,]) + rnorm(p * (p + 1) / 2, 
                                           sd = fac_noise * sqrt(2 / (p * (p + 1))))
    }

    X_nless = array(NA, dim = c(n, p, p))
    X = array(NA, dim = c(n, p, p))
    for (t in 1:n) {
      V = vector_to_symmetric(Z[t,], p)
      X_nless[t,,] = Exp_BWS_core(V, mu)
      
      noise_P = rnorm(p * (p + 1) / 2, sd = inj_noise * sqrt(2 / (p * (p + 1))))
      noise = vector_to_symmetric(noise_P, p)

      X[t,,] = Exp_BWS_core(noise, X_nless[t,,])
    }
    
    return (list("X" = X, "X_nless" = X_nless, "Z" = Z, "Z_nless" = Z_nless, 
                 "A" = A, "Factors" = Factors, "mu" = mu))
  }
  
}


###########################
# Fix p = 10, r = 5, vary n = 50, 100, 200
set.seed(5566)  
num_sim = 100
type = 1

n = 50
p = 10
n_test = 200

FVUs = array(0, dim = c(num_sim, 3, 10))
PEs = array(0, dim = c(num_sim, 3, 10))

for (zz in 1:num_sim) {
  dta = dta_gen(n + n_test, type, "BWS", fac_sd = 2)
  RFM = main_BWS(dta$X, 10, verbose = F, test_size = n_test, mu_tol = 1e-3, h = 6)
  
  oracle_noise = 0
  for (i in 1:n_test) {
    oracle_noise = oracle_noise + 
      sum((dta$X_nless[n + i,,] - Exp_BWS_core(vector_to_symmetric(dta$Z_nless[i + n,], p), dta$mu)))^2
  }
  oracle_noise = sqrt(oracle_noise / n_test)
  
  pred_err_RFM = rep(0, 10)
  for (i in 1:10) {
    for (t in 1:n_test) {
      pred_err_RFM[i] = pred_err_RFM[i] + norm(RFM$x_hat_RFM[t,i,,] - dta$X_nless[-c(1:n),,][t,,], "F")^2
    }
  }
  pred_err_RFM = sqrt(pred_err_RFM / n_test)
  
  pred_err_LFM = rep(0, 10)
  for (i in 1:10) {
    for (t in 1:n_test) {
      pred_err_LFM[i] = pred_err_LFM[i] + norm(RFM$x_hat_LFM[t,i,,] - dta$X_nless[-c(1:n),,][t,,], "F")^2
    }
  }
  pred_err_LFM = sqrt(pred_err_LFM / n_test)
  
  PEs[zz,1,] = pred_err_RFM
  PEs[zz,2,] = pred_err_LFM
  PEs[zz,3,] = oracle_noise
  
  oracle_noise = 0
  for (i in 1:n) {
    oracle_noise = oracle_noise + 
      sum((dta$X[i,,] - Exp_BWS_core(vector_to_symmetric(dta$Z_nless[i,], p), dta$mu)))^2
  }
  
  FVUs[zz,1,] = RFM$FVU_e
  FVUs[zz,2,] = RFM$FVU_e_linear
  FVUs[zz,3,] = oracle_noise / RFM$TV_e
  
  cat("iteration", zz, "\n")
}

plot_assist(PEs, "Prediction errors", paste("n =", n, "; p =", p), T)
plot_assist(FVUs, "Fraction of errors unexplained", "", F)

################################

n = 100
p = 10
n_test = 200

FVUs = array(0, dim = c(num_sim, 3, 10))
PEs = array(0, dim = c(num_sim, 3, 10))

for (zz in 1:num_sim) {
  dta = dta_gen(n + n_test, type, "BWS", fac_sd = 2)
  RFM = main_BWS(dta$X, 10, verbose = F, test_size = n_test, mu_tol = 1e-3, h = 6)
  
  oracle_noise = 0
  for (i in 1:n_test) {
    oracle_noise = oracle_noise + 
      sum((dta$X_nless[n + i,,] - Exp_BWS_core(vector_to_symmetric(dta$Z_nless[i + n,], p), dta$mu)))^2
  }
  oracle_noise = sqrt(oracle_noise / n_test)
  
  pred_err_RFM = rep(0, 10)
  for (i in 1:10) {
    for (t in 1:n_test) {
      pred_err_RFM[i] = pred_err_RFM[i] + norm(RFM$x_hat_RFM[t,i,,] - dta$X_nless[-c(1:n),,][t,,], "F")^2
    }
  }
  pred_err_RFM = sqrt(pred_err_RFM / n_test)
  
  pred_err_LFM = rep(0, 10)
  for (i in 1:10) {
    for (t in 1:n_test) {
      pred_err_LFM[i] = pred_err_LFM[i] + norm(RFM$x_hat_LFM[t,i,,] - dta$X_nless[-c(1:n),,][t,,], "F")^2
    }
  }
  pred_err_LFM = sqrt(pred_err_LFM / n_test)
  
  PEs[zz,1,] = pred_err_RFM
  PEs[zz,2,] = pred_err_LFM
  PEs[zz,3,] = oracle_noise
  
  oracle_noise = 0
  for (i in 1:n) {
    oracle_noise = oracle_noise + 
      sum((dta$X[i,,] - Exp_BWS_core(vector_to_symmetric(dta$Z_nless[i,], p), dta$mu)))^2
  }
  
  FVUs[zz,1,] = RFM$FVU_e
  FVUs[zz,2,] = RFM$FVU_e_linear
  FVUs[zz,3,] = oracle_noise / RFM$TV_e
  
  cat("iteration", zz, "\n")
}

plot_assist(PEs, "Prediction errors", paste("n =", n, "; p =", p), T)
plot_assist(FVUs, "Fraction of errors unexplained", "", F)


################################

n = 200
p = 10
n_test = 200

FVUs = array(0, dim = c(num_sim, 3, 10))
PEs = array(0, dim = c(num_sim, 3, 10))

for (zz in 1:num_sim) {
  dta = dta_gen(n + n_test, type, "BWS", fac_sd = 2)
  RFM = main_BWS(dta$X, 10, verbose = F, test_size = n_test, mu_tol = 1e-3, h = 6)
  
  oracle_noise = 0
  for (i in 1:n_test) {
    oracle_noise = oracle_noise + 
      sum((dta$X_nless[n + i,,] - Exp_BWS_core(vector_to_symmetric(dta$Z_nless[i + n,], p), dta$mu)))^2
  }
  oracle_noise = sqrt(oracle_noise / n_test)
  
  pred_err_RFM = rep(0, 10)
  for (i in 1:10) {
    for (t in 1:n_test) {
      pred_err_RFM[i] = pred_err_RFM[i] + norm(RFM$x_hat_RFM[t,i,,] - dta$X_nless[-c(1:n),,][t,,], "F")^2
    }
  }
  pred_err_RFM = sqrt(pred_err_RFM / n_test)
  
  pred_err_LFM = rep(0, 10)
  for (i in 1:10) {
    for (t in 1:n_test) {
      pred_err_LFM[i] = pred_err_LFM[i] + norm(RFM$x_hat_LFM[t,i,,] - dta$X_nless[-c(1:n),,][t,,], "F")^2
    }
  }
  pred_err_LFM = sqrt(pred_err_LFM / n_test)
  
  PEs[zz,1,] = pred_err_RFM
  PEs[zz,2,] = pred_err_LFM
  PEs[zz,3,] = oracle_noise
  
  oracle_noise = 0
  for (i in 1:n) {
    oracle_noise = oracle_noise + 
      sum((dta$X[i,,] - Exp_BWS_core(vector_to_symmetric(dta$Z_nless[i,], p), dta$mu)))^2
  }
  
  FVUs[zz,1,] = RFM$FVU_e
  FVUs[zz,2,] = RFM$FVU_e_linear
  FVUs[zz,3,] = oracle_noise / RFM$TV_e
  
  cat("iteration", zz, "\n")
}

plot_assist(PEs, "Prediction errors", paste("n =", n, "; p =", p), T)
plot_assist(FVUs, "Fraction of errors unexplained", "", F)

# cmPEs = colMeans(PEs, dims = 1)
# tsp = plot_time_series(x = data.frame(Time = c(1:10),
#                                       "Riemannian Factor Model" = cmPEs[1,],
#                                       "Linear Factor Model" = cmPEs[2,],
#                                       "Oracle prediction error" = cmPEs[3,],
#                                       check.names = F),
#                        time = "Time",
#                        series_cols = c("Riemannian Factor Model", "Linear Factor Model", "Oracle prediction error"),
#                        title = "",
#                        y_label = "(Euclidean) Prediction errors",
#                        x_label = "Number of factors",
#                        l_size = 0.5,
#                        p_size = 3.5,
#                        legend_title = "",
#                        legend_rows = 2,
#                        y_lim = c(0, 16))
# print(tsp)
# 
# 
# cmFVUs = colMeans(FVUs, dims = 1)
# tsp = plot_time_series(x = data.frame(Time = c(1:10),
#                                       "Riemannian Factor Model" = cmFVUs[1,],
#                                       "Linear Factor Model" = cmFVUs[2,],
#                                       "Oracle FVU" = cmFVUs[3,],
#                                       check.names = F),
#                        time = "Time",
#                        series_cols = c("Riemannian Factor Model", "Linear Factor Model", "Oracle FVU"),
#                        title = "",
#                        y_label = "Fraction of variance unexplained",
#                        x_label = "Number of factors",
#                        l_size = 0.5,
#                        p_size = 3.5,
#                        legend_title = "",
#                        legend_rows = 2,
#                        y_lim = c(0, 1))
# print(tsp)


################################
