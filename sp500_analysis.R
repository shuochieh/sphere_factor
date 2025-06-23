library(MASS)
library(reshape2)
library(gridExtra)
library(ggplot2)
source("./main_func.R")
source("./BWS_util.R")

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
    theme_minimal(base_size = 15) +  # Minimal theme with larger base font
    theme(
      legend.position = legend.pos, 
      axis.text.x = element_text(angle = 0),  
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    ) +
    guides(color = guide_legend(nrow = legend_rows), shape = guide_legend(nrow = legend_rows))
  if (!is.null(y_lim)) {
    p <- p + scale_y_continuous(limits = y_lim)
  }
  
  return (p)
}

plot_assist = function (res1, res2 = NULL, oracle = NULL, x_vals = NULL,
                        labs = NULL, ylim = NULL, main = NULL,
                        fraction = TRUE) {
  mean1 = res1

  if (!is.null(res2)) {
    mean2 = res2
  }
  
  if (is.null(ylim)) {
    if (fraction) {
      ylim = c(0, 1)
    } else {
      ylim = range(c(c(res1), c(res2)))
    }
  }
  
  if (is.null(main)) {
    main = ""
  }
  
  if (!is.null(labs)) {
    xlab = labs[1]
    ylab = labs[2]
  } else {
    xlab = ""
    ylab = ""
  }
  
  if (is.null(x_vals)) {
    x_vals = 1:length(mean1)
  }
  plot(x = x_vals, y = mean1, type = "n",
       ylim = ylim, 
       xlab = xlab, ylab = ylab, main = main, bty = "L")
  grid(col = "lightgray", lty = "dotted", lwd = 1)
  
  # polygon(c(x_vals, rev(x_vals)), c(q05_1, rev(q95_1)),
  #         col = adjustcolor("lightblue", alpha.f = 0.6), border = NA)
  if (!is.null(res2)) {
    # polygon(c(x_vals, rev(x_vals)), c(q05_2, rev(q95_2)),
    #         col = adjustcolor("lightsalmon", alpha.f = 0.6), border = NA)
  }
  lines(x_vals, mean1, col = "steelblue", lwd = 2)
  points(x_vals, mean1, pch = 19, col = "lightblue")
  if (!is.null(res2)) {
    lines(x_vals, mean2, col = "firebrick", lwd = 2, lty = 2)
    points(x_vals, mean2, pch = 17, col = "lightsalmon")
  }
  
  if (!is.null(oracle)) {
    abline(h = oracle, col = "black", lty = 4, lwd = 1.5)
  }
} 

heat_plot = function (A, i, lims) {
  df = melt(A)
  colnames(df) = c("Row", "Column", "Value")
  
  df$Row <- factor(df$Row, levels = rev(rownames(A)))
  df$Column <- factor(df$Column, levels = colnames(A))
  
  mid_point = mean(lims)
  
  hm = ggplot(df, aes(Column, Row, fill = Value)) +
    geom_tile() +
    # scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", 
    #                      midpoint = mid_point, limits = lims) +
    scale_fill_gradient2(low = "white", high = "darkred", limits = lims) +
    theme_minimal() +
    coord_fixed() +  
    labs(x = "", y = "", fill = "", 
         title = paste("Loading matrix", i)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return (hm)
}


dta = load("./sp500_covariance/sp500_13.RData")
dta = covariances
rm("covariances")

temp = array(aperm(dta, c(2, 1, 3, 4)), 
             c(dim(dta)[1] * dim(dta)[2], dim(dta)[c(3,4)]))

# Uncomment if wish to deal with correlation matrix
# for (i in 1:dim(temp)[1]) {
#   temp[i,,] = cov2cor(temp[i,,])
# }

dta = temp
rm("temp")

dta = dta * 10000 # Convert to percentage points
overall_covariance_training = overall_covariance_training * 10000

# mu = mean_on_BWS(dta, tau = 0.5, max.iter = 100, batch_size = NULL, verbose = TRUE)
# tsp = plot_time_series(x = data.frame(Time = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month"),
#                                       "Geodesic distance to Frechet mean" = geod_BWS(dta, mu),
#                                       check.names = F),
#                        time = "Time",
#                        series_cols = "Geodesic distance to Frechet mean",
#                        title = "Geodesic distance to the Frechet mean",
#                        y_label = "",
#                        x_label = "Time",
#                        l_size = 0.7,
#                        p_size = 0,
#                        legend.pos = "none",
#                        y_lim = c(0, 20))
# print(tsp)

# plot(x = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month"),
#      y = geod_BWS(dta, mu), 
#      type = "l", ylab = "", xlab="", 
#      main = "Geodesic distance to the Frechet mean")

results = main_BWS(dta, r = 15, test_size = 36, h = 6, batch_size = 30, 
                   max.iter = 100, return_predictions = TRUE)
plot_assist(results$FVU_RFM_BWS, res2 = results$FVU_LYB_BWS,
            labs = c("number of factors", "FUGV"), ylim = c(0.2, 1))
legend("topright",
       legend = c("RFM", "LFM"),
       col = c("steelblue", "firebrick"),
       lty = c(1, 2), lwd = c(2, 2), bty = "n")
legend("topright",
       legend = c("RFM", "LFM"),
       col = c("lightblue", "lightsalmon"),
       lty = c(1, 1), lwd = 0.01,
       pch = c(19, 17), bty = "n")

plot_assist(results$FVU_RFM_Euc, res2 = results$FVU_LYB_Euc,
            labs = c("number of factors", "FUV"), ylim = c(0.2, 1.0))
legend("topright",
       legend = c("RFM", "LFM"),
       col = c("steelblue", "firebrick"),
       lty = c(1, 2), lwd = c(2, 2), bty = "n")
legend("topright",
       legend = c("RFM", "LFM"),
       col = c("lightblue", "lightsalmon"),
       lty = c(1, 1), lwd = 0.01,
       pch = c(19, 17), bty = "n")

RFM_xhat = array(NA, dim = c(15, 36, 13, 13))
LFM_xhat = array(NA, dim = c(15, 36, 13, 13))
for (k in 1:15) {
  results = main_BWS(dta, r = k, test_size = 36, h = 6, batch_size = 30, max.iter = 100,
                     return_predictions = TRUE)
  RFM_xhat[k,,,] = results$RFM_xhat
  LFM_xhat[k,,,] = results$LYB_xhat
  
  
  cat("iteration", k, "\n")
}


diag_loss = array(0, dim = c(2, 15))
offdiag_loss = array(0, dim = c(2, 15))
for (k in 1:15) {
  for (i in 1:36) {
    truth = dta[(300 - 36 + i),,]
    prediction = LFM_xhat[k,i,,] # project_to_SPD(LFM_xhat[k,i,,])
    e = truth - prediction
    
    diag_loss[1, k] = diag_loss[1, k] + sum(diag(e)^2)
    offdiag_loss[1, k] = offdiag_loss[1, k] + sum(e[upper.tri(e)]^2)
    
    prediction = RFM_xhat[k,i,,]
    e = truth - prediction

    diag_loss[2, k] = diag_loss[2, k] + sum(diag(e)^2)
    offdiag_loss[2, k] = offdiag_loss[2, k] + sum(e[upper.tri(e)]^2)
  }
}

par(mfrow = c(1, 2))
plot_assist(diag_loss[2,] / diag_loss[1,], res2 = offdiag_loss[2,] / offdiag_loss[1,],
            ylim = c(0.7, 1.3), x_vals = 1:15,
            labs = c("number of factors", "ratio"))
legend("topright",
       legend = c("Diagonal", "Off-diagonal"),
       col = c("steelblue", "firebrick"),
       lty = c(1, 2), lwd = c(2, 2), bty = "n")
legend("topright",
       legend = c("Diagonal", "Off-diagonal"),
       col = c("lightblue", "lightsalmon"),
       lty = c(1, 1), lwd = 0.01,
       pch = c(19, 17), bty = "n")

plot_assist(sqrt((diag_loss[2,] + 2 * offdiag_loss[2,]) / 2), 
            res2 = sqrt((diag_loss[1,] + 2 * offdiag_loss[1,]) / 2),
            ylim = c(40, 75), x_vals = 1:15,
            labs = c("number of factors", "Frobenius errors"))
legend("topright",
       legend = c("RFM", "LFM"),
       col = c("steelblue", "firebrick"),
       lty = c(1, 2), lwd = c(2, 2), bty = "n")
legend("topright",
       legend = c("RFM", "LFM"),
       col = c("lightblue", "lightsalmon"),
       lty = c(1, 1), lwd = 0.01,
       pch = c(19, 17), bty = "n")


par(mfrow = c(2, 3))
for (j in 1:6) {
  plot(x = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month")[1:264],
       y = results$Factors[,j],
       type = "l",
       xlab = "", ylab = "",
       col = "steelblue",
       main = paste("Factor", j), 
       bty = "L")
  year_starts = seq(from = as.Date("2000-01-01"), to = as.Date("2024-01-01"), by = "5 years")
  abline(v = year_starts, col = "lightgray", lty = "dotted", lwd = 1)
}

ps = NULL
for (i in 1:2) {
  if (i == 1) {
    Exp_V = Exp_BWS(log_to_tangent(results$V[,i], results$E), results$mu_hat)
  } else {
    Exp_V = Exp_BWS(-log_to_tangent(results$V[,i], results$E), results$mu_hat)
  }
  Exp_V = Exp_V[c(1, 2, 3, 9, 4, 13, 5, 8, 11, 6, 7, 10, 12),]
  Exp_V = Exp_V[,c(1, 2, 3, 9, 4, 13, 5, 8, 11, 6, 7, 10, 12)]
  # diag(Exp_V) = 0
  rownames(Exp_V) = selected_companies[c(1, 2, 3, 9, 4, 13, 5, 8, 11, 6, 7, 10, 12)]
  colnames(Exp_V) = selected_companies[c(1, 2, 3, 9, 4, 13, 5, 8, 11, 6, 7, 10, 12)]
  ps[[i]] = heat_plot(Exp_V, i, lims = c(0, 15))
}

grid.arrange(ps[[1]], ps[[2]], ncol = 2)

### Out-of-sample factor forecasting

RFM_res = dyn_RFM(dta, r = 1, test_size = 36, h = 6, batch_size = 30)
LFM_res = dyn_LFM(dta, r = 1, test_size = 36, h = 6)

cos_dist = array(NA, dim = c(13, 3, 36))
BWS_errors = array(NA, dim = c(3, 36))
for (m in 1:36) {
  
  truth0 = dta[(300 - 36 + m),,]
  truth = eigen(truth0)

  for (k in 1:13) {
    # Riemannian factor model
    Sigma_hat = RFM_res[m,,]
    temp = eigen(Sigma_hat)
    v = as.matrix(temp$vectors[,1:k])
    cos_dist[k,1,m] = subspace_d(v, truth$vectors[,1:k])
    if (k == 1) {
      BWS_errors[1,m] = geod_BWS_core(Sigma_hat, truth0)^2
    }

    # Linear factor model
    Sigma_hat = LFM_res[m,,]
    temp = eigen(Sigma_hat)
    v = as.matrix(temp$vectors[,1:k])
    cos_dist[k,2,m] = subspace_d(v, truth$vectors[,1:k])
    if (k == 1) {
      BWS_errors[2,m] = geod_BWS_core(project_to_SPD(Sigma_hat), truth0)^2
    }

    # Random walk model
    Sigma_hat = dta[(300- 36 + m - 1),,]
    temp = eigen(Sigma_hat)
    v = as.matrix(temp$vectors[,1:k])
    cos_dist[k,3,m] = subspace_d(v, truth$vectors[,1:k])
    if (k == 1) {
      BWS_errors[3,m] = geod_BWS_core(Sigma_hat, truth0)^2
    }

  }
}

par(mfrow = c(1, 3))
for (k in c(1:3)) {
  x_dates = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month")[265:300]
  plot(x = x_dates, 
       y = cos_dist[k,1,], type = "n",
       ylim = c(0.0, 1.0), main = paste0("leading ", k, "-dimensional eigenspace"),
       xlab = "", ylab = "Subspace distance", bty = "L")
  grid(col = "lightgray", lty = "dotted", lwd = 1)
  lines(x_dates, cos_dist[k,1,], col = "steelblue", lwd = 2)
  points(x_dates, cos_dist[k,1,], pch = 19, col = "lightblue")
  lines(x_dates, cos_dist[k,2,], col = "firebrick", lwd = 1.5, lty = 2)
  points(x_dates, cos_dist[k,2,], pch = 17, col = "lightsalmon")
  lines(x_dates, cos_dist[k,3,], col = "darkseagreen4", lwd = 1.5, lty = 4)
  points(x_dates, cos_dist[k,3,], pch = 23, col = "darkseagreen1",
         bg = "darkseagreen1")
  
  if (k == 1) {
    legend("bottomright",
           legend = c(paste0("RFM (Mean = ", round(rowMeans(cos_dist[k,,])[1], 2), ")"), 
                      paste0("LFM (Mean = ", round(rowMeans(cos_dist[k,,])[2], 2), ")"), 
                      paste0("LOCF (Mean = ", round(rowMeans(cos_dist[k,,])[3], 2), ")")),
           col = c("steelblue", "firebrick", "darkseagreen4"),
           lty = c(1, 2, 4), lwd = c(2, 1.5, 1.5), bty = "n")
    legend("bottomright",
           legend = c(paste0("RFM (Mean = ", round(rowMeans(cos_dist[k,,])[1], 2), ")"), 
                      paste0("LFM (Mean = ", round(rowMeans(cos_dist[k,,])[2], 2), ")"), 
                      paste0("LOCF (Mean = ", round(rowMeans(cos_dist[k,,])[3], 2), ")")),
           col = c("lightblue", "lightsalmon", "darkseagreen1"),
           lty = c(1, 2, 4), lwd = 0.01,
           pch = c(19, 17, 23), bty = "n",
           pt.bg = c(NULL, NULL, "darkseagreen1"))
  } else {
    legend("bottomright",
           legend = c(paste0("RFM (Mean = ", round(rowMeans(cos_dist[k,,])[1], 2), ")"), 
                      paste0("LFM (Mean = ", round(rowMeans(cos_dist[k,,])[2], 2), ")"), 
                      paste0("LOCF (Mean = ", round(rowMeans(cos_dist[k,,])[3], 2), ")")),
           col = c("steelblue", "firebrick", "darkseagreen4"),
           lty = c(1, 2, 4), lwd = c(2, 1.5, 1.5), bty = "n")
    legend("bottomright",
           legend = c(paste0("RFM (Mean = ", round(rowMeans(cos_dist[k,,])[1], 2), ")"), 
                      paste0("LFM (Mean = ", round(rowMeans(cos_dist[k,,])[2], 2), ")"), 
                      paste0("LOCF (Mean = ", round(rowMeans(cos_dist[k,,])[3], 2), ")")),
           col = c("lightblue", "lightsalmon", "darkseagreen1"),
           lty = c(1, 2, 4), lwd = 0.01,
           pch = c(19, 17, 23), bty = "n",
           pt.bg = c(NULL, NULL, "darkseagreen1"))
  }

}

par(mfrow = c(1, 1))
x_dates = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month")[265:300]
plot(x = x_dates, 
     y = sqrt(BWS_errors[1,]), type = "n",
     ylim = c(0, 5), 
     xlab = "", ylab = "BW distance", bty = "L")
grid(col = "lightgray", lty = "dotted", lwd = 1)
lines(x_dates, sqrt(BWS_errors[1,]), col = "steelblue", lwd = 2)
points(x_dates, sqrt(BWS_errors[1,]), pch = 19, col = "lightblue")
lines(x_dates, sqrt(BWS_errors[2,]), col = "firebrick", lwd = 1.5, lty = 2)
points(x_dates, sqrt(BWS_errors[2,]), pch = 17, col = "lightsalmon")
lines(x_dates, sqrt(BWS_errors[3,]), col = "darkseagreen4", lwd = 1.5, lty = 4)
points(x_dates, sqrt(BWS_errors[3,]), pch = 23, col = "darkseagreen1",
       bg = "darkseagreen1")

legend("bottomright",
       legend = c(paste0("RFM (Mean = ", round(rowMeans(sqrt(BWS_errors))[1], 2), ")"), 
                  paste0("LFM (Mean = ", round(rowMeans(sqrt(BWS_errors))[2], 2), ")"), 
                  paste0("LOCF (Mean = ", round(rowMeans(sqrt(BWS_errors))[3], 2), ")")),
       col = c("steelblue", "firebrick", "darkseagreen4"),
       lty = c(1, 2, 4), lwd = c(2, 1.5, 1.5), bty = "n")
legend("bottomright",
       legend = c(paste0("RFM (Mean = ", round(rowMeans(sqrt(BWS_errors))[1], 2), ")"), 
                  paste0("LFM (Mean = ", round(rowMeans(sqrt(BWS_errors))[2], 2), ")"), 
                  paste0("LOCF (Mean = ", round(rowMeans(sqrt(BWS_errors))[3], 2), ")")),
       col = c("lightblue", "lightsalmon", "darkseagreen1"),
       lty = c(1, 2, 4), lwd = 0.00,
       pch = c(19, 17, 23), bty = "n",
       pt.bg = c(NULL, NULL, "darkseagreen1"))

risk_error = array(NA, dim = c(3, 36))
for (m in 1:36) {
  truth = dta[(300 - 36 + m),,]
  truth_lag = dta[(300 - 36 + m - 1),,]
  w = solve(truth_lag) %*% rep(1, 13)
  w = w / sum(w)
  true_risk = c(t(w) %*% truth %*% w)
  
  # Riemannian factor model
  Sigma_hat = RFM_res[m,,]
  predicted_risk = c(t(w) %*% Sigma_hat %*% w)
  risk_error[1,m] = (predicted_risk - true_risk)^2
  
  # Linear factor model
  Sigma_hat = LFM_res[m,,]
  predicted_risk = c(t(w) %*% Sigma_hat %*% w)
  risk_error[2,m] = (predicted_risk - true_risk)^2
  
  # Random walk model
  Sigma_hat = dta[(300- 36 + m - 1),,]
  predicted_risk = c(t(w) %*% Sigma_hat %*% w)
  risk_error[3,m] = (predicted_risk - true_risk)^2
}
risk_error = sqrt(risk_error)
x_dates = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month")[265:300]
plot(x = x_dates,
     y = risk_error[1,],
     type = "n", ylim = c(0, 12),
     xlab = "", ylab = "risk prediction error", bty = "L")
grid(col = "lightgray", lty = "dotted", lwd = 1)

lines(x_dates, risk_error[1,], col = "steelblue", lwd = 2)
points(x_dates, risk_error[1,], col = "lightblue", pch = 19, cex = 0.8)
lines(x_dates, risk_error[2,], col = "firebrick", lwd = 1.5, lty = 2)
points(x_dates, risk_error[2,], col = "lightsalmon", pch = 17, cex = 0.8)
lines(x_dates, risk_error[3,], col = "darkseagreen4", lwd = 1.5, lty = 4)
points(x_dates, risk_error[3,], col = "darkseagreen1", bg = "darkseagreen1",
       pch = 23, cex = 0.8)

legend("topright",
       legend = c(paste0("RFM (Mean = ", round(mean(risk_error[1,]), 2), ")"), 
                  paste0("LFM (Mean = ", round(mean(risk_error[2,]), 2), ")"), 
                  paste0("LOCF (Mean = ", round(mean(risk_error[3,]), 2), ")")),
       col = c("steelblue", "firebrick", "darkseagreen4"),
       lty = c(1, 2, 4), lwd = c(2, 1.5, 1.5), bty = "n")
legend("topright",
       legend = c(paste0("RFM (Mean = ", round(mean(risk_error[1,]), 2), ")"), 
                  paste0("LFM (Mean = ", round(mean(risk_error[2,]), 2), ")"), 
                  paste0("LOCF (Mean = ", round(mean(risk_error[3,]), 2), ")")),
       col = c("lightblue", "lightsalmon", "darkseagreen1"),
       lty = c(1, 2, 4), lwd = 0.00,
       pch = c(19, 17, 23), bty = "n",
       pt.bg = c(NULL, NULL, "darkseagreen1"))











