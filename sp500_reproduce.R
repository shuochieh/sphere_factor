library(MASS)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(lubridate)
source("./main_func.R")
source("./BWS_util.R")

set.seed(1)

# Compute monthly VIX for later
vix = read.csv("./sp500_covariance/VIXCLS.csv")

vix_monthly = rep(NA, 240)
years = year(as.Date(vix[,1]))
months = month(as.Date(vix[,1]))
counter = 1
for (yr in 2000:2019) {
  for (m in 1:12) {
    idx = intersect(which(years == yr), which(months == m))
    vix_monthly[counter] = mean(vix[idx,2], na.rm = T)
    counter = counter + 1
  }
}

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

heat_plot = function (A, main, lims) {
  df = melt(A)
  colnames(df) = c("Row", "Column", "Value")
  
  df$Row <- factor(df$Row, levels = rev(rownames(A)))
  df$Column <- factor(df$Column, levels = colnames(A))
  
  mid_point = 0 # mean(lims)
  
  hm = ggplot(df, aes(Column, Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                         midpoint = mid_point, limits = lims) +
    # scale_fill_gradient2(low = "white", high = "darkred", limits = lims) +
    theme_minimal() +
    coord_fixed() +  
    labs(x = "", y = "", fill = "", 
         title = main) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return (hm)
}

dta = load("./sp500_covariance/sp500_12bySector.RData")
dta = covariances
rm("covariances")

temp = array(aperm(dta, c(2, 1, 3, 4)), 
             c(dim(dta)[1] * dim(dta)[2], dim(dta)[c(3,4)]))

dta = temp
rm("temp")

dta = dta * 10000 # Convert to percentage points
overall_covariance_training = overall_covariance_training * 10000
q = length(selected_companies)

# Use observations from 2000--2019
dta_dates = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month")[1:240]
dta = dta[1:240,,]

RFM_xhat = array(NA, dim = c(15, 36, q, q))
LFM_xhat = array(NA, dim = c(15, 36, q, q))
for (k in 1:15) {
  results = main_BWS(dta, r = k, test_size = 36, h = 6, batch_size = 30, max.iter = 100,
                     return_predictions = TRUE)
  RFM_xhat[k,,,] = results$RFM_xhat
  LFM_xhat[k,,,] = results$LYB_xhat
  
  
  cat("iteration", k, "\n")
}

par(mfrow = c(1, 2))
plot(x = dta_dates[1:(length(dta_dates) - 36)],
     y = results$Factors[,1],
     type = "l",
     xlab = "", ylab = "",
     col = "steelblue",
     main = "Factor 1")
year_starts = seq(from = as.Date(dta_dates[1]), to = as.Date(dta_dates[length(dta_dates) - 36]), by = "year")
abline(v = year_starts, col = "lightgray", lty = "dotted", lwd = 1)
abline(h = 0, col = "lightgray", lty = "dotted", lwd = 1)
legend("topleft", legend = c("Factor", "VIX closing price"),
       col = c("steelblue", "darkorange"),
       lty = c(1, 2), bty = "n")
par(new = TRUE)
plot(x = dta_dates[1:(length(dta_dates) - 36)],
     y = vix_monthly[1:(240 - 36)],
     type = "l",
     lty = 2,
     xlab = "", ylab = "",
     col = "darkorange",
     axes = FALSE)
axis(side = 4)
# mtext("VIX closing price", side = 4, line = 3)
plot(x = dta_dates[1:(length(dta_dates) - 36)],
     y = -results$Factors[,2],
     type = "l",
     xlab = "", ylab = "",
     col = "steelblue",
     main = "Factor 2")
year_starts = seq(from = as.Date(dta_dates[1]), to = as.Date(dta_dates[length(dta_dates) - 36]), by = "year")
abline(v = year_starts, col = "lightgray", lty = "dotted", lwd = 1)
abline(h = 0, col = "lightgray", lty = "dotted", lwd = 1)


ps = NULL
for (i in 1:14) {
  grid = seq(from = -6.0, to = 6.0, length.out = 7)
  if (i <= 7) {
    inc = grid[i]
    Exp_V = Exp_BWS(log_to_tangent(inc * results$V[,1], results$E), results$mu_hat)
  } else {
    inc = grid[i - 7]
    Exp_V = Exp_BWS(log_to_tangent(-inc * results$V[,2], results$E), results$mu_hat)
  }
  
  rownames(Exp_V) = selected_companies
  colnames(Exp_V) = selected_companies
  lim = 35
  if (i < 7) {
    ps[[i]] = heat_plot(Exp_V, paste("s =", round(inc,2)), lims = c(-lim, lim)) + 
      theme(legend.position = "none") +
      theme(plot.margin = unit(c(0.0, 0.2, 0.0, 0.0), "cm")) + 
      theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold"))
    
  } else if (i == 7) {
    ps[[i]] = heat_plot(Exp_V, paste("s =", round(inc,2)), lims = c(-lim, lim)) +
      theme(legend.position = "none") +
      theme(plot.margin = unit(c(0.0, 0.2, 0.0, 0.0), "cm")) + 
      theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold"))
  } else if (i < 14) {
    ps[[i]] = heat_plot(Exp_V, paste("s =", round(inc,2)), lims = c(-lim, lim)) + 
      theme(legend.position = "none") +
      theme(plot.margin = unit(c(0.0, 0.2, 0.0, 0.0), "cm")) + 
      theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold"))
  } else {
    ps[[i]] = heat_plot(Exp_V, paste("s =", round(inc,2)), lims = c(-lim, lim)) +
      theme(legend.position = "none") +
      theme(plot.margin = unit(c(0.0, 0.2, 0.0, 0.0), "cm")) + 
      theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 12, face = "bold"))
  }
}

grid.arrange(grobs = ps[1:7], nrow = 1, ncol = 7)
grid.arrange(grobs = ps[8:14], nrow = 1, ncol = 7)


### Out-of-sample factor forecasting

EWMA = function (dta, test_size = 36, lambda = 0.94) {
  res = array(NA, dim = dim(dta))
  res[1,,] = matrix(0, nrow = dim(dta)[2], ncol = dim(dta)[3])
  sigma_tilde = matrix(0, nrow = dim(dta)[2], ncol = dim(dta)[3])
  
  for (i in 2:(dim(dta)[1])) {
    sigma_tilde = lambda * sigma_tilde + (1 - lambda) * dta[i - 1,,]
    # sigma_tilde = ((lambda - lambda^i) / (1 - lambda^i)) * sigma_tilde + ((1 - lambda) / (1 - lambda^i)) * dta[i - 1,,]
    res[i,,] = sigma_tilde
  }
  
  return (res)
}

RFM_res = dyn_RFM(dta, r = 2, test_size = 36, h = 6, batch_size = 30)
LFM_res = dyn_LFM(dta, r = 1, test_size = 36, h = 6)
EWMA_res = EWMA(dta, 36, 0.94)


BWS_errors = array(NA, dim = c(4, 36))
Euc_errors = array(NA, dim = c(4, 36))
for (m in 1:36) {
  
  truth0 = dta[(length(dta_dates) - 36 + m),,]
  truth = eigen(truth0)
  
  for (k in 1:q) {
    # Riemannian factor model
    Sigma_hat = RFM_res[m,,]
    temp = eigen(Sigma_hat)
    v = as.matrix(temp$vectors[,1:k])
    if (k == 1) {
      BWS_errors[1,m] = geod_BWS_core(Sigma_hat, truth0)^2
      Euc_errors[1,m] = norm(Sigma_hat - truth0, "F")
    }
    
    # Linear factor model
    Sigma_hat = LFM_res[m,,]
    if (k == 1) {
      BWS_errors[2,m] = geod_BWS_core(project_to_SPD(Sigma_hat, 1e-6), truth0)^2
      Euc_errors[2,m] = norm(Sigma_hat - truth0, "F")
    }
    Sigma_hat = project_to_SPD(Sigma_hat, 1e-6)
    temp = eigen(Sigma_hat)
    v = as.matrix(temp$vectors[,1:k])

    # Random walk model
    Sigma_hat = dta[(length(dta_dates) - 36 + m - 1),,]
    temp = eigen(Sigma_hat)
    v = as.matrix(temp$vectors[,1:k])
    if (k == 1) {
      BWS_errors[3,m] = geod_BWS_core(Sigma_hat, truth0)^2
      Euc_errors[3,m] = norm(Sigma_hat - truth0, "F")
    }
    
    # EWMA
    Sigma_hat = EWMA_res[(length(dta_dates) - 36 + m),,]
    temp = eigen(Sigma_hat)
    v = as.matrix(temp$vectors[,1:k])
    if (k == 1) {
      BWS_errors[4,m] = geod_BWS_core(Sigma_hat, truth0)^2
      Euc_errors[4,m] = norm(Sigma_hat - truth0, "F")
    }
  }
}


par(mfrow = c(1, 2))
x_dates = tail(dta_dates, 36)
plot(x = x_dates, 
     y = sqrt(BWS_errors[1,]), type = "n",
     ylim = c(0, 6), 
     xlab = "", ylab = "BW distance", bty = "L")
grid(col = "lightgray", lty = "dotted", lwd = 1)
lines(x_dates, sqrt(BWS_errors[2,]), col = "firebrick", lwd = 1.5, lty = 2)
points(x_dates, sqrt(BWS_errors[2,]), pch = 17, col = "lightsalmon")
lines(x_dates, sqrt(BWS_errors[3,]), col = "darkseagreen4", lwd = 1.5, lty = 4)
points(x_dates, sqrt(BWS_errors[3,]), pch = 23, col = "darkseagreen1",
       bg = "darkseagreen1")
lines(x_dates, sqrt(BWS_errors[4,]), col = "mediumpurple4", lwd = 1.5, lty = 6)
points(x_dates, sqrt(BWS_errors[4,]), pch = 15, col = "mediumpurple1")
lines(x_dates, sqrt(BWS_errors[1,]), col = "steelblue", lwd = 2)
points(x_dates, sqrt(BWS_errors[1,]), pch = 19, col = "lightblue")

model_names <- c("RFM", "LFM", "LOCF", "EWMA")
legend_labels <- sapply(1:4, function(i) {
  sprintf("%s (Mean = %.2f; Median = %.2f)",
          model_names[i],
          round(rowMeans(sqrt(BWS_errors))[i], 2),
          round(median(sqrt(BWS_errors[i,])), 2))
})

legend(
  x = x_dates[1], 
  y = 1.25,
  legend = legend_labels,
  col = c("steelblue", "firebrick", "darkseagreen4", "mediumpurple4"), # Line colors
  lty = c(1, 2, 4, 6),
  lwd = c(2, 1.5, 1.5, 1.5),
  pch = NA,             # Hide points
  bty = "n",
  cex = 0.9,
  text.font = 2         # Bold font
)
legend(
  x = x_dates[1],       # <-- Identical x
  y = 1.25,               # <-- Identical y
  legend = legend_labels,   # <-- Identical labels (for spacing)
  col = c("lightblue", "lightsalmon", "darkseagreen1", "mediumpurple1"), # Point colors
  pt.bg = c(NA, NA, "darkseagreen1", NA),
  pch = c(19, 17, 23, 15),
  lty = 0,                # Hide lines
  bty = "n",
  cex = 0.9,            # <-- Identical cex
  text.col = "transparent" # <-- Makes text invisible
)

x_dates = tail(dta_dates, 36)
plot(x = x_dates, 
     y = Euc_errors[1,], type = "n",
     ylim = c(0, 60), 
     xlab = "", ylab = "Frobenius distance", bty = "L")
grid(col = "lightgray", lty = "dotted", lwd = 1)
lines(x_dates, Euc_errors[2,], col = "firebrick", lwd = 1.5, lty = 2)
points(x_dates, Euc_errors[2,], pch = 17, col = "lightsalmon")
lines(x_dates, Euc_errors[3,], col = "darkseagreen4", lwd = 1.5, lty = 4)
points(x_dates, Euc_errors[3,], pch = 23, col = "darkseagreen1",
       bg = "darkseagreen1")
lines(x_dates, Euc_errors[4,], col = "mediumpurple4", lwd = 1.5, lty = 6)
points(x_dates, Euc_errors[4,], pch = 15, col = "mediumpurple1")
lines(x_dates, Euc_errors[1,], col = "steelblue", lwd = 2)
points(x_dates, Euc_errors[1,], pch = 19, col = "lightblue")

model_names <- c("RFM", "LFM", "LOCF", "EWMA")
legend_labels <- sapply(1:4, function(i) {
  sprintf("%s (Mean = %.2f; Median = %.2f)",
          model_names[i],
          round(rowMeans(Euc_errors)[i], 2),
          round(median(Euc_errors[i,]), 2))
})

legend(
  x = x_dates[1], 
  y = 60,
  legend = legend_labels,
  col = c("steelblue", "firebrick", "darkseagreen4", "mediumpurple4"), # Line colors
  lty = c(1, 2, 4, 6),
  lwd = c(2, 1.5, 1.5, 1.5),
  pch = NA,             # Hide points
  bty = "n",
  cex = 0.9,
  text.font = 2         # Bold font
)

legend(
  x = x_dates[1],       # <-- Identical x
  y = 60,               # <-- Identical y
  legend = legend_labels,   # <-- Identical labels (for spacing)
  col = c("lightblue", "lightsalmon", "darkseagreen1", "mediumpurple1"), # Point colors
  pt.bg = c(NA, NA, "darkseagreen1", NA),
  pch = c(19, 17, 23, 15),
  lty = 0,                # Hide lines
  bty = "n",
  cex = 0.9,            # <-- Identical cex
  text.col = "transparent" # <-- Makes text invisible
)

risk_error = array(NA, dim = c(4, 36))
for (m in 1:36) {
  truth = dta[(length(dta_dates) - 36 + m),,]
  truth_lag = dta[(length(dta_dates) - 36 + m - 1),,]
  w = solve(truth_lag) %*% rep(1, q)
  w = w / sum(w)
  true_risk = c(t(w) %*% truth %*% w)
  
  # Riemannian factor model
  Sigma_hat = RFM_res[m,,]
  predicted_risk = c(t(w) %*% Sigma_hat %*% w)
  risk_error[1,m] = (predicted_risk - true_risk)^2
  
  # Linear factor model
  Sigma_hat = LFM_res[m,,] # project_to_SPD(LFM_res[m,,], 1e-6)
  predicted_risk = c(t(w) %*% Sigma_hat %*% w)
  risk_error[2,m] = (predicted_risk - true_risk)^2
  
  # Random walk model
  Sigma_hat = dta[(length(dta_dates) - 36 + m - 1),,]
  predicted_risk = c(t(w) %*% Sigma_hat %*% w)
  risk_error[3,m] = (predicted_risk - true_risk)^2
  
  # EWMA
  Sigma_hat = EWMA_res[(length(dta_dates) - 36 + m),,]
  predicted_risk = c(t(w) %*% Sigma_hat %*% w)
  risk_error[4,m] = (predicted_risk - true_risk)^2
}

risk_error = sqrt(risk_error)
x_dates = tail(dta_dates, 36)

par(mfrow = c(1, 1))
plot(x = x_dates,
     y = risk_error[1,],
     type = "n", ylim = c(0, 18),
     xlab = "", ylab = "risk prediction error", bty = "L")
grid(col = "lightgray", lty = "dotted", lwd = 1)
lines(x_dates, risk_error[2,], col = "firebrick", lwd = 1.5, lty = 2)
points(x_dates, risk_error[2,], col = "lightsalmon", pch = 17, cex = 0.8)
lines(x_dates, risk_error[3,], col = "darkseagreen4", lwd = 1.5, lty = 4)
points(x_dates, risk_error[3,], col = "darkseagreen1", bg = "darkseagreen1",
       pch = 23, cex = 0.8)
lines(x_dates, risk_error[4,], col = "mediumpurple4", lwd = 1.5, lty = 6)
points(x_dates, risk_error[4,], col = "mediumpurple1", pch = 15, cex = 0.8)
lines(x_dates, risk_error[1,], col = "steelblue", lwd = 2)
points(x_dates, risk_error[1,], col = "lightblue", pch = 19, cex = 0.8)

model_names <- c("RFM", "LFM", "LOCF", "EWMA")
legend_labels <- sapply(1:4, function(i) {
  sprintf("%s (Mean = %.2f; Median = %.2f)",
          model_names[i],
          round(mean(risk_error[i,]), 2),
          round(median(risk_error[i,]), 2))
})

legend(
  x = x_dates[1], 
  y = 18,
  legend = legend_labels,
  col = c("steelblue", "firebrick", "darkseagreen4", "mediumpurple4"), # Line colors
  lty = c(1, 2, 4, 6),
  lwd = c(2, 1.5, 1.5, 1.5),
  pch = NA,             # Hide points
  bty = "n",
  cex = 0.9,
  text.font = 2         # Bold font
)

# 3. Second Legend Call (Draws POINTS on top)
legend(
  x = x_dates[1],       # <-- Identical x
  y = 18,               # <-- Identical y
  legend = legend_labels,   # <-- Identical labels (for spacing)
  col = c("lightblue", "lightsalmon", "darkseagreen1", "mediumpurple1"), # Point colors
  pt.bg = c(NA, NA, "darkseagreen1", NA),
  pch = c(19, 17, 23, 15),
  lty = 0,                # Hide lines
  bty = "n",
  cex = 0.9,            # <-- Identical cex
  text.col = "transparent" # <-- Makes text invisible
)
