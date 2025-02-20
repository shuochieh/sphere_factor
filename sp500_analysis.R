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

mu = mean_on_BWS(dta, tau = 0.1, verbose = TRUE)

tsp = plot_time_series(x = data.frame(Time = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month"),
                                      "Geodesic distance to Frechet mean" = geod_BWS(dta, mu),
                                      check.names = F),
                       time = "Time",
                       series_cols = "Geodesic distance to Frechet mean",
                       title = "Geodesic distance to the Frechet mean",
                       y_label = "",
                       x_label = "Time",
                       l_size = 0.7,
                       p_size = 0,
                       legend.pos = "none",
                       y_lim = c(0, 20))
print(tsp)

# plot(x = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month"),
#      y = geod_BWS(dta, mu), 
#      type = "l", ylab = "", xlab="", 
#      main = "Geodesic distance to the Frechet mean")

results = main_BWS(dta, 20, 30)

tsp = plot_time_series(x = data.frame(Time = c(1:15),
                                      "Riemannian Factor Model" = results$FVU_e[1:15],
                                      "Linear Factor Model" = results$FVU_e_linea[1:15],
                                      check.names = F),
                       time = "Time",
                       series_cols = c("Riemannian Factor Model", "Linear Factor Model"),
                       title = "",
                       y_label = "Fraction of (Euclidean) variations unexplained",
                       x_label = "Number of factors",
                       l_size = 0.5,
                       p_size = 2.5,
                       legend_title = "",
                       y_lim = c(0, 0.55))
print(tsp)

tsp = plot_time_series(x = data.frame(Time = c(1:15),
                                      "Riemannian Factor Model" = results$pe_e[1:15],
                                      "Linear Factor Model" = results$pe_e_linear[1:15],
                                      check.names = F),
                       time = "Time",
                       series_cols = c("Riemannian Factor Model", "Linear Factor Model"),
                       title = "",
                       y_label = "(Euclidean) Prediction errors",
                       x_label = "Number of factors",
                       l_size = 0.5,
                       p_size = 2.5,
                       legend_title = "",
                       y_lim = c(0, 16))
print(tsp)

tsp = plot_time_series(x = data.frame(Time = seq(from = as.Date("2000-01-01"), to = as.Date("2024-12-01"), by = "month")[1:270],
                                      "1st Factor" = results$Factors[,1],
                                      "2nd Factor" = results$Factors[,2],
                                      check.names = F),
                       time = "Time",
                       series_cols = c("1st Factor", "2nd Factor"),
                       title = "",
                       y_label = "",
                       x_label = "Time",
                       l_size = 0.7,
                       p_size = 0,
                       legend = "top",
                       legend_title = "")
print(tsp)


for (i in 1:2) {
  A = vector_to_symmetric(results$V[,i], 13)
  A = A[c(1, 2, 3, 9, 4, 13, 5, 8, 11, 6, 7, 10, 12),] 
  A = A[,c(1, 2, 3, 9, 4, 13, 5, 8, 11, 6, 7, 10, 12)]
  rownames(A) = selected_companies[c(1, 2, 3, 9, 4, 13, 5, 8, 11, 6, 7, 10, 12)]
  colnames(A) = selected_companies[c(1, 2, 3, 9, 4, 13, 5, 8, 11, 6, 7, 10, 12)]
  
  df = melt(A)
  colnames(df) = c("Row", "Column", "Value")
  
  df$Row <- factor(df$Row, levels = rev(rownames(A)))
  df$Column <- factor(df$Column, levels = colnames(A))
  
  hm = ggplot(df, aes(Column, Row, fill = Value)) +
    geom_tile() +
    # scale_fill_distiller(palette = "RdBu", limits = c(-0.5, 0.5), direction = 1) +
    # scale_fill_gradient(low = "black", high = "white", limits = c(-0.5, 0.4)) +  # Grayscale
    scale_fill_gradient2(low = "#92C5DE", mid = "white", high = "red", 
                         midpoint = 0.0, limits = c(-0.5, 0.4)) +
    theme_minimal() +
    coord_fixed() +  # Ensure square tiles
    labs(x = "", y = "", fill = "Loading", 
         title = paste("Loading matrix", i)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(hm)
}

for (i in 1:3) {
  A = abs(vector_to_symmetric(results$V[,i], 13))
  A = A[c(1, 8, 9, 10, 2, 5, 3, 6, 4, 12, 7, 11, 13),] 
  A = A[,c(1, 8, 9, 10, 2, 5, 3, 6, 4, 12, 7, 11, 13)]
  rownames(A) = selected_companies[c(1, 8, 9, 10, 2, 5, 3, 6, 4, 12, 7, 11, 13)]
  colnames(A) = selected_companies[c(1, 8, 9, 10, 2, 5, 3, 6, 4, 12, 7, 11, 13)]
  
  df = melt(A)
  colnames(df) = c("Row", "Column", "Value")
  
  df$Row <- factor(df$Row, levels = rev(rownames(A)))
  df$Column <- factor(df$Column, levels = colnames(A))
  
  hm = ggplot(df, aes(Column, Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black", limits = c(0.0, 0.5)) +  # Grayscale
    theme_minimal() +
    coord_fixed() +  # Ensure square tiles
    labs(x = "", y = "", fill = "Loading", 
         title = paste("Loading matrix", i, "(in absolute value)")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(hm)
}

# Baselines
proj_SPD = function (A) {
  if (!isSymmetric(A)) {
    stop("proj_SPD: A is not symmetric")
  }
  
  temp = eigen(A)
  pos.ev = temp$values[which(temp$values > 0)]
  V = temp$vectors[,which(temp$values > 0)]
  
  return (V %*% diag(pos.ev) %*% t(V))
}

diag_loss = array(0, dim = c(2, 15))
offdiag_loss = array(0, dim = c(2, 15))
overall_loss = array(0, dim = c(2, 15))

overall_cov_loss = 0
random_walk_loss = 0

for (i in 1:30) {
  truth = dta[dim(dta)[1] - 30 + i,,]
  for (k in 1:15) {
    temp = results$x_hat_RFM[i,k,,]
    offdiag_loss[1,k] = offdiag_loss[1,k] + sum((truth[upper.tri(truth)] - temp[upper.tri(temp)])^2)
    diag_loss[1,k] = diag_loss[1,k] + sum((diag(temp) - diag(truth))^2)
    overall_loss[1,k] = overall_loss[1,k] + norm(temp - truth, "F")^2
    
    temp = proj_SPD(results$x_hat_LFM[i,k,,])
    offdiag_loss[2,k] = offdiag_loss[2,k] + sum((truth[upper.tri(truth)] - temp[upper.tri(temp)])^2)
    diag_loss[2,k] = diag_loss[2,k] + sum((diag(temp) - diag(truth))^2)
    overall_loss[2,k] = overall_loss[2,k] + norm(temp - truth, "F")^2
  }
  
  overall_cov_loss = overall_cov_loss + norm(overall_covariance_training - truth, "F")^2
  rw = dta[dim(dta)[1] - 30 + i - 1,,]
  random_walk_loss = random_walk_loss + norm(rw - truth, "F")^2
}
diag_loss = sqrt(diag_loss / 30)
offdiag_loss = sqrt(offdiag_loss / 30)
overall_loss = sqrt(overall_loss / 30)
overall_cov_loss = rep(sqrt(overall_cov_loss / 30), 15)
random_walk_loss = rep(sqrt(random_walk_loss / 30), 15)



tsp = plot_time_series(x = data.frame(Time = c(1:15),
                                      "Riemannian Factor Model" = diag_loss[1,],
                                      "Linear Factor Model (after projection)" = diag_loss[2,],
                                      check.names = F),
                       time = "Time",
                       series_cols = c("Riemannian Factor Model", "Linear Factor Model (after projection)"),
                       title = "",
                       y_label = "Prediction errors of diagonal elements",
                       x_label = "Number of factors",
                       l_size = 0.5,
                       p_size = 2.5,
                       legend_title = "",
                       y_lim = c(0, 9))
print(tsp)

tsp = plot_time_series(x = data.frame(Time = c(1:15),
                                      "Riemannian Factor Model" = offdiag_loss[1,],
                                      "Linear Factor Model (after projection)" = offdiag_loss[2,],
                                      check.names = F),
                       time = "Time",
                       series_cols = c("Riemannian Factor Model", "Linear Factor Model (after projection)"),
                       title = "",
                       y_label = "Prediction errors of off-diagonal elements",
                       x_label = "Number of factors",
                       l_size = 0.5,
                       p_size = 2.5,
                       legend_title = "",
                       y_lim = c(0, 9))
print(tsp)

tsp = plot_time_series(x = data.frame(Time = c(1:15),
                                      "Riemannian Factor Model" = overall_loss[1,],
                                      "Linear Factor Model" = results$pe_e_linear[1:15],
                                      "Linear Factor Model (after projection)" = overall_loss[2,],
                                      "Random Walk" = random_walk_loss,
                                      "Constant Covariance" = overall_cov_loss,
                                      check.names = F),
                       time = "Time",
                       series_cols = c("Riemannian Factor Model", "Linear Factor Model", "Linear Factor Model (after projection)",
                                       "Random Walk", "Constant Covariance"),
                       title = "",
                       y_label = "Overall prediction errors",
                       x_label = "Number of factors",
                       l_size = 0.5,
                       p_size = 2.5,
                       legend_title = "",
                       legend.pos = "top",
                       legend_rows = 2,
                       y_lim = c(0, 25))
print(tsp)






















