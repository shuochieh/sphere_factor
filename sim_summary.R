plot_assist = function (res1, res2 = NULL, oracle = NULL,
                        labs = NULL, ylim = NULL, main = NULL) {
  mean1 = colMeans(res1)
  q05_1 = apply(res1, 2, quantile, 0.05)
  q95_1 = apply(res1, 2, quantile, 0.95)
  
  if (!is.null(res2)) {
    mean2 = colMeans(res2)
    q05_2 = apply(res2, 2, quantile, 0.05)
    q95_2 = apply(res2, 2, quantile, 0.95)
  }
  
  if (is.null(ylim)) {
    ylim = c(0, 1)
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
  
  x_vals = 1:length(mean1)
  plot(x = x_vals, y = mean1, type = "n",
       ylim = ylim, 
       xlab = xlab, ylab = ylab, main = main, bty = "L")
  grid(col = "lightgray", lty = "dotted", lwd = 1)
  
  polygon(c(x_vals, rev(x_vals)), c(q05_1, rev(q95_1)),
          col = adjustcolor("lightblue", alpha.f = 0.6), border = NA)
  if (!is.null(res2)) {
    polygon(c(x_vals, rev(x_vals)), c(q05_2, rev(q95_2)),
            col = adjustcolor("lightsalmon", alpha.f = 0.6), border = NA)
  }
  lines(x_vals, mean1, col = "steelblue", lwd = 2)
  if (!is.null(res2)) {
    lines(x_vals, mean2, col = "firebrick", lwd = 2, lty = 2)
  }
  
  if (!is.null(oracle)) {
    abline(h = oracle, col = "black", lty = 4, lwd = 1.5)
  }
} 

cat_assist = function (x, type = 1) {
  res = NULL
  for (e in x) {
    if (type == 1) {
      res = c(res, paste(e, "&"))
    } else {
      res = c(res, paste0("(", e, ") & "))
    }
  }
  return (res)
}

par(mfrow = c(2, 3))
for (case in c(2, 3, 5, 6)) {
  for (p in c(5, 10)) {
    for (n in c(50, 100, 200)) {
      e <- new.env()
      load(paste0("./save/FVU_RFM_BWS_n", n, "_p", p, "_case", case, ".RData"), envir = e)
      FVU_RFM_BWS = Re(e$FVU_RFM_BWS)
      
      e <- new.env()
      load(paste0("./save/FVU_LYB_BWS_n", n, "_p", p, "_case", case, ".RData"), envir = e)
      FVU_LYB_BWS = Re(e$FVU_LYB_BWS)
      
      e <- new.env()
      load(paste0("./save/oracle_BWS_n", n, "_p", p, "_case", case, ".RData"), envir = e)
      oracle_BWS = Re(e$oracle_BWS)
      
      cat("Case", case, "n =", n, "p =", p, "\n")
      cat("Geodesic FVU (RFM):", cat_assist(round(colMeans(FVU_RFM_BWS) * 100, 1)), "\n")
      cat("Geodesic FVU (LYB):", cat_assist(round(colMeans(FVU_LYB_BWS) * 100, 1)), "\n")
      cat("Geodesic FVU (RFM) SD:", cat_assist(round(apply(FVU_RFM_BWS * 100, 2, sd), 1), 2), "\n")
      cat("Geodesic FVU (LYB) SD:", cat_assist(round(apply(FVU_LYB_BWS * 100, 2, sd), 1), 2), "\n")
      cat("Oracle:", round(mean(oracle_BWS * 100), 1), 
          "SD:", round(sd(oracle_BWS * 100), 1), "\n\n")
      
      plot_assist(FVU_RFM_BWS, res2 = FVU_LYB_BWS, oracle = mean(oracle_BWS), 
                  labs = c("number of factors", "rGMSE"), ylim = c(0.1, 1),
                  main = paste0("n = ", n, "; q = ", p))
    }
  }
}

dev.off()
par(mfrow = c(3, 3))
for (case in 1:6) {
  ns = c(50, 100, 200)
  ps = c(5, 10, 20)
  for (n_i in 1:3) {
    for (p_i in 1:3) {
      n = ns[n_i] ; p = ps[p_i]
      
      e <- new.env()
      load(paste0("./save/FVU_RFM_Euc_n", n, "_p", p, "_case", case, ".RData"), envir = e)
      FVU_RFM_Euc = Re(e$FVU_RFM_Euc)
      
      e <- new.env()
      load(paste0("./save/FVU_LYB_Euc_n", n, "_p", p, "_case", case, ".RData"), envir = e)
      FVU_LYB_Euc = Re(e$FVU_LYB_Euc)
      
      ratio = FVU_RFM_Euc / FVU_LYB_Euc
      
      cat("Case", case, "n =", n, "p =", p, "\n")
      cat("Ratio of Euclidean FVU:", round(colMeans(ratio), 2), "\n")
      cat("Ratio SD:", round(apply(ratio, 2, sd), 2), "\n\n")

      plot_assist(ratio, oracle = mean(oracle_BWS), 
                  labs = c("number of factors", "Euclidean FVU ratio"), ylim = c(0.5, 1.5),
                  main = paste0("n = ", n, "; q = ", p))
    }
  }
}

dev.off()
par(mfrow = c(2, 2))
for (case in c(2, 3, 5, 6)) {
  ns = c(50, 100, 200)
  ps = c(5, 10)
  for (p_i in 1:2) {
    loading_d_all = array(NA, dim = c(300, 3))
    for (n_i in 1:3) {
      n = ns[n_i] ; p = ps[p_i]
      
      e <- new.env()
      load(paste0("./save/loading_d_n", n, "_p", p, "_case", case, ".RData"), envir = e)
      loading_d = Re(e$loading_d)
      
      
      
      cat("Case", case, "n =", n, "p =", p, "\n")
      cat("Loading dist:", round(mean(loading_d[,5]), 2), "\n")
      cat("Loading dist SD:", round(sd(loading_d[,5]), 2), "\n\n")
      
      loading_d_all[,n_i] = loading_d[,5]
    }
    
    if (case %% 3 == 1) {
      alpha = 0.5
    } else if (case %% 3 == 2) {
      alpha = 0.2
    } else if (case %% 3 == 0) {
      alpha = 0.8
    }
    
    boxplot(loading_d_all, 
            names = c(paste0("n = ", 50), paste0("n = ", 100), paste0("n = ", 200)),
            main = paste0("alpha = ", alpha , "; q = ", p), 
            xlab = "", ylab = "", ylim = c(0.1, 0.7),
            col = c("skyblue", "brown", "lightgreen"))
  }
}

dev.off()
for (case in c(2, 3, 5, 6)) {
  for (p in c(5, 10)) {
    for (n in c(50, 100, 200)) {
      e <- new.env()
      load(paste0("./save/r_hat_RFM_n", n, "_p", p, "_case", case, ".RData"), envir = e)
      r_hat_RFM = Re(e$r_hat_RFM)
      
      e <- new.env()
      load(paste0("./save/r_hat_LYB_n", n, "_p", p, "_case", case, ".RData"), envir = e)
      r_hat_LYB = Re(e$r_hat_LYB)
      
      cat("Case", case, "n =", n, "p =", p, "\n")
      cat("RFM frequency of correct rank:", round(mean(r_hat_RFM == 5), 2), "\n")
      cat("LYB frequency of correct rank:", round(mean(r_hat_LYB == 5), 2), "\n\n")
    }
  }
}





