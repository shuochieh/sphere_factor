
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

par(mfcol = c(3, 3))
for (case in c(1:1)) {
  for (n in c(50, 100, 200)) {
    for (p in c(5, 10, 20)) {
      e <- new.env()
      load(paste0("./save/FVU_RFM_Sphere_n", n, "_d", p, "_case", case, ".RData"), envir = e)
      FVU_RFM_Sphere = Re(e$FVU_RFM_Sphere)
      
      e <- new.env()
      load(paste0("./save/FVU_LYB_Sphere_n", n, "_d", p, "_case", case, ".RData"), envir = e)
      FVU_LYB_Sphere = Re(e$FVU_LYB_Sphere)
      
      e <- new.env()
      load(paste0("./save/oracle_Sphere_n", n, "_d", p, "_case", case, ".RData"), envir = e)
      oracle_Sphere = Re(e$oracle_Sphere)
      
      cat("Case", case, "n =", n, "p =", p, "\n")
      cat("Geodesic FVU (RFM):", cat_assist(round(colMeans(FVU_RFM_Sphere) * 100, 1)), "\n")
      cat("Geodesic FVU (LYB):", cat_assist(round(colMeans(FVU_LYB_Sphere) * 100, 1)), "\n")
      cat("Geodesic FVU (RFM) SD:", cat_assist(round(apply(FVU_RFM_Sphere * 100, 2, sd), 1), 2), "\n")
      cat("Geodesic FVU (LYB) SD:", cat_assist(round(apply(FVU_LYB_Sphere * 100, 2, sd), 1), 2), "\n")
      cat("Oracle:", round(mean(oracle_Sphere * 100), 1), 
          "SD:", round(sd(oracle_Sphere * 100), 1), "\n\n")
      
      plot_assist(FVU_RFM_Sphere, res2 = FVU_LYB_Sphere, oracle = mean(oracle_Sphere), 
                  labs = c("number of factors", "rMSE"), ylim = c(0.0, 1),
                  main = paste0("n = ", n, "; q = ", p))
    }
  }
}

dev.off()
par(mfrow = c(3, 3))
for (case in 1:1) {
  ns = c(50, 100, 200)
  ps = c(5, 10, 20)
  for (n_i in 1:3) {
    for (p_i in 1:length(ps)) {
      n = ns[n_i] ; p = ps[p_i]
      
      e <- new.env()
      load(paste0("./save/FVU_RFM_Euc_Sphere_n", n, "_d", p, "_case", case, ".RData"), envir = e)
      FVU_RFM_Euc = Re(e$FVU_RFM_Euc)
      
      e <- new.env()
      load(paste0("./save/FVU_LYB_Euc_Sphere_n", n, "_d", p, "_case", case, ".RData"), envir = e)
      FVU_LYB_Euc = Re(e$FVU_LYB_Euc)
      
      plot_assist(FVU_RFM_Euc, res2 = FVU_LYB_Euc, 
                  labs = c("number of factors", "Euclidean FVU"), ylim = c(0.1, 1),
                  main = paste0("n = ", n, "; d = ", p))
      
      
      ratio = FVU_RFM_Euc / FVU_LYB_Euc
      
      cat("Case", case, "n =", n, "d =", p, "\n")
      cat("Euclidean FVU (RFM):", cat_assist(round(colMeans(FVU_RFM_Euc) * 100, 1)), "\n")
      cat("Euclidean FVU (LYB):", cat_assist(round(colMeans(FVU_LYB_Euc) * 100, 1)), "\n")
      cat("Euclidean FVU (RFM) SD:", cat_assist(round(apply(FVU_RFM_Euc * 100, 2, sd), 1), 2), "\n")
      cat("Euclidean FVU (LYB) SD:", cat_assist(round(apply(FVU_LYB_Euc * 100, 2, sd), 1), 2), "\n")
      cat("Ratio of Euclidean FVU:", round(colMeans(ratio), 2), "\n\n")
      # cat("Ratio SD:", round(apply(ratio, 2, sd), 2), "\n\n")
      
      # plot_assist(ratio, oracle = mean(oracle_BWS), 
      #             labs = c("number of factors", "Euclidean FVU ratio"), ylim = c(0.5, 1.5),
      #             main = paste0("n = ", n, "; q = ", p))
    }
  }
}

dev.off()
par(mfrow = c(1, 3))
for (case in 1:1) {
  ns = c(50, 100, 200)
  ps = c(5, 10, 20)
  for (p_i in 1:3) {
    loading_d_all = array(NA, dim = c(100, 3))
    for (n_i in 1:3) {
      n = ns[n_i] ; p = ps[p_i]
      
      e <- new.env()
      load(paste0("./save/loading_d_Sphere_n", n, "_d", p, "_case", case, ".RData"), envir = e)
      loading_d = Re(e$loading_d)
      
      
      
      cat("Case", case, "n =", n, "d =", p, "\n")
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
            main = paste0("alpha = ", alpha , "; d = ", p), 
            xlab = "", ylab = "", ylim = c(0.1, 0.7),
            col = c("skyblue", "brown", "lightgreen"))
  }
}

dev.off()
for (case in 1:1) {
  for (p in c(5, 10, 20)) {
    for (n in c(50, 100, 200)) {
      e <- new.env()
      load(paste0("./save/r_hat_RFM_Sphere_n", n, "_d", p, "_case", case, ".RData"), envir = e)
      r_hat_RFM = Re(e$r_hat_RFM)
      
      e <- new.env()
      load(paste0("./save/r_hat_LYB_Sphere_n", n, "_d", p, "_case", case, ".RData"), envir = e)
      r_hat_LYB = Re(e$r_hat_LYB)
      
      cat("Case", case, "n =", n, "d =", p, "\n")
      cat("RFM frequency of correct rank:", round(mean(r_hat_RFM == 5), 2), "\n")
      cat("LYB frequency of correct rank:", round(mean(r_hat_LYB == 5), 2), "\n\n")
    }
  }
}





