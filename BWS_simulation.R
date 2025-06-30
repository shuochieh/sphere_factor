library(foreach)
library(doParallel)
source("./BWS_util.R")
source("./main_func.R")

t_isonormal_sampler = function (n, p, norm_cut = Inf, sd = 1) {
  # draws truncated isotropic p-dimensional normal
  counter = 0
  res = array(NA, dim = c(n, p))
  while (counter < n) {
    n_sampler = n * 10
    draw = matrix(rnorm(n_sampler * p, sd = sd), nrow = n_sampler, ncol = p)
    idx_effective = which(sqrt(rowSums(draw^2)) < norm_cut)
    n_eff = length(idx_effective)
    if (n_eff == 0) {
      next
    }
    if (counter + n_eff < n) {
      res[(counter + 1):(counter + n_eff),] = draw[idx_effective,]
      counter = counter + n_eff
    } else {
      res[(counter + 1):n,] = draw[idx_effective[1:(n - counter)],]
      counter = counter + n_eff
    }
  }
  
  if (p > 1) {
    return (res)
  } else {
    return (c(res))
  } 
}

#' generate synthetic data
#' 
#' @param n sample size
#' @param mu_type type of base point (1 = identity, 2 = Toeplitz)
#' @param s controls size of factor
#' @param z_noise controls size of noise
#' @param alpha AR coefficient for factors
dta_gen_BWS = function (n, p, mu_type, r = 5, 
                        s = 1, z_noise = 1, alpha = 0.7) {
  
  if (mu_type == 1) {
    mu = diag(rep(1, p))
  } else if (mu_type == 2) {
    mu = 5 * toeplitz(0.8^c(0:(p - 1)))
  } else {
    stop("dta_gen_BWS: unsupported mu_type")
  }
  
  Factors = array(0, dim = c(n + 100, r))
  for (t in 2:(n + 100)) {
    Factors[t,] = alpha * Factors[t - 1,] + 
      s * t_isonormal_sampler(1, r, 1, sd = 1)
  }
  Factors = Factors[-c(1:100),]
  
  A = array(rnorm(p * (p + 1) * r / 2), dim = c(p * (p + 1) / 2, r))
  A = qr.Q(qr(A))
  
  Z = array(NA, dim = c(n, p * (p + 1) / 2))
  Z_nless = array(NA, dim = c(n, p * (p + 1) / 2))
  for (t in 1:n) {
    Z_nless[t,] = A %*% c(Factors[t,])
    pertb = rnorm(p * (p + 1) / 2)
    pertb = pertb / sqrt(sum(pertb^2))
    Z[t,] = A %*% c(Factors[t,]) + z_noise * pertb
  }
  
  coord = tan_basis_bws(mu)
  
  X = array(NA, dim = c(n, p, p))
  X_nless = array(NA, dim = c(n, p, p))
  
  reported = 0
  for (t in 1:n) {
    V = log_to_tangent(Z[t,], coord$E)
    X[t,,] = Exp_BWS_core(V, mu)
    
    V_nless = log_to_tangent(Z_nless[t,], coord$E)
    X_nless[t,,] = Exp_BWS_core(V_nless, mu)
    
    # pct <- floor(10 * t / n)  
    # if (pct > reported) {
    #   cat("  dta_gen:", paste0(pct * 10, "% complete\n"))
    #   reported <- pct
    # }
  }
  
  return (list("X" = X, "X_nless" = X_nless, "Z" = Z, "Z_nless" = Z_nless, 
               "A" = A, "Factors" = Factors, "mu" = mu))
}

### Test: Do not run
# p = 10
# n = 200
# dta = dta_gen(n, p = p)
# res = main_BWS(dta$X, 10, batch_size = 30, max.iter = 10, true_A = dta$A, 
#                true_mu = dta$mu, test_size = 0)
# round(res$FVU_RFM_BWS, 3)
# round(res$FVU_LYB_BWS, 3)
# round(res$FVU_RFM_Euc / res$FVU_LYB_Euc, 3)
# round(res$loading_dist[1:5], 3)
# res$r_hat_RFM
# res$r_hat_LYB

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


###########################
num_sim = 300

# CASE SWITCHING HELPER
# n = 50, p = 10
# Case 1: alpha = 0.5, z_noise = 1.5; mu_type = 1 (base line)
# Case 2: alpha = 0.2, z_noise = 1.5; mu_type = 1 (low signal)
# Case 3: alpha = 0.8, z_noise = 1.5; mu_type = 1 (high signal)
# Case 4: alpha = 0.5, z_noise = 1.5; mu_type = 2
# Case 5: alpha = 0.2, z_noise = 1.5; mu_type = 2
# Case 6: alpha = 0.8, z_noise = 1.5; mu_type = 2
case_param = function (case) {
  if (case == 1) {
    alpha = 0.5; z_noise = 1.0; s = 1.5; mu_type = 1
  } else if (case == 2) {
    alpha = 0.2; z_noise = 1.0; s = 1.5; mu_type = 1
  } else if (case == 3) {
    alpha = 0.8; z_noise = 1.0; s = 1.5; mu_type = 1
  } else if (case == 4) {
    alpha = 0.5; z_noise = 1.0; s = 1.5; mu_type = 2
  } else if (case == 5) {
    alpha = 0.2; z_noise = 1.0; s = 1.5; mu_type = 2
  } else if (case == 6) {
    alpha = 0.8; z_noise = 1.0; s = 1.5; mu_type = 2
  } else {
    stop("case_param: unsupported case")
  } 
  
  return (list("alpha" = alpha, "z_noise" = z_noise, "s" = s,
               "mu_type" = mu_type))
}

for (case in c(2, 3, 5, 6)) {
  cat("Case:", case, "...\n")
  parms = case_param(case)
  s = parms$s
  alpha = parms$alpha
  z_noise = parms$z_noise
  mu_type = parms$mu_type
  par(mfrow = c(2, 3))
  for (p in c(5, 10)) {
    for (n in c(50, 100, 200)) {
      n_cores = 12 #detectCores() - 1
      cl = makeCluster(n_cores)
      registerDoParallel(cl)
      
      results = foreach(i = 1:num_sim, .packages = c("maotai", "expm", "deSolve"),
                        .inorder = FALSE) %dopar% {
                          dta = dta_gen_BWS(n = n + 200, p = p, mu_type = mu_type, r = 5,
                                            s = s, z_noise = z_noise, alpha = alpha)
                          sim_res = main_BWS(dta$X, 10, test_size = 200, h = 6, batch_size = 16,
                                             max.iter = 16, true_A = dta$A, true_mu = dta$mu)
                          oracle_BWS = sum(geod_BWS(dta$X_nless, dta$X)^2) / sum(geod_BWS(dta$mu, dta$X)^2)
                          
                          sink()
                          cat("  iteration", i, "\n")
                          list("FVU_RFM_BWS" = sim_res$FVU_RFM_BWS, "FVU_LYB_BWS" = sim_res$FVU_LYB_BWS,
                               "FVU_RFM_Euc" = sim_res$FVU_RFM_Euc, "FVU_LYB_Euc" = sim_res$FVU_LYB_Euc,
                               "loading_d" = sim_res$loading_dist, 
                               "r_hat_RFM" = sim_res$r_hat_RFM, "r_hat_LYB" = sim_res$r_hat_LYB,
                               "oracle_BWS" = oracle_BWS)
                        }
      stopCluster(cl)
      
      FVU_RFM_BWS = array(NA, dim = c(num_sim, 10))
      FVU_LYB_BWS = array(NA, dim = c(num_sim, 10))
      FVU_RFM_Euc = array(NA, dim = c(num_sim, 10))
      FVU_LYB_Euc = array(NA, dim = c(num_sim, 10))
      loading_d = array(NA, dim = c(num_sim, 10))
      r_hat_RFM = rep(NA, num_sim)
      r_hat_LYB = rep(NA, num_sim)
      oracle_BWS = rep(NA, num_sim)
      
      for (i in 1:num_sim) {
        FVU_RFM_BWS[i,] = results[[i]]$FVU_RFM_BWS
        FVU_LYB_BWS[i,] = results[[i]]$FVU_LYB_BWS
        FVU_RFM_Euc[i,] = results[[i]]$FVU_RFM_Euc
        FVU_LYB_Euc[i,] = results[[i]]$FVU_LYB_Euc
        loading_d[i,] = results[[i]]$loading_d
        r_hat_RFM[i] = results[[i]]$r_hat_RFM
        r_hat_LYB[i] = results[[i]]$r_hat_LYB
        oracle_BWS[i] = results[[i]]$oracle_BWS
      }
      
      save(FVU_RFM_BWS, 
           file = paste0("./save/FVU_RFM_BWS_n", n, "_p", p, "_case", case, ".RData"))
      save(FVU_LYB_BWS, 
           file = paste0("./save/FVU_LYB_BWS_n", n, "_p", p, "_case", case, ".RData"))
      save(FVU_RFM_Euc, 
           file = paste0("./save/FVU_RFM_Euc_n", n, "_p", p, "_case", case, ".RData"))
      save(FVU_LYB_Euc, 
           file = paste0("./save/FVU_LYB_Euc_n", n, "_p", p, "_case", case, ".RData"))
      save(loading_d, 
           file = paste0("./save/loading_d_n", n, "_p", p, "_case", case, ".RData"))
      save(r_hat_RFM, 
           file = paste0("./save/r_hat_RFM_n", n, "_p", p, "_case", case, ".RData"))
      save(r_hat_LYB, 
           file = paste0("./save/r_hat_LYB_n", n, "_p", p, "_case", case, ".RData"))
      save(oracle_BWS,
           file = paste0("./save/oracle_BWS_n", n, "_p", p, "_case", case, ".RData"))
      
      
      plot_assist(FVU_RFM_BWS, res2 = FVU_LYB_BWS, oracle = mean(oracle_BWS), 
                  labs = c("number of factors", "Geodesic FVU"), 
                  main = paste0("n = ", n, "; q = ", p),
                  ylim = c(0.2, 1))
      
      cat("sub-case closed\n\n")
      
    }
  }
}
