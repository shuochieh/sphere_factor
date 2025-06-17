library(foreach)
library(doParallel)

source("./sphere_util.R")
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

#' generate synthetic (product-)sphere-valued data
#' 
#' @param n sample size
#' @param qs a vector of sphere ambient dimensions
#' @param r number of factors
#' @param s controls size of factor
#' @param z_noise controls size of noise
#' @param alpha AR coefficient for factors
dta_gen_sphere = function (n, qs, r = 5, s = 1, z_noise = 1, alpha = 0.8) {
  
  # if (r > min(qs) - 1) {
  #   stop("dta_gen_sphere: does not support r > min(qs) - 1")
  # }
  
  d = length(qs)
  mus = vector("list", d)
  for (j in 1:d) {
    temp = rnorm(qs[j])
    mus[[j]] = temp / norm(temp, "2")
  }
  
  Factors = array(0, dim = c(n + 100, r))
  for (t in 2:(n + 100)) {
    Factors[t,] = alpha * Factors[t - 1,] + s * t_isonormal_sampler(1, r, 1, sd = 1)
  }
  Factors = Factors[-c(1:100),]
  
  A = array(rnorm((sum(qs) - d) * r), dim = c(sum(qs) - d, r))
  A = qr.Q(qr(A))
  
  Z = array(NA, dim = c(n, sum(qs) - d))
  Z_nless = array(NA, dim = c(n, sum(qs) - d))
  
  Z_nless = Factors %*% t(A)
  pertb = array(rnorm(n * (sum(qs) - d)), dim = c(n, sum(qs) - d))
  for (j in 1:d) {
    idx = q_index(qs, j, "intrinsic")
    pertb[, idx] = pertb[, idx] / apply(pertb[, idx], 1, norm, "2")
  }
  # pertb = pertb / apply(pertb, 1, norm, "2")
  Z = Z_nless + z_noise * pertb
  
  E = vector("list", d)
  for (j in 1:d) {
    E[[j]] = tan_basis_sphere(mus[[j]])
  }
  
  X_nless = vector("list", d)
  X = vector("list", d)
  
  for (j in 1:d) {
    idx2 = q_index(qs, j, "intrinsic")
    
    V_nl = Z_nless[,idx2] %*% t(E[[j]]) # n by q_j
    X_nless[[j]] = Exp_sphere(V_nl, mus[[j]])
    
    V = Z[,idx2] %*% t(E[[j]])
    X[[j]] = Exp_sphere(V, mus[[j]])
  }
  
  return (list("X" = X, "X_nless" = X_nless, "Z" = Z, "Z_nless" = Z_nless,
               "A" = A, "Factors" = Factors, "mus" = mus, "E" = E))
  
}

### Test: Do not run
# d = 4
# qs = rep(20, d)
# n = 200
# dta = dta_gen_sphere(n, qs, r = 5, s = 1.5, z_noise = 0.5, alpha = 0.8)
# res = main_sphere(dta$X, 10, max.iter = 100, tau = 0.1, true_A = dta$A,
#                true_mu = dta$mu, true_E = dta$E, test_size = 50)
# round(res$FVU_RFM_Sphere, 3)
# round(res$FVU_LYB_Sphere, 3)
# round(res$FVU_RFM_Euc, 3)
# round(res$FVU_LYB_Euc, 3)
# # round(res$FVU_RFM_Euc / res$FVU_LYB_Euc, 3)
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
num_sim = 100

# CASE SWITCHING HELPER
# Case 1: alpha = 0.5, z_noise = 0.5 (base line)
# Case 2: alpha = 0.2, z_noise = 0.5 (low signal)
# Case 3: alpha = 0.8, z_noise = 0.5 (high signal)
case_param = function (case) {
  if (case == 1) {
    alpha = 0.5; z_noise = 0.5; s = 2.5
  } else if (case == 2) {
    alpha = 0.2; z_noise = 0.5; s = 3 * sqrt(1 - alpha^2)
  } else if (case == 3) {
    alpha = 0.8; z_noise = 0.5; s = 3 * sqrt(1 - alpha^2)
  } else {
    stop("case_param: unsupported case")
  } 
  
  return (list("alpha" = alpha, "z_noise" = z_noise, "s" = s))
}

# ORACLE HELPER
oracle_geod = function (x_nless, x, mu) {
  res = 0
  total_var = 0
  d = length(x_nless)
  for (j in 1:d) {
    res = res + sum(geod_sphere(x_nless[[j]], x[[j]])^2)
    total_var = total_var + sum(geod_sphere(mu[[j]], x[[j]])^2)
  }
  
  return (res / total_var)
}

for (case in 1:3) {
  cat("Case:", case, "...\n")
  parms = case_param(case)
  s = parms$s
  alpha = parms$alpha
  z_noise = parms$z_noise
  par(mfrow = c(3, 3))
  for (d in c(5, 10, 20)) {
    for (n in c(50, 100, 200)) {
      n_cores = 12 #detectCores() - 1
      cl = makeCluster(n_cores)
      registerDoParallel(cl)
      
      results = foreach(i = 1:num_sim, .packages = c("Matrix"),
                        .inorder = FALSE) %dopar% {
                          dta = dta_gen_sphere(n = n + 200, qs = rep(5, d), 
                                               r = 5, s = s * sqrt(d / 5), z_noise = z_noise, 
                                               alpha = alpha)
                          sim_res = main_sphere(dta$X, 10, test_size = 200, h = 6, tau = 0.1,
                                                max.iter = 100, true_A = dta$A, true_mu = dta$mu,
                                                true_E = dta$E, fraction = TRUE)
                          oracle_Sphere = oracle_geod(dta$X_nless, dta$X, sim_res$mu_hat)
                          
                          sink()
                          cat("  iteration", i, "\n")
                          list("FVU_RFM_Sphere" = sim_res$FVU_RFM_Sphere, 
                               "FVU_LYB_Sphere" = sim_res$FVU_LYB_Sphere,
                               "FVU_RFM_Euc" = sim_res$FVU_RFM_Euc, 
                               "FVU_LYB_Euc" = sim_res$FVU_LYB_Euc,
                               "loading_d" = sim_res$loading_dist, 
                               "r_hat_RFM" = sim_res$r_hat_RFM, "r_hat_LYB" = sim_res$r_hat_LYB,
                               "oracle_Sphere" = oracle_Sphere)
                        }
      stopCluster(cl)
      
      FVU_RFM_Sphere = array(NA, dim = c(num_sim, 10))
      FVU_LYB_Sphere = array(NA, dim = c(num_sim, 10))
      FVU_RFM_Euc = array(NA, dim = c(num_sim, 10))
      FVU_LYB_Euc = array(NA, dim = c(num_sim, 10))
      loading_d = array(NA, dim = c(num_sim, 10))
      r_hat_RFM = rep(NA, num_sim)
      r_hat_LYB = rep(NA, num_sim)
      oracle_Sphere = rep(NA, num_sim)
      
      for (i in 1:num_sim) {
        FVU_RFM_Sphere[i,] = results[[i]]$FVU_RFM_Sphere
        FVU_LYB_Sphere[i,] = results[[i]]$FVU_LYB_Sphere
        FVU_RFM_Euc[i,] = results[[i]]$FVU_RFM_Euc
        FVU_LYB_Euc[i,] = results[[i]]$FVU_LYB_Euc
        loading_d[i,] = results[[i]]$loading_d
        r_hat_RFM[i] = results[[i]]$r_hat_RFM
        r_hat_LYB[i] = results[[i]]$r_hat_LYB
        oracle_Sphere[i] = results[[i]]$oracle_Sphere
      }
      
      save(FVU_RFM_Sphere, 
           file = paste0("./save/FVU_RFM_Sphere_n", n, "_d", d, "_case", case, ".RData"))
      save(FVU_LYB_Sphere, 
           file = paste0("./save/FVU_LYB_Sphere_n", n, "_d", d, "_case", case, ".RData"))
      save(FVU_RFM_Euc, 
           file = paste0("./save/FVU_RFM_Euc_Sphere_n", n, "_d", d, "_case", case, ".RData"))
      save(FVU_LYB_Euc, 
           file = paste0("./save/FVU_LYB_Euc_Sphere_n", n, "_d", d, "_case", case, ".RData"))
      save(loading_d, 
           file = paste0("./save/loading_d_Sphere_n", n, "_d", d, "_case", case, ".RData"))
      save(r_hat_RFM, 
           file = paste0("./save/r_hat_RFM_Sphere_n", n, "_d", d, "_case", case, ".RData"))
      save(r_hat_LYB, 
           file = paste0("./save/r_hat_LYB_Sphere_n", n, "_d", d, "_case", case, ".RData"))
      save(oracle_Sphere,
           file = paste0("./save/oracle_Sphere_n", n, "_d", d, "_case", case, ".RData"))
      
      
      plot_assist(FVU_RFM_Euc, res2 = FVU_LYB_Euc, # oracle = mean(oracle_Sphere), 
                  labs = c("number of factors", "Euclidean FVU"), 
                  main = paste0("n = ", n, "; d = ", d))
      
      cat("sub-case closed\n\n")
      
    }
  }
}















