library(tidyr)
library(ggplot2)
library(reshape2)
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
  
  if (type == 1) {
    mu = diag(rep(1, p))
  } else if (type == 2) {
    mu = 5 * toeplitz(0.6^c(0:(p - 1)))
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
    
    pct <- floor(10 * t / n)  # integer division gives 0 to 10
    if (pct > reported) {
      cat("  dta_gen:", paste0(pct * 10, "% complete\n"))
      reported <- pct
    }
  }
  
  return (list("X" = X, "X_nless" = X_nless, "Z" = Z, "Z_nless" = Z_nless, 
               "A" = A, "Factors" = Factors, "mu" = mu))
}

### Test: Do not run
# p = 10
# n = 200
# dta = dta_gen(n, "prototype", p = p)
# res = main_BWS(dta$X, 10, batch_size = 30, max.iter = 10, true_A = dta$A, 
#                true_mu = dta$mu, test_size = 0)
# round(res$FVU_RFM_BWS, 3)
# round(res$FVU_LYB_BWS, 3)
# round(res$FVU_RFM_Euc / res$FVU_LYB_Euc, 3)
# round(res$loading_dist[1:5], 3)
# res$r_hat_RFM
# res$r_hat_LYB


###########################
set.seed(5566)
num_sim = 100

ns = c(50, 100, 200)
p = 10
n_test = 200

for (type in types) {
  for (n in ns) {
    results = vector("list", length = num_sim)
    
    for (zz in 1:num_sim) {
      dta = dta_gen(n + n_test, type)
      res = main_BWS(dta$X, 10, n_test, true_A = dta$A, true_mu = dta$mu,
                     fraction = TRUE)
    }
  }
}


p = 10
n = 1000
dta = dta_gen(n, "prototype", p = p)
temp = mean_on_BWS(dta$X[1:200,,], batch_size = 30, max.iter = 100, tau = 0.5, verbose = T, tol = 0)
print(temp)
geod_BWS(temp, dta$mu)
mean(geod_BWS(dta$X, dta$mu))

z_log = array(NA, dim = c(n, p * (p + 1) / 2))
z_log_true = array(NA, dim = c(n, p * (p + 1) / 2))
for (m in 1:n) {
  v = Log_BWS(dta$X[m,,], temp)
  v_true = Log_BWS(dta$X[m,,], dta$mu)
  z_log[m,] = tangent_in_E(v, temp)
  z_log_true[m,] = tangent_in_E(v_true, dta$mu)
}
colMeans(z_log)
colMeans(z_log_true)

n = 300
p = 5
mean_dist1 = rep(0, 100)
loading_dist1 = rep(0, 100)
for (zz in 1:100) {
  dta = dta_gen(n, "prototype", p = p)
  res = main_BWS(dta$X, 10, batch_size = 20, max.iter = 50,
                 true_A = dta$A, true_mu = dta$mu, fraction = TRUE)
  mean_dist1[zz] = geod_BWS_core(res$mu_hat, dta$mu)
  loading_dist1[zz] = res$loading_dist[5]
  cat("iteration", zz, "\n")
}
n = 300
p = 5
mean_dist2 = rep(0, 100)
loading_dist2 = rep(0, 100)
for (zz in 1:100) {
  dta = dta_gen(n, 0, p = p)
  res = main_BWS(dta$X, 10, batch_size = 20, max.iter = 50,
                 true_A = dta$A, true_mu = dta$mu, fraction = TRUE)
  mean_dist2[zz] = geod_BWS_core(res$mu_hat, dta$mu)
  loading_dist2[zz] = res$loading_dist[5]
  cat("iteration", zz, "\n")
}
plot(x = mean_dist1, y = loading_dist1, col = "darkblue", 
     xlim = range(c(mean_dist1, mean_dist2)),
     ylim = range(c(loading_dist1, loading_dist2)))
points(x = mean_dist2, y = loading_dist2, col = "steelblue", pch = 19)
plot(x = mean_dist1, y = mean_dist2, xlim = c(0, 5), ylim = c(0, 5))
abline(a = 0, b = 1, col = "blue")
plot(x = loading_dist1, y = loading_dist2, xlim = c(0, 1), ylim = c(0, 1))
abline(a = 0, b = 1, col = "blue")


###########################
# Fix p = 10, r = 5, vary n = 50, 100, 200
#set.seed(5566)  
num_sim = 100
#type = 1

n = 50
p = 10
n_test = 200

results = vector("list", length = num_sim)

FVUs = array(0, dim = c(num_sim, 3, 10))
gFVUs = array(0, dim = c(num_sim, 3, 10))

PEs = array(0, dim = c(num_sim, 3, 10))
gPEs = array(0, dim = c(num_sim, 3, 10))
gfit = array(0, dim = c(num_sim, 10))      # only for RFM

all_gfit = array(NA, dim = c(num_sim, 3, 10))

for (zz in 1:num_sim) {
  dta = dta_gen(n + n_test, type)
  RFM = main_BWS(dta$X, 10, verbose = F, test_size = n_test, mu_tol = 1e-3, h = 6)
  
  results[[zz]] = RFM
  
  oracle_PE = 0
  for (i in 1:n_test) {
    oracle_PE = oracle_PE +
      sum((dta$X[n + i,,] - dta$X_nless[n + i,,])^2)
  }
  oracle_PE = sqrt(oracle_PE / n_test)
  PEs[zz,1,] = RFM$pe_e
  PEs[zz,2,] = RFM$pe_e_linear
  PEs[zz,3,] = oracle_PE
  
  oracle_gPE = 0
  for (i in 1:n_test) {
    oracle_gPE = oracle_gPE +
      geod_BWS_core(dta$X[n + i,,], dta$X_nless[n + i,,])^2
  }
  oracle_gPE = sqrt(oracle_gPE / n_test)
  gPEs[zz,1,] = RFM$pe_g
  gPEs[zz,2,] = RFM$pe_g_linear
  gPEs[zz,3,] = oracle_gPE

  oracle_FVU = 0
  for (i in 1:n) {
    oracle_FVU = oracle_FVU + 
      sum((dta$X[i,,] - dta$X_nless[i,,])^2)
  }
  
  oracle_gFVU = 0
  for (i in 1:n) {
    oracle_gFVU = oracle_gFVU + 
      geod_BWS_core(dta$X[i,,], dta$X_nless[i,,])^2
  }

  FVUs[zz,1,] = RFM$FVU_e
  FVUs[zz,2,] = RFM$FVU_e_linear
  FVUs[zz,3,] = oracle_FVU / RFM$TV_e
  
  gFVUs[zz,1,] = RFM$FVU_g
  gFVUs[zz,2,] = RFM$FVU_g_linear
  gFVUs[zz,3,] = oracle_gFVU / RFM$TV_g
  
  for (i in 1:10) {
    A_hat = RFM$V[,1:i]
    for (t in 1:n_test) {
      xhat = Exp_BWS_core(vector_to_symmetric(A_hat %*% t(A_hat) %*% dta$Z_nless[n + t,], p),
                          RFM$mu_hat)
      gfit[zz,i] = gfit[zz,i] + geod_BWS_core(xhat, dta$X_nless[n + t,,])^2
    }
  }
  if (zz == num_sim) {
    gfit = sqrt(gfit / n_test)
  }
  cat("iteration", zz, "\n")
}

all_gfit[,1,] = gfit

pdf(paste0("./save/fPE_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(PEs), 
            label.y = "Pseudo-prediction errors (Frobenius distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()
pdf(paste0("./save/gPE_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(gPEs), 
            label.y = "Pseudo-prediction errors (geodesic distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()
pdf(paste0("./save/fFVU_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = FVUs, 
            label.y = "Fraction of (Euclidean) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()
pdf(paste0("./save/gFVU_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(gFVUs), 
            label.y = "Fraction of (geodesic) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()


pdf(paste0("./save/fPE_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(PEs), 
            label.y = "Pseudo-prediction errors (Frobenius distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()
pdf(paste0("./save/gPE_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(gPEs), 
            label.y = "Pseudo-prediction errors (geodesic distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()
pdf(paste0("./save/fFVU_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = FVUs, 
            label.y = "Fraction of (Euclidean) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()
pdf(paste0("./save/gFVU_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(gFVUs), 
            label.y = "Fraction of (geodesic) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()

save(results, file = paste0("./save/results_n", n, "_type", type, ".RData"))

################################

n = 100
n_test = 200

results = vector("list", length = num_sim)

FVUs = array(0, dim = c(num_sim, 3, 10))
gFVUs = array(0, dim = c(num_sim, 3, 10))

PEs = array(0, dim = c(num_sim, 3, 10))
gPEs = array(0, dim = c(num_sim, 3, 10))
gfit = array(0, dim = c(num_sim, 10))      # only for RFM

for (zz in 1:num_sim) {
  dta = dta_gen(n + n_test, type)
  RFM = main_BWS(dta$X, 10, verbose = F, test_size = n_test, mu_tol = 1e-3, h = 6)
  
  results[[zz]] = RFM
  
  oracle_PE = 0
  for (i in 1:n_test) {
    oracle_PE = oracle_PE +
      sum((dta$X[n + i,,] - dta$X_nless[n + i,,])^2)
  }
  oracle_PE = sqrt(oracle_PE / n_test)
  PEs[zz,1,] = RFM$pe_e
  PEs[zz,2,] = RFM$pe_e_linear
  PEs[zz,3,] = oracle_PE
  
  oracle_gPE = 0
  for (i in 1:n_test) {
    oracle_gPE = oracle_gPE +
      geod_BWS_core(dta$X[n + i,,], dta$X_nless[n + i,,])^2
  }
  oracle_gPE = sqrt(oracle_gPE / n_test)
  gPEs[zz,1,] = RFM$pe_g
  gPEs[zz,2,] = RFM$pe_g_linear
  gPEs[zz,3,] = oracle_gPE
  
  oracle_FVU = 0
  for (i in 1:n) {
    oracle_FVU = oracle_FVU + 
      sum((dta$X[i,,] - dta$X_nless[i,,])^2)
  }
  
  oracle_gFVU = 0
  for (i in 1:n) {
    oracle_gFVU = oracle_gFVU + 
      geod_BWS_core(dta$X[i,,], dta$X_nless[i,,])^2
  }
  
  FVUs[zz,1,] = RFM$FVU_e
  FVUs[zz,2,] = RFM$FVU_e_linear
  FVUs[zz,3,] = oracle_FVU / RFM$TV_e
  
  gFVUs[zz,1,] = RFM$FVU_g
  gFVUs[zz,2,] = RFM$FVU_g_linear
  gFVUs[zz,3,] = oracle_gFVU / RFM$TV_g
  
  for (i in 1:10) {
    A_hat = RFM$V[,1:i]
    for (t in 1:n_test) {
      xhat = Exp_BWS_core(vector_to_symmetric(A_hat %*% t(A_hat) %*% dta$Z_nless[n + t,], p),
                          RFM$mu_hat)
      gfit[zz,i] = gfit[zz,i] + geod_BWS_core(xhat, dta$X_nless[n + t,,])^2
    }
  }
  if (zz == num_sim) {
    gfit = sqrt(gfit / n_test)
  }
  cat("iteration", zz, "\n")
}

all_gfit[,2,] = gfit

pdf(paste0("./save/fPE_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(PEs), 
            label.y = "Pseudo-prediction errors (Frobenius distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()
pdf(paste0("./save/gPE_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(gPEs), 
            label.y = "Pseudo-prediction errors (geodesic distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()
pdf(paste0("./save/fFVU_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = FVUs, 
            label.y = "Fraction of (Euclidean) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()
pdf(paste0("./save/gFVU_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(gFVUs), 
            label.y = "Fraction of (geodesic) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()


pdf(paste0("./save/fPE_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(PEs), 
            label.y = "Pseudo-prediction errors (Frobenius distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()
pdf(paste0("./save/gPE_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(gPEs), 
            label.y = "Pseudo-prediction errors (geodesic distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()
pdf(paste0("./save/fFVU_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = FVUs, 
            label.y = "Fraction of (Euclidean) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()
pdf(paste0("./save/gFVU_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(gFVUs), 
            label.y = "Fraction of (geodesic) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()

save(results, file = paste0("./save/results_n", n, "_type", type, ".RData"))

################################

n = 200
n_test = 200

results = vector("list", length = num_sim)

FVUs = array(0, dim = c(num_sim, 3, 10))
gFVUs = array(0, dim = c(num_sim, 3, 10))

PEs = array(0, dim = c(num_sim, 3, 10))
gPEs = array(0, dim = c(num_sim, 3, 10))
gfit = array(0, dim = c(num_sim, 10))      # only for RFM

for (zz in 1:num_sim) {
  dta = dta_gen(n + n_test, type)
  RFM = main_BWS(dta$X, 10, verbose = F, test_size = n_test, mu_tol = 1e-3, h = 6)
  
  results[[zz]] = RFM
  
  oracle_PE = 0
  for (i in 1:n_test) {
    oracle_PE = oracle_PE +
      sum((dta$X[n + i,,] - dta$X_nless[n + i,,])^2)
  }
  oracle_PE = sqrt(oracle_PE / n_test)
  PEs[zz,1,] = RFM$pe_e
  PEs[zz,2,] = RFM$pe_e_linear
  PEs[zz,3,] = oracle_PE
  
  oracle_gPE = 0
  for (i in 1:n_test) {
    oracle_gPE = oracle_gPE +
      geod_BWS_core(dta$X[n + i,,], dta$X_nless[n + i,,])^2
  }
  oracle_gPE = sqrt(oracle_gPE / n_test)
  gPEs[zz,1,] = RFM$pe_g
  gPEs[zz,2,] = RFM$pe_g_linear
  gPEs[zz,3,] = oracle_gPE
  
  oracle_FVU = 0
  for (i in 1:n) {
    oracle_FVU = oracle_FVU + 
      sum((dta$X[i,,] - dta$X_nless[i,,])^2)
  }
  
  oracle_gFVU = 0
  for (i in 1:n) {
    oracle_gFVU = oracle_gFVU + 
      geod_BWS_core(dta$X[i,,], dta$X_nless[i,,])^2
  }
  
  FVUs[zz,1,] = RFM$FVU_e
  FVUs[zz,2,] = RFM$FVU_e_linear
  FVUs[zz,3,] = oracle_FVU / RFM$TV_e
  
  gFVUs[zz,1,] = RFM$FVU_g
  gFVUs[zz,2,] = RFM$FVU_g_linear
  gFVUs[zz,3,] = oracle_gFVU / RFM$TV_g
  
  for (i in 1:10) {
    A_hat = RFM$V[,1:i]
    for (t in 1:n_test) {
      xhat = Exp_BWS_core(vector_to_symmetric(A_hat %*% t(A_hat) %*% dta$Z_nless[n + t,], p),
                          RFM$mu_hat)
      gfit[zz,i] = gfit[zz,i] + geod_BWS_core(xhat, dta$X_nless[n + t,,])^2
    }
  }
  if (zz == num_sim) {
    gfit = sqrt(gfit / n_test)
  }
  cat("iteration", zz, "\n")
}

all_gfit[,3,] = gfit

pdf(paste0("./save/fPE_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(PEs), 
            label.y = "Pseudo-prediction errors (Frobenius distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()
pdf(paste0("./save/gPE_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(gPEs), 
            label.y = "Pseudo-prediction errors (geodesic distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()
pdf(paste0("./save/fFVU_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = FVUs, 
            label.y = "Fraction of (Euclidean) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()
pdf(paste0("./save/gFVU_n", n, "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(gFVUs), 
            label.y = "Fraction of (geodesic) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = T)
dev.off()


pdf(paste0("./save/fPE_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(PEs), 
            label.y = "Pseudo-prediction errors (Frobenius distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()
pdf(paste0("./save/gPE_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(gPEs), 
            label.y = "Pseudo-prediction errors (geodesic distance)",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()
pdf(paste0("./save/fFVU_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = FVUs, 
            label.y = "Fraction of (Euclidean) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()
pdf(paste0("./save/gFVU_n", n, "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(gFVUs), 
            label.y = "Fraction of (geodesic) variance unexplained",
            Labels = c("Riemannian factor model", "Linear factor model", "Oracle"),
            legend.if = F)
dev.off()

save(results, file = paste0("./save/results_n", n, "_type", type, ".RData"))

################################
pdf(paste0("./save/gfit", "_type", type, ".pdf"), width = 8, height = 6)
plot_assist(x = Re(all_gfit), 
            label.y = "Goodness-of-fit (geodesic distance)",
            colors = colors()[c(468, 510, 552)],
            Labels = c("n = 50", "n = 100", "n = 200"),
            legend.if = T)
dev.off()
pdf(paste0("./save/gfit", "_type", type, "_nolegend.pdf"), width = 8, height = 6)
plot_assist(x = Re(all_gfit), 
            label.y = "Goodness-of-fit (geodesic distance)",
            colors = colors()[c(468, 510, 552)],
            Labels = c("n = 50", "n = 100", "n = 200"),
            legend.if = F)
dev.off()


