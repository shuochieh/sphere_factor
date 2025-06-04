library(tidyr)
library(ggplot2)
library(reshape2)
source("./main_func.R")

dta_gen = function (n, type) {
  # type: specifications

  if (type == 1) {
    # BWS space p by p SPD matrices
    manifold.type = "BWS"
    p = 10
    r = 5
    fac_noise = sqrt(15)
    fac_sd = 2.0         # SNR~3.73
    alpha = 0.8          # AR coefficient for factor process
    
    mu = diag(5, p)
  }

  if (type == 2) {
    # BWS space p by p SPD matrices
    manifold.type = "BWS"
    p = 10
    r = 5
    fac_noise = sqrt(10)
    fac_sd = 2.0         # SNR~5.56
    alpha = 0.8          # AR coefficient for factor process
    
    mu = diag(5, p)
  }

  if (type == 3) {
    # BWS space p by p SPD matrices
    manifold.type = "BWS"
    p = 10
    r = 5
    fac_noise = sqrt(5)
    fac_sd = 2.0         # SNR~11.11
    alpha = 0.8          # AR coefficient for factor process
    
    mu = diag(5, p)
  }
  
  if (type == 4) {
    # BWS space p by p SPD matrices
    manifold.type = "BWS"
    p = 10
    r = 5
    fac_noise = sqrt(15)
    fac_sd = 2.0         # SNR~3.73
    alpha = 0.8 # AR coefficient for factor process
    
    mu = 10 * toeplitz(0.6^c(0:(p - 1)))
  }
  
  if (type == 5) {
    # BWS space p by p SPD matrices
    manifold.type = "BWS"
    p = 10
    r = 5
    fac_noise = sqrt(10)
    fac_sd = 2.0         # SNR~5.56
    alpha = 0.8 # AR coefficient for factor process
    
    mu = 10 * toeplitz(0.6^c(0:(p - 1)))
  }

  if (type == 6) {
    # BWS space p by p SPD matrices
    manifold.type = "BWS"
    p = 10
    r = 5
    fac_noise = sqrt(5)
    fac_sd = 2.0         # SNR~11.11
    alpha = 0.8 # AR coefficient for factor process
    
    mu = 10 * toeplitz(0.6^c(0:(p - 1)))
  }
  
  if (manifold.type == "BWS") {

    Factors = array(0, dim = c(n + 100, r))
    for (t in 2:(n + 100)) {
      Factors[t,] = alpha * Factors[t - 1,] + rnorm(r, sd = fac_sd)
    }
    Factors = Factors[-c(1:100),]
    
    A = array(rnorm(p * (p + 1) * r / 2), dim = c(p * (p + 1) / 2, r))
    A = qr.Q(qr(A))
    # A = t(t(A) / apply(A, 2, norm, "2"))

    Z = array(NA, dim = c(n, p * (p + 1) / 2))
    Z_nless = array(NA, dim = c(n, p * (p + 1) / 2))
    for (t in 1:n) {
      Z_nless[t,] = A %*% c(Factors[t,])
      Z[t,] = A %*% c(Factors[t,]) + rnorm(p * (p + 1) / 2, 
                                           sd = fac_noise * sqrt(2 / (p * (p + 1))))
    }

    X = array(NA, dim = c(n, p, p))
    X_nless = array(NA, dim = c(n, p, p))
    for (t in 1:n) {
      V = vector_to_symmetric(Z[t,], p)
      X[t,,] = Exp_BWS_core(V, mu)
      
      V_nless = vector_to_symmetric(Z_nless[t,], p)
      X_nless[t,,] = Exp_BWS_core(V_nless, mu)
    }
    
    return (list("X" = X, "X_nless" = X_nless, "Z" = Z, "Z_nless" = Z_nless, 
                 "A" = A, "Factors" = Factors, "mu" = mu))
  }
  
}


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


