library(vMF)
library(ggplot2)
source("./main_func.R")

noise_inject = function (x, type = "uniform", b = NULL, kappa = NULL) {
  # x: a unit vector
  # b: parameter of uniform noise 
  # type: {uniform, vMF} type of noise distribution
  #       If vMF is used, intensity parameter kappa should be carefully supplied
  # Inject noise around a vector x on the unit sphere
  
  d = length(x)
  if (type == "uniform") {
    if (is.null(b)) {
      stop("noise_inject: b should be supplied when uniform distribution is used")
    }
    res = x + runif(d, -b, b)
    res = Exp_sphere(perp_proj(res, x), x)
  } else if (type == "vMF") {
    if (is.null(kappa)) {
      stop("noise_inject: kappa should be supplied when vMF distribution is used")
    }
    
    res = rvMF(1, kappa * x)
  }
  
  return (res)
}

dta_gen = function (n, d, p, r, lambda = 0.1, spec = 0,
                    noise_type = "uniform", b = NULL, kappa = NULL) {
  # n: sample size
  # d: number of spherical components
  # p: dimension of each sphere = p - 1 (encoded as R^p vectors)
  # r: number of true factors
  # lambda: scaling of the latent variable (controls signal strength, E|z_t|^2)
  # spec: model specification {0, 1, 2}
  # noise_type: type of noise distribution ("uniform" or "vMF"). vMF =  von Mises-Fisher distribution
  # b, kappa: parameters for noise distributions
  ### Generate synthetic data
  
  # generate reference points
  mus = array(rnorm(d * p), dim = c(d, p))
  mus = mus / rowSums(abs(mus))
  mus = sign(mus) * sqrt(abs(mus))
  
  # generate latent factor process
  Factors = array(0, dim = c(n + 100, r))
  for (t in 2:(n + 100)) {
    Factors[t,] = 0.8 * Factors[t - 1,] + rnorm(r)
  }
  Factors = Factors[-c(1:100),] # n by r
  
  # generate factor loadings
  if (spec == 0) {
    # default setting
    A = array(runif(d * (p - 1) * r, min = -1, max = 1), 
              dim = c(d * (p - 1), r))
    EA2 = 1 / 3
  } else if (spec == 1) {
    A = array(runif(d * (p - 1) * r, min = -2, max = 2), 
              dim = c(d * (p - 1), r))
    EA2 = 4 / 3
  } else if (spec == 2) {
    A = array(runif(d * (p - 1) * r, min = -0.5, max = 0.5), 
              dim = c(d * (p - 1), r))
    EA2 = 1 / 12
  } else {
    stop(paste("dta_gen: invalid specification:", spec))
  }
  # A is (p - 1 + p - 1 + ... + p - 1) by r
  
  # generate latent observations (embeded in the ambient space)
  Vs = vector("list", length = d)
  for (j in 1:d) {
    Vs[[j]] = svd(mus[j,], nu = p, nv = p)$u[,-1] # Coordinate system on the tangent space encoded in R^p
  }
  V = bdiag(Vs) # (p + ... + p) by (p - 1 + ... + p - 1)
  
  Z = array(NA, dim = c(n, d * p))
  A_tilde = as.matrix(V %*% A) # (p + ... + p) by r
  lambda_tilde = sqrt((1 - 0.8^2) * lambda / ((p - 1) * r * EA2))
  for (t in 1:n) {
    Z[t,] = lambda_tilde * (A_tilde %*% c(Factors[t,]))
  }
  
  # generate observations
  X = array(NA, dim = c(n, d, p))
  X_noiseless = array(NA, dim = c(n, d, p))
  for (j in 1:d) {
    X_noiseless[,j,] = Exp_sphere(Z[,((j - 1) * p + 1):(j * p)], mus[j,])
    for (t in 1:n) {
      X[t,j,] = noise_inject(X_noiseless[t,j,], type = noise_type,
                             b = b, kappa = kappa)
    }
  }
  
  return (list("x" = X, "A_tilde" = A_tilde, "x_nless" = X_noiseless, "
               Factors" = Factors,
               "Z" = Z, "mus" = mus))
}

simulation_core = function (n, d, p, r_true, r_model, lambda, spec,
                            noise_type, b = NULL, kappa = NULL) {
  N = n + 300
  dta = dta_gen(n = N, d = d, p = p, r = r_true, lambda = lambda, spec = spec,
                noise_type = noise_type, b = b, kappa = kappa)
  
  x = dta$x
  A_tilde = dta$A_tilde
  x_nless = dta$X_nless

  res = main(x, r = r_model, test_size = 300)
  
  FVU_e = res$FVU_e
  FVU_e_linear = res$FVU_e_linear
  FVU_e_sprt = res$FVU_e_sprt
  pe_e = res$pe_e
  pe_e_linear = res$pe_e_linear
  pe_e_sprt = res$pe_e_sprt
  pe_g = res$pe_g
  pe_g_linear = res$pe_g_linear
  pe_g_sprt = res$pe_g_sprt
  

  return (list("FVU" = c(FVU_e, FVU_e_linear, FVU_e_sprt),
               "pe" = c(pe_e, pe_e_linear, pe_e_sprt),
               "pe_g" = c(pe_g, pe_g_linear, pe_g_sprt)))
}

simulation = function (n_sim, n, d, p, r_true, r_model, lambda, spec,
                       noise_type, b = NULL, kappa = NULL) {
  for (zz in 1:n_sim) {
    if (zz == 1) {
      temp = simulation_core(n, d, p, r_true, r_model, lambda, spec,
                             noise_type, b, kappa)
      FVU = array(NA, c(n_sim, length(temp$FVU)))
      pe = array(NA, c(n_sim, length(temp$pe)))
      pe_g = array(NA, c(n_sim, length(temp$pe)))
      
      FVU[1,] = temp$FVU
      pe[1,] = temp$pe
      pe_g[1,] = temp$pe_g
    } else {
      temp = simulation_core(n, d, p, r_true, r_model, lambda, spec,
                             noise_type, b, kappa)
      FVU[zz,] = temp$FVU
      pe[zz,] = temp$pe
      pe_g[zz,] = temp$pe_g
    }
    
    cat("iteration", zz, "\n")
  }
  
  return (list("FVU" = FVU, "pe" = pe, "pe_g" = pe_g))
}

n_sim = 300
n = 300
d = 5
p = 10
r_true = 3
r_model = 6 
lambda = 0.01
spec = 0
noise_type = "uniform"
b = 0.00 * sqrt(1/p)
kappa = 5 * p

res = simulation(n_sim, n, d, p, r_true, r_model, lambda, spec,
                 noise_type, b, kappa)

latex_print = function (v) {
  m = length(v)
  for (i in 1:m) {
    if (i == m) {
      cat(" ", v[i], "\n")
    } else {
      cat(" ", v[i], " &")
    }
  }
}

# sphere (FVU)
latex_print(round(colMeans(res$FVU[,1:r_model]), 3))
latex_print(round(apply(res$FVU[,1:r_model], 2, sd), 3))
apply(round(apply(res$FVU[,1:r_model], 2, quantile, probs = c(0.05, 0.95)), 3), 1, latex_print)

# linear (FVU)
latex_print(round(colMeans(res$FVU[,(r_model + 1):(2 * r_model)]), 3))
latex_print(round(apply(res$FVU[,(r_model + 1):(2 * r_model)], 2, sd), 3))
apply(round(apply(res$FVU[,(r_model + 1):(2 * r_model)], 2, quantile, probs = c(0.05, 0.95)), 3), 1, latex_print)

# separate (FVU)
latex_print(round(colMeans(res$FVU[,(2 * r_model + 1):(3 * r_model)]), 3))
latex_print(round(apply(res$FVU[,(2 * r_model + 1):(3 * r_model)], 2, sd), 3))
apply(round(apply(res$FVU[,(2 * r_model + 1):(3 * r_model)], 2, quantile, probs = c(0.05, 0.95)), 3), 1, latex_print)

# sphere 
latex_print(round(colMeans(res$pe[,1:r_model]), 3))
latex_print(round(apply(res$pe[,1:r_model], 2, sd), 3))
apply(round(apply(res$pe[,1:r_model], 2, quantile, probs = c(0.05, 0.95)), 3), 1, latex_print)

# linear
latex_print(round(colMeans(res$pe[,(r_model + 1):(2 * r_model)]), 3))
latex_print(round(apply(res$pe[,(r_model + 1):(2 * r_model)], 2, sd), 3))
apply(round(apply(res$pe[,(r_model + 1):(2 * r_model)], 2, quantile, probs = c(0.05, 0.95)), 3), 1, latex_print)

# separate
latex_print(round(colMeans(res$pe[,(2 * r_model + 1):(3 * r_model)]), 3))
latex_print(round(apply(res$pe[,(2 * r_model + 1):(3 * r_model)], 2, sd), 3))
apply(round(apply(res$pe[,(2 * r_model + 1):(3 * r_model)], 2, quantile, probs = c(0.05, 0.95)), 3), 1, latex_print)









