library(vMF)

source("./sphere_util.R")

perp_proj = function (x, mu) {
  # x: d-dimensional vector
  # mu: d-dimensional vector
  # Projects x to the perpendicular subspace of mu
  
  res = mu * c(t(mu) %*% x) / (norm(mu, type = "2")^2)
  res = x - res
  
  return (res)
}

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
  # spec: model specification
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
  lambda_tilde = sqrt((1 - 0.8^2) / ((p - 1) * r * EA2))
  for (t in 1:n) {
    Z[t,] = lambda_tilde * A_tilde %*% c(Factors[t,])
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
  
  return (list("x" = X, "A_tilde" = A_tilde, "X_nless" = X_noiseless, "
               Factors" = Factors,
               "Z" = Z, "mus" = mus))
}

main = function (x, r, test_size = 0) {
  # x: data that has been taken square root transformations
  # test_size: proportion of data that will be reserved as test set
  # streamlined function for synthetic and real data analysis
  
  n = dim(x)[1]
  d = dim(x)[2]
  p = dim(x)[3]
  
  if (test_size > 0) {
    x_test = x[-c(1:floor(n * (1 - test_size))),,]
    x = x[1:floor(n * (1 - test_size)),,]
    n = dim(x)[1]
    n_test = dim(x_test)[1]
    
    log_x_test = array(NA, dim = c(n_test, d * p))
  }
  
  # estimate mu
  mu_hat = array(NA, dim = c(d, p))
  log_x = array(NA, dim = c(n, d * p))
  geod_to_mean = array(NA, dim = c(n, d))
  for (i in 1:d) {
    mu_hat[i,] = mean_on_sphere(x[,i,])
    log_x[,((i - 1) * p + 1):(i * p)] = Log_sphere(x[,i,], mu_hat[i,])
    if (test_size > 0) {
      log_x_test[,((i - 1) * p + 1):(i * p)] = Log_sphere(x_test[,i,], mu_hat[i,])
    }
    
    geod_to_mean[,i] = geod_sphere(mu_hat[i,], x[,i,])
  }
  
  # estimate factor model
  model = LYB_fm(log_x, r, h = 6)
  V = model$V
  Factors = model$f_hat
  
  # fractions of variance explained
  FVU_g = rep(0, r)  # geodesic distance
  FVU_e = rep(0, r)  # Euclidean distance
  TV_g = 0
  TV_e = 0
  for (k in 1:r) {
    z_hat = predict_fm(V[,1:k], model$mean, log_x)
    for (i in 1:d) {
      x_hat = Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,])
      FVU_g[k] = FVU_g[k] + sum(geod_sphere(x[,i,], x_hat)^2)
      FVU_e[k] = FVU_e[k] + sum((x[,i,] - x_hat)^2)
      
      if (k == 1) {
        TV_g = TV_g + sum(geod_sphere(x[,i,], mu_hat[i,])^2)
        TV_e = TV_e + sum((t(x[,i,]) - mu_hat[i,])^2)
      }
    }
  }
  FVU_g = FVU_g / TV_g
  FVU_e = FVU_e / TV_e
  
  if (test_size > 0) {
    pred_err_g = rep(0, r)
    pred_err_e = rep(0, r)
    for (k in 1:r) {
      z_hat = predict_fm(V[,1:k], model$mean, log_x_test)
      for (i in 1:d) {
        x_hat = Exp_sphere(z_hat[,((i - 1) * p + 1):(i * p)], mu_hat[i,])
        pred_err_g[k] = pred_err_g[k] + sum(geod_sphere(x_test[,i,], x_hat)^2)
        pred_err_e[k] = pred_err_e[k] + sum((x_test[,i,] - x_hat)^2)
      }
    }
    pred_err_g = sqrt(pred_err_g / n_test)
    pred_err_e = sqrt(pred_err_e / n_test)
  }
  
  # compare with linear factor model
  FVU_g_linear = rep(0, r)
  FVU_e_linear = rep(0, r)
  
  temp_x = array(NA, dim = c(n, d * p))
  for (t in 1:n) {
    for (i in 1:d) {
      temp_x[t,((i - 1) * p + 1):(i * p)] = x[t,i,]
    }
  }
  
  model_linear = LYB_fm(temp_x, r, h = 6)
  V_linear = model_linear$V
  for (k in 1:r) {
    x_hat = predict_fm(V_linear[,1:k], model_linear$mean, temp_x)
    FVU_e_linear[k] = norm(x_hat - temp_x, "F")^2 / TV_e
    for (i in 1:d) {
      proj_xhat = x_hat[,((i - 1) * p + 1):(i * p)] / 
                  apply(x_hat[,((i - 1) * p + 1):(i * p)], 1, norm, "2")
      FVU_g_linear[k] = FVU_g_linear[k] + sum(geod_sphere(x[,i,], proj_xhat)^2)
    }
    FVU_g_linear[k] = FVU_g_linear[k] / TV_g
  }
  
  if (test_size > 0) {
    pred_err_g_linear = rep(0, r)
    pred_err_e_linear = rep(0, r)
    temp_x_test = array(NA, dim = c(n_test, d * p))
    for (t in 1:n_test) {
      for (i in 1:d) {
        temp_x_test[t,i,] = x_test[t,j,]
      }
    }
    
    for (k in 1:r) {
      x_hat = predict_fm(V_linear[,1:k], model_linear$mean, temp_x_test)
      pred_err_e_linear[k] = sqrt(norm(x_hat - temp_x_test, "F")^2 / n_test)
      
      for (i in 1:d) {
        proj_xhat = x_hat[,((i - 1) * p + 1):(i * p)] /
                    apply(x_hat[,((i - 1) * p + 1):(i * p)], 1, norm, "2")
        pred_err_g_linear[k] = pred_err_g_linear[k] + sum(geod_sphere(x_test[,i,],
                                                                      proj_xhat)^2)
      }
      pred_err_g_linear[k] = sqrt(pred_err_g_linear[k] / n_test)
    }
  }
  
  if (test_size > 0) {
    return (list("FVU_g" = FVU_g, "FVU_e" = FVU_e, "pe_g" = pred_err_g, "pe_e" = pred_err_e,
                 "FVU_g_linear" = FVU_g_linear, "FVU_e_linear" = FVU_e_linear, 
                 "pe_g_linear" = pred_err_g_linear,
                 "pe_e_linear" = pred_err_e_linear))
  } 
  return (list("FVU_g" = FVU_g, "FVU_e" = FVU_e, 
               "FVU_g_linear" = FVU_g_linear, "FVU_e_linear" = FVU_e_linear))

}

















