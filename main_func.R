source("./sphere_util.R")

main = function (x, r, test_size = 0, oracle_mu = NULL) {
  # x: data that has been taken square root transformations
  # test_size: number of data that will be reserved as test set
  # oracle_mu: only used in simulation
  # streamlined function for synthetic and real data analysis
  
  n = dim(x)[1]
  d = dim(x)[2]
  p = dim(x)[3]
  
  if (test_size > 0) {
    x_test = x[-c(1:(n - test_size)),,]
    x = x[1:(n - test_size),,]
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
        temp_x_test[t,((i - 1) * p + 1):(i * p)] = x_test[t,i,]
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
  
  # compare with separate factor model
  FVU_g_sprt = rep(0, r)
  FVU_e_sprt = rep(0, r)
  V_sprt = array(NA, dim = c(d * p, r))
  for (i in 1:d) {
    model_sprt = LYB_fm(log_x[,((i - 1) * p + 1):(i * p)], r, h = 6)
    V_sprt[((i - 1) * p + 1):(i * p),1:r] = model_sprt$V
    for (k in 1:r) {
      z_hat = predict_fm(model_sprt$V[,1:k], model_sprt$mean, log_x[,((i - 1) * p + 1):(i * p)])
      x_hat = Exp_sphere(z_hat, mu_hat[i,])
      
      FVU_g_sprt[k] = FVU_g_sprt[k] + sum(geod_sphere(x_hat, x[,i,])^2)
      FVU_e_sprt[k] = FVU_e_sprt[k] + sum((x_hat - x[,i,])^2)
    }
  }
  FVU_g_sprt = FVU_g_sprt / TV_g
  FVU_e_sprt = FVU_e_sprt / TV_e
  
  if (test_size > 0) {
    pred_err_g_sprt = rep(0, r)
    pred_err_e_sprt = rep(0, r)
    for (i in 1:d) {
      model_sprt = LYB_fm(log_x[,((i - 1) * p + 1):(i * p)], r, h = 6)
      for (k in 1:r) {
        z_hat = predict_fm(V_sprt[((i - 1) * p + 1):(i * p),1:k], model_sprt$mean, 
                           log_x_test[,((i - 1) * p + 1):(i * p)])
        x_hat = Exp_sphere(z_hat, mu_hat[i,])
        pred_err_g_sprt[k] = pred_err_g_sprt[k] + sum(geod_sphere(x_test[,i,], x_hat)^2)
        pred_err_e_sprt[k] = pred_err_e_sprt[k] + sum((x_test[,i,] - x_hat)^2)
      }
    }
    pred_err_g_sprt = sqrt(pred_err_g_sprt / n_test)
    pred_err_e_sprt = sqrt(pred_err_e_sprt / n_test)
  }
  
  
  # compare with factor model with mu given by oracle
  
  if (test_size > 0) {
    return (list("FVU_g" = FVU_g, "FVU_e" = FVU_e, 
                 "pe_g" = pred_err_g, "pe_e" = pred_err_e,
                 "V" = V,
                 "FVU_g_linear" = FVU_g_linear, "FVU_e_linear" = FVU_e_linear, 
                 "pe_g_linear" = pred_err_g_linear, "pe_e_linear" = pred_err_e_linear,
                 "V_linear" = V_linear,
                 "FVU_g_sprt" = FVU_g_sprt, "FVU_e_sprt" = FVU_e_sprt,
                 "pe_g_sprt" = pred_err_g_sprt, "pe_e_sprt" = pred_err_e_sprt,
                 "V_sprt" = V_sprt))
  } 
  return (list("FVU_g" = FVU_g, "FVU_e" = FVU_e, 
               "V" = V,
               "FVU_g_linear" = FVU_g_linear, "FVU_e_linear" = FVU_e_linear,
               "V_linear" = V_linear,
               "FVU_g_sprt" = FVU_g_sprt, "FVU_e_sprt" = FVU_e_sprt,
               "V_sprt" = V_sprt))

}

















