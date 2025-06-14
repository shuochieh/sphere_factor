source("./sphere_util.R")
source("./BWS_util.R")

#' This function is used internally to evaluate the factor models
#' 
#' Computes evaluation metric for RFM
#' 
#' @param x_test raw test data
#' @param factor_model output from rfm_bws; an LYB_fm object
#' @param evaluation_type to compute BWS distance or Euclidean (Frobenius distance)
#' @param  fraction whether to return fraction of variance unexplained or squared prediction errors
#' 
Frac_Var_bws = function (x_test, RFM_model, evaluation_type = "BWS", fraction = TRUE) {
  
  mu_hat = RFM_model$mu_hat
  E = RFM_model$E
  E_lyapunov = RFM_model$E_lyapunov
  factor_model = RFM_model$factor_model
  
  V = factor_model$V
  z_mean = factor_model$mean
  r = dim(V)[2]
  
  if (length(dim(x_test)) == 3) {
    n = dim(x_test)[1]
    p = dim(x_test)[2]
    x.is.array = TRUE
  } else if (length(dim(x_test)) == 2) {
    n = 1
    p = dim(x_test)[1]
    x.is.array = FALSE
  }
  
  # construct log-mapped data
  log_x_vec = log_vec_construct(x_test, mu_hat, E_lyapunov)
  
  # make predictions and evaluate
  res = rep(0, r)
  total_var = 0
  for (i in 1:r) {
    z_hat = predict_fm(V[,1:i], z_mean, log_x_vec)
    
    # predict
    if (x.is.array) {
      x_hat = array(NA, dim = c(n, p, p))
      for (m in 1:n) {
        temp = log_to_tangent(z_hat[m,], E)
        x_hat[m,,] = Exp_BWS(temp, mu_hat)
      }
    } else {
      temp = log_to_tangent(z_hat, E)
      x_hat = Exp_BWS(temp, mu_hat)
    }
    
    # evaluate
    if (evaluation_type == "BWS") {
      if (x.is.array) {
        for (m in 1:n) {
          res[i] = res[i] + geod_BWS_core(x_hat[m,,], x_test[m,,])^2
          if (i == 1) {
            total_var = total_var + geod_BWS_core(mu_hat, x_test[m,,])^2
          }
        }
      } else {
        res[i] = geod_BWS_core(x_hat, x_test)^2
      }
    } else if (evaluation_type == "Euclidean") {
      if (x.is.array) {
        Euclidean_mu = colMeans(x_test, dims = 1)
        for (m in 1:n) {
          res[i] = res[i] + norm(x_hat[m,,] - x_test[m,,], type = "F")^2
          if (i == 1) {
            total_var = total_var + norm(Euclidean_mu - x_test[m,,], type = "F")^2
          }
        }
      } else {
        res = norm(x_hat - x_test, type = "F")^2
      }
    } else {
      stop("Frac_Var_bws: unsupported evaluation type")
    }
  }
  
  if (fraction) {
    res = res / total_var
  }
  
  return (res)
}

#' This function is used internally to evaluate the factor models
#' 
#' Computes evaluation metric for linear factor model
#' 
#' @param x_test raw test data
#' @param factor_model output from LYB_fm
#' @param mu_hat BWS Frechet mean, externally given
#' @param evaluation_type to compute BWS distance or Euclidean (Frobenius distance)
#' @param fraction whether to return fraction of variance unexplained or squared prediction errors
#' @param epsilon used in projecting predictions to SPD (only when evaluation type == "BWS")
#' 
Frac_Var_LYB = function (x_test, factor_model, mu_hat,
                         evaluation_type = "BWS",
                         fraction = TRUE,
                         epsilon = 1e-6) {
  V = factor_model$V
  z_mean = factor_model$mean
  r = dim(V)[2]
  
  if (length(dim(x_test)) == 3) {
    n = dim(x_test)[1]
    p = dim(x_test)[2]
    x.is.array = TRUE
  } else if (length(dim(x_test)) == 2) {
    n = 1
    p = dim(x_test)[1]
    x.is.array = FALSE
  }
  
  # encode SPD as vectors
  if (x.is.array) {
    z = array(NA, dim = c(n, p * (p + 1) / 2))
    for (m in 1:n) {
      z[m,] = symmetric_to_vector(x_test[m,,])
    }
  } else {
    z = symmetric_to_vector(x_test)
  }
  
  # make predictions and evaluate
  res = rep(0, r)
  total_var = 0
  for (i in 1:r) {
    z_hat = predict_fm(V[,1:i], z_mean, z)
    
    # predict
    if (x.is.array) {
      x_hat = array(NA, dim = c(n, p, p))
      for (m in 1:n) {
        x_hat[m,,] = vector_to_symmetric(z_hat[m,], p)
      }
    } else {
      x_hat = vector_to_symmetric(z_hat, p)
    }
    
    # evaluate
    if (evaluation_type == "BWS") {
      if (x.is.array) {
        for (m in 1:n) {
          temp = project_to_SPD(x_hat[m,,], epsilon)
          res[i] = res[i] + (Re(geod_BWS_core(temp, x_test[m,,])))^2
          if (i == 1) {
            total_var = total_var + (Re(geod_BWS_core(mu_hat, x_test[m,,])))^2
          }
        }
      } else {
        res[i] = (Re(geod_BWS_core(project_to_SPD(x_hat, epsilon), x_test)))^2
      }
    } else if (evaluation_type == "Euclidean") {
      if (x.is.array) {
        Euclidean_mu = colMeans(x_test, dims = 1)
        for (m in 1:n) {
          res[i] = res[i] + norm(x_hat[m,,] - x_test[m,,], type = "F")^2
          if (i == 1) {
            total_var = total_var + norm(Euclidean_mu - x_test[m,,], type = "F")^2
          }
        }
      } else {
        res = norm(x_hat - x_test, type = "F")^2
      }
    } else {
      stop("Frac_Var_bws: unsupported evaluation type")
    }
  }
  
  if (fraction) {
    res = res / total_var
  }
  
  return (res)
}

subspace_d = function (U, V) {
  if (is.vector(U)) {
    U = U / sqrt(sum(U^2))
    P1 = U %*% t(U)
    q1 = 1
  } else {
    Q = qr.Q(qr(U))
    P1 = Q %*% t(Q)
    q1 = ncol(U)
  }
  if (is.vector(V)) {
    V = V / sqrt(sum(V^2))
    P2 = V %*% t(V)
    q2 = 1
  } else {
    Q = qr.Q(qr(V))
    P2 = Q %*% t(Q)
    q2 = ncol(V)
  }
  
  res = 1 - sum(diag(P1 %*% P2)) / max(q1, q2)
  
  return (sqrt(res))
}

#' Streamlined function for synthetic and real data analysis with BWS-data
#' 
#' @param dta SPD data generated by dta_gen
#' @param r number of factors used
#' @param test_size the number of data reserved as test set
#' @param h number of lags used in estimating factor model
#' @param true_A true factor loading used in simulation
#' @param true_mu true Frechet mean used in simulation
#' @param fraction returns fraction of variance unexplained; returns squared prediction errors if false
#' 
main_BWS = function (x, r, test_size = 0, h = 6, batch_size = NULL, max.iter = 100,
                     true_A = NULL, true_mu = NULL, fraction = TRUE) {
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  if (!is.null(true_A)) {
    A = true_A
  }

  if (test_size > 0) {
    x_test = x[-c(1:(n - test_size)),,]
    x = x[1:(n - test_size),,]
    n = dim(x)[1]
    n_test = dim(x_test)[1]
  }
  
  # estimate RFM
  model = rfm_bws(x, r = r, h = h, batch_size = batch_size, max.iter = max.iter)
  V = model$A
  Factors = model$f_hat
  mu_hat = model$mu_hat
  E = model$E
  E_lyapunov = model$E_lyapunov
  r_hat_RFM = model$r_hat
  
  # compare estimated loading space
  if (!is.null(A)) {
    transported_V = array(NA, dim = c(p * (p + 1) / 2, r))
    subspace_dist = rep(NA, r)
    for (i in 1:r) {
      v_to_matrix = log_to_tangent(V[,i], E)
      transport_matrix = pt_bws(mu_hat, true_mu, v_to_matrix)
      transported_V[,i] = tangent_in_E(transport_matrix, dta$mu)
      
      subspace_dist[i] = subspace_d(A, transported_V[,1:i])
    }
  } else {
    subspace_dist = NULL
  }

  # results
  if (test_size > 0) {
    res1 = Frac_Var_bws(x_test, model, "BWS", fraction)
    res2 = Frac_Var_bws(x_test, model, "Euclidean", fraction)
  } else {
    res1 = Frac_Var_bws(x, model, "BWS", fraction)
    res2 = Frac_Var_bws(x, model, "Euclidean", fraction)
  }
  
  # estimate linear factor model
  # encode data as vectors
  x_vector = array(NA, dim = c(n, p * (p + 1) / 2))
  for (m in 1:n) {
    x_vector[m,] = symmetric_to_vector(x[m,,])
  }
  model = LYB_fm(x_vector, r = r, h = h)
  r_hat_LYB = model$r_hat
  
  if (test_size > 0) {
    res3 = Frac_Var_LYB(x_test, model, mu_hat, "BWS", fraction)
    res4 = Frac_Var_LYB(x_test, model, mu_hat, "Euclidean", fraction)
  } else {
    res3 = Frac_Var_LYB(x, model, mu_hat, "BWS", fraction)
    res4 = Frac_Var_LYB(x, model, mu_hat, "Euclidean", fraction)
  }
  
  return (list("mu_hat" = mu_hat, 
               "FVU_RFM_BWS" = res1, "FVU_RFM_Euc" = res2,
               "FVU_LYB_BWS" = res3, "FVU_LYB_Euc" = res4,
               "loading_dist" = subspace_dist,
               "r_hat_RFM" = r_hat_RFM, "r_hat_LYB" = r_hat_LYB))
}

main_uneven_sphere = function (x, r, test_size = 0, h = 6) {
  # x: data that has been taken square root transformations
  #    which is a list of n by p_i matrices
  #    with length of the list = d
  # test_size: number of data that will be reserved as test set
  # streamlined function for synthetic and real data analysis
  
  n = nrow(x[[1]])
  d = length(x)
  ps = rep(NA, d)
  for (i in 1:d) {
    ps[i] = ncol(x[[i]])
  }
  
  if (test_size > 0) {
    x_test = vector("list", d)
    for (i in 1:d) {
      x_test[[i]] = x[[i]][-c(1:(n - test_size)),, drop = FALSE]
      x[[i]] = x[[i]][1:(n - test_size),, drop = FALSE]
    }
    n = nrow(x[[1]])
    n_test = nrow(x_test[[1]])
    
    log_x_test = matrix(nrow = n_test, ncol = sum(ps))
  }
  
  # estimate mu
  mu_hat = vector("list", d)
  mu_hat_Euclidean = vector("list", d)
  log_x = matrix(NA, nrow = n, ncol = sum(ps))
  geod_to_mean = array(NA, dim = c(n, d))
  for (i in 1:d) {
    mu_hat[[i]] = mean_on_sphere(x[[i]])
    mu_hat_Euclidean[[i]] = colMeans(x[[i]])
    if (i == 1) {
      log_x[,1:ps[1]] = Log_sphere(x[[i]], mu_hat[[i]])
      if (test_size > 0) {
        log_x_test[,1:ps[1]] = Log_sphere(x_test[[i]], mu_hat[[i]])
      }
    } else {
      log_x[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))] = Log_sphere(x[[i]], mu_hat[[i]])
      if (test_size > 0) {
        log_x_test[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))] = Log_sphere(x_test[[i]], mu_hat[[i]])
      }
    }

    geod_to_mean[,i] = geod_sphere(mu_hat[[i]], x[[i]])
  }
  
  # estimate factor model
  model = LYB_fm(log_x, r, h = h)
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
      if (i == 1) {
        x_hat = Exp_sphere(z_hat[,1:ps[1]], mu_hat[[i]])
      } else {
        x_hat = Exp_sphere(z_hat[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))], mu_hat[[i]])
      }
      FVU_g[k] = FVU_g[k] + sum(geod_sphere(x[[i]], x_hat)^2)
      FVU_e[k] = FVU_e[k] + sum((x[[i]] - x_hat)^2)

      if (k == 1) {
        TV_g = TV_g + sum(geod_sphere(x[[i]], mu_hat[[i]])^2)
        TV_e = TV_e + sum((t(x[[i]]) - mu_hat_Euclidean[[i]])^2)
      }
    }
  }
  FVU_g = FVU_g / TV_g
  FVU_e = FVU_e / TV_e
  
  if (test_size > 0) {
    pred_err_g = matrix(0, ncol = r, nrow = d)
    pred_err_e = matrix(0, ncol = r, nrow = d)
    for (k in 1:r) {
      z_hat = predict_fm(V[,1:k], model$mean, log_x_test)
      for (i in 1:d) {
        if (i == 1) {
          x_hat = Exp_sphere(z_hat[,1:ps[1]], mu_hat[[i]])
        } else {
          x_hat = Exp_sphere(z_hat[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))], mu_hat[[i]])
        }
        pred_err_g[i,k] = pred_err_g[i,k] + sum(geod_sphere(x_test[[i]], x_hat)^2)
        
        # Euclidean prediction error for the raw data
        pred_err_e[i,k] = pred_err_e[i,k] + sum((x_test[[i]]^2 - x_hat^2)^2)
      }
    }
    pred_err_g = sqrt(pred_err_g / n_test)
    pred_err_e = sqrt(pred_err_e / n_test)
  }
  
  # compare with linear factor model
  FVU_g_linear = rep(0, r)
  FVU_e_linear = rep(0, r)
  
  temp_x = NULL
  for (i in 1:d) {
    temp_x = cbind(temp_x, x[[i]])
  }

  model_linear = LYB_fm(temp_x, r, h = h)
  V_linear = model_linear$V
  for (k in 1:r) {
    x_hat = predict_fm(V_linear[,1:k], model_linear$mean, temp_x)
    FVU_e_linear[k] = norm(x_hat - temp_x, "F")^2 / TV_e
    for (i in 1:d) {
      if (i == 1) {
        proj_xhat = x_hat[,1:ps[1]] / 
          apply(x_hat[,1:ps[1]], 1, norm, "2")
      } else {
        proj_xhat = x_hat[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))] / 
          apply(x_hat[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))], 1, norm, "2")
      }
      FVU_g_linear[k] = FVU_g_linear[k] + sum(geod_sphere(x[[i]], proj_xhat)^2)
    }
    FVU_g_linear[k] = FVU_g_linear[k] / TV_g
  }
  
  if (test_size > 0) {
    pred_err_g_linear = matrix(0, ncol = r, nrow = d)
    pred_err_e_linear = matrix(0, ncol = r, nrow = d)
    temp_x_test = array(NA, dim = c(n_test, sum(ps)))
    for (i in 1:d) {
      if (i == 1) {
        temp_x_test[,1:sum(ps[1])] = x_test[[1]]
      } else {
        temp_x_test[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))] = x_test[[i]]
      }
    }

    for (k in 1:r) {
      x_hat = predict_fm(V_linear[,1:k], model_linear$mean, temp_x_test)
      # pred_err_e_linear[k] = sqrt(norm(x_hat^2 - temp_x_test^2, "F")^2 / n_test)
      
      for (i in 1:d) {
        if (i == 1) {
          proj_xhat = x_hat[,1:ps[1]] / 
            apply(x_hat[,1:ps[1]], 1, norm, "2")
          pred_err_e_linear[1,k] = sqrt(norm(x_hat[,1:ps[1]]^2 - temp_x_test[,1:ps[1]]^2, "F")^2 / n_test)
        } else {
          proj_xhat = x_hat[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))] / 
            apply(x_hat[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))], 1, norm, "2")
          pred_err_e_linear[i,k] = sqrt(norm(x_hat[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))]^2 - 
                                               temp_x_test[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))]^2, "F")^2 / n_test)
        }
        pred_err_g_linear[i,k] = pred_err_g_linear[i,k] + sum(geod_sphere(x_test[[i]], proj_xhat)^2)
      }
      
      pred_err_g_linear[,k] = sqrt(pred_err_g_linear[,k] / n_test)
    }
  }
  
  # Compare with linear factor model applied directly to raw data
  temp_x = NULL
  for (i in 1:d) {
    temp_x = cbind(temp_x, x[[i]]^2)
  }
  
  model_linear_direct = LYB_fm(temp_x, r, h = h)
  V_linear = model_linear_direct$V
  
  if (test_size > 0) {
    pred_err_e_linear_direct = matrix(0, ncol = r, nrow = d)
    temp_x_test = array(NA, dim = c(n_test, sum(ps)))
    for (i in 1:d) {
      if (i == 1) {
        temp_x_test[,1:sum(ps[1])] = x_test[[1]]^2
      } else {
        temp_x_test[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))] = x_test[[i]]^2
      }
    }
    
    for (k in 1:r) {
      x_hat = predict_fm(V_linear[,1:k], model_linear_direct$mean, temp_x_test)
      for (i in 1:d) {
        if (i == 1) {
          pred_err_e_linear_direct[1,k] = sqrt(norm(x_hat[,1:ps[1]] - temp_x_test[,1:ps[1]], "F")^2 / n_test)
        } else {
          pred_err_e_linear_direct[i,k] = sqrt(norm(x_hat[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))] - 
                                                      temp_x_test[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))], "F")^2 / n_test)
        }
      } 
    }
  }
  
  # compare with separate factor model
  r = min(min(ps) - 1, r)
  FVU_g_sprt = rep(0, r)
  FVU_e_sprt = rep(0, r)
  V_sprt = array(NA, dim = c(sum(ps), r))
  
  for (i in 1:d) {
    if (i == 1) {
      model_sprt = LYB_fm(log_x[,1:sum(ps[1])], r, h = h)
      V_sprt[1:ps[1], 1:r] = model_sprt$V
    } else {
      model_sprt = LYB_fm(log_x[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))], r, h = 6)
      V_sprt[(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i])),1:r] = model_sprt$V
    }
    for (k in 1:r) {
      if (i == 1) {
        z_hat = predict_fm(model_sprt$V[,1:k], model_sprt$mean, log_x[,1:ps[1]])
      } else {
        z_hat = predict_fm(model_sprt$V[,1:k], model_sprt$mean, log_x[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))])
      }
      x_hat = Exp_sphere(z_hat, mu_hat[[i]])

      FVU_g_sprt[k] = FVU_g_sprt[k] + sum(geod_sphere(x_hat, x[[i]])^2)
      FVU_e_sprt[k] = FVU_e_sprt[k] + sum((x_hat - x[[i]])^2)
    }
  }
  FVU_g_sprt = FVU_g_sprt / TV_g
  FVU_e_sprt = FVU_e_sprt / TV_e

  if (test_size > 0) {
    pred_err_g_sprt = matrix(0, ncol = r, nrow = d)
    pred_err_e_sprt = matrix(0, ncol = r, nrow = d)
    for (i in 1:d) {
      if (i == 1) {
        model_sprt = LYB_fm(log_x[,1:ps[1]], r, h = 6)
      } else {
        model_sprt = LYB_fm(log_x[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))], r, h = 6)
      }
      for (k in 1:r) {
        if (i == 1) {
          z_hat = predict_fm(V_sprt[1:ps[1],1:k], model_sprt$mean,
                             log_x_test[,1:ps[1]])
        } else {
          z_hat = predict_fm(V_sprt[(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i])),1:k], model_sprt$mean,
                             log_x_test[,(sum(ps[1:(i - 1)]) + 1):(sum(ps[1:i]))])
        }
        x_hat = Exp_sphere(z_hat, mu_hat[[i]])
        pred_err_g_sprt[i,k] = pred_err_g_sprt[i,k] + sum(geod_sphere(x_test[[i]], x_hat)^2)
        
        pred_err_e_sprt[i,k] = pred_err_e_sprt[i,k] + sum((x_test[[i]]^2 - x_hat^2)^2)
      }
    }
    pred_err_g_sprt = sqrt(pred_err_g_sprt / n_test)
    pred_err_e_sprt = sqrt(pred_err_e_sprt / n_test)
  }
  
  
  if (test_size > 0) {
    return (list("geod_to_mean" = geod_to_mean,
                 "FVU_g" = FVU_g, "FVU_e" = FVU_e, 
                 "pe_g" = pred_err_g, "pe_e" = pred_err_e,
                 "V" = V,
                 "FVU_g_linear" = FVU_g_linear, "FVU_e_linear" = FVU_e_linear, 
                 "pe_g_linear" = pred_err_g_linear, "pe_e_linear" = pred_err_e_linear,
                 "V_linear" = V_linear,
                 "FVU_g_sprt" = FVU_g_sprt, "FVU_e_sprt" = FVU_e_sprt,
                 "pe_g_sprt" = pred_err_g_sprt, "pe_e_sprt" = pred_err_e_sprt,
                 "pe_e_linear_direct" = pred_err_e_linear_direct,
                 "TV_e" = TV_e))
  } 
  return (list("geod_to_mean" = geod_to_mean,
               "FVU_g" = FVU_g, "FVU_e" = FVU_e, 
               "V" = V,
               "FVU_g_linear" = FVU_g_linear, "FVU_e_linear" = FVU_e_linear, 
               "V_linear" = V_linear,
               "FVU_g_sprt" = FVU_g_sprt, "FVU_e_sprt" = FVU_e_sprt,
               "TV_e" = TV_e))
}















