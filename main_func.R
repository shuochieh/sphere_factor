# load whichever is needed (but not both)
# source("./sphere_util.R")
# source("./BWS_util.R")

subspace_d = function (U, V, type = "sine-theta") {
  if (type == "trace-projection") {
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
    
    res = sqrt(res)
  }
  if (type == "sine-theta") {
    if (is.vector(U)) {
      if (!is.vector(V)) {
        # stop("subspace_d: U and V must have the same dimension for sine-theta distance")
        return (Inf)
      }
      res = acos(c(t(U) %*% V))
      res = sin(res)
    } else if (is.vector(V)) {
      if (!is.vector(U)) {
        # stop("subspace_d: U and V must have the same dimension for sine-theta distance")
        return (Inf)
      }
      res = acos(c(t(U) %*% V))
      res = sin(res)
    } else {
      if (ncol(U) != ncol(V)) {
        # stop("subspace_d: U and V must have the same dimension for sine-theta distance")
        return (Inf)
      }
      Q1 = qr.Q(qr(U))
      Q2 = qr.Q(qr(V))
      
      res = pmin(pmax(svd(t(Q1) %*% Q2)$d, 0), 1) # for numerical stability
      
      res = max(sin(acos(res)))
    }
  }

  return (res)
}

#' Streamlined function for synthetic and real data analysis with BWS-data
#' 
#' @param x SPD data 
#' @param r number of factors used
#' @param test_size the number of data reserved as test set
#' @param h number of lags used in estimating factor model
#' @param true_A true factor loading used in simulation
#' @param true_mu true Frechet mean used in simulation
#' @param fraction returns fraction of variance unexplained; returns squared prediction errors if false
#' 
main_BWS = function (x, r, test_size = 0, h = 6, batch_size = NULL, max.iter = 100,
                     mu_hat = NULL, true_A = NULL, true_mu = NULL, fraction = TRUE,
                     return_predictions = FALSE) {
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  Euclidean_mean = colMeans(x[1:(n - test_size),,], dims = 1)
  
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
  model = rfm_bws(x, r = r, h = h, batch_size = batch_size, max.iter = max.iter,
                  mu_hat = mu_hat)
  V = model$A
  Factors = model$f_hat
  mu_hat = model$mu_hat
  E = model$E
  E_lyapunov = model$E_lyapunov
  r_hat_RFM = model$r_hat
  z_bar = model$factor_model$mean
  
  # compare estimated loading space
  if (!is.null(true_A)) {
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
    res1 = Frac_Var_bws(x_test, model, "BWS", fraction, 
                        return_predictions = return_predictions,
                        Euclidean_mean = Euclidean_mean)
    res2 = Frac_Var_bws(x_test, model, "Euclidean", fraction, 
                        return_predictions = return_predictions,
                        Euclidean_mean = Euclidean_mean)
  } else {
    res1 = Frac_Var_bws(x, model, "BWS", fraction, 
                        return_predictions = return_predictions,
                        Euclidean_mean = Euclidean_mean)
    res2 = Frac_Var_bws(x, model, "Euclidean", fraction, 
                        return_predictions = return_predictions,
                        Euclidean_mean = Euclidean_mean)
  }
  
  # estimate linear factor model
  # encode data as vectors
  x_vector = array(NA, dim = c(n, p * (p + 1) / 2))
  for (m in 1:n) {
    x_vector[m,] = symmetric_to_vector(x[m,,])
  }
  model = LYB_fm(x_vector, r = r, h = h)
  r_hat_LYB = model$r_hat
  V_LYB = model$V
  
  if (test_size > 0) {
    res3 = Frac_Var_LYB(x_test, model, mu_hat, "BWS", fraction, 
                        return_predictions = return_predictions,
                        Euclidean_mean = Euclidean_mean)
    res4 = Frac_Var_LYB(x_test, model, mu_hat, "Euclidean", fraction, 
                        return_predictions = return_predictions,
                        Euclidean_mean = Euclidean_mean)
  } else {
    res3 = Frac_Var_LYB(x, model, mu_hat, "BWS", fraction, 
                        return_predictions = return_predictions,
                        Euclidean_mean = Euclidean_mean)
    res4 = Frac_Var_LYB(x, model, mu_hat, "Euclidean", fraction, 
                        return_predictions = return_predictions,
                        Euclidean_mean = Euclidean_mean)
  }
  
  if (return_predictions) {
    return (list("mu_hat" = mu_hat, 
                 "FVU_RFM_BWS" = res1$res, "FVU_RFM_Euc" = res2$res,
                 "FVU_LYB_BWS" = res3$res, "FVU_LYB_Euc" = res4$res,
                 "loading_dist" = subspace_dist,
                 "r_hat_RFM" = r_hat_RFM, "r_hat_LYB" = r_hat_LYB,
                 "Factors" = Factors, "z_bar" = z_bar,
                 "V" = V, "V_LYB" = V_LYB,
                 "E" = E, "E_lyapunov" = E_lyapunov,
                 "RFM_xhat" = res1$xhat, 
                 "LYB_xhat" = res3$xhat))
    
  }
  
  return (list("mu_hat" = mu_hat, 
               "FVU_RFM_BWS" = res1, "FVU_RFM_Euc" = res2,
               "FVU_LYB_BWS" = res3, "FVU_LYB_Euc" = res4,
               "loading_dist" = subspace_dist,
               "r_hat_RFM" = r_hat_RFM, "r_hat_LYB" = r_hat_LYB,
               "Factors" = Factors, "z_bar" = z_bar,
               "V" = V,  "V_LYB" = V_LYB,
               "E" = E, "E_lyapunov" = E_lyapunov))
}

#' Streamlined function for synthetic and real data analysis with (product-)sphere-valued
#' data
#' 
#' @param x a list of (n by q_j) array of data (q_j: ambient dimension)
#' @param r number of factors used
#' @param test_size the number of data reserved as test set
#' 
main_sphere = function (x, r, test_size = 0, h = 6, tau = 0.5, max.iter = 100,
                        true_A = NULL, true_mu = NULL, true_E = NULL, fraction = TRUE) {

  n = nrow(x[[1]])
  d = length(x)
  qs = rep(NA, d)
  for (i in 1:d) {
    qs[i] = ncol(x[[i]])
  }
  
  if (test_size > 0) {
    x_test = vector("list", d)
    for (i in 1:d) {
      x_test[[i]] = x[[i]][-c(1:(n - test_size)),, drop = FALSE]
      x[[i]] = x[[i]][1:(n - test_size),, drop = FALSE]
    }
    n = nrow(x[[1]])
    n_test = test_size
  }
  
  Euclidean_mean = vector("list", d)
  for (j in 1:d) {
    Euclidean_mean[[j]] = colMeans(x[[j]])
  }
  
  # estimate RFM
  model = rfm_sphere(x, r = r, h = h, tau = tau, max.iter = max.iter)
  V = model$A
  Factors = model$f_hat
  mu_hat = model$mu_hat
  E = model$E
  r_hat_RFM = model$r_hat
  
  # compare estimated loading space
  if (!is.null(true_A)) {
    transported_V = array(NA, dim = c(dim(true_A)[1], r))
    subspace_dist = rep(NA, r)
    
    for (i in 1:r) {
      for (j in 1:d) {
        idx = q_index(qs, j, "intrinsic")
        temp = pt_sphere(mu_hat[[j]], true_mu[[j]], c(E[[j]] %*% V[idx,i])) # q_j by 1
        transported_V[idx, i] = c(temp %*% true_E[[j]])
      }
      subspace_dist[i] = subspace_d(true_A, transported_V[,1:i])
    }
  } else {
    subspace_dist = NULL
  }
  
  # results
  if (test_size > 0) {
    res1 = Frac_Var_sphere(x_test, model, "Sphere", fraction, Euclidean_mean = Euclidean_mean)
    res2 = Frac_Var_sphere(x_test, model, "Euclidean", fraction, Euclidean_mean = Euclidean_mean)
  } else {
    res1 = Frac_Var_sphere(x, model, "Sphere", fraction, Euclidean_mean = Euclidean_mean)
    res2 = Frac_Var_sphere(x, model, "Euclidean", fraction, Euclidean_mean = Euclidean_mean)
  }
  
  # estimate linear model
  x_vector = NULL
  for (j in 1:d) {
    x_vector = cbind(x_vector, x[[j]]) # n by (q_1+...+q_d)
  }
  model = LYB_fm(x_vector, r = r, h = h)
  r_hat_LYB = model$r_hat
  
  if (test_size > 0) {
    res3 = Frac_Var_LYB(x_test, model, mu_hat, "Sphere", fraction, Euclidean_mean = Euclidean_mean)
    res4 = Frac_Var_LYB(x_test, model, mu_hat, "Euclidean", fraction, Euclidean_mean = Euclidean_mean)
  } else {
    res3 = Frac_Var_LYB(x, model, mu_hat, "Sphere", fraction, Euclidean_mean = Euclidean_mean)
    res4 = Frac_Var_LYB(x, model, mu_hat, "Euclidean", fraction, Euclidean_mean = Euclidean_mean)
  }
  
  return (list("mu_hat" = mu_hat,
               "FVU_RFM_Sphere" = res1, "FVU_RFM_Euc" = res2,
               "FVU_LYB_Sphere" = res3, "FVU_LYB_Euc" = res4, 
               "loading_dist" = subspace_dist,
               "r_hat_RFM" = r_hat_RFM, "r_hat_LYB" = r_hat_LYB))
}















