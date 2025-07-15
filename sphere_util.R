library(Matrix)

#' check if x is on the sphere (within some tolerance)
is_on_sphere = function (x, tol = 1e-6) {
  if (abs(norm(x, "2") - 1) > tol) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}

#' computes the geodesic distance between two (array) of points on the spheres
#' 
#' If x and y are both arrays, a vector of n distances corresponding to pointwise
#' geodesic distances is returned.
#' 
#' @param x an $(n \times q)$ array of points or a $q$-dimensional vector
#' @param y an $(n \times q)$ array of points or a $q$-dimensional vector
geod_sphere = function (x, y) {
  if (is.vector(x) && is.vector(y)) {
    
    if (is_on_sphere(x)) {
      x = x / norm(x, "2")   # for numerical stability
    } else {
      stop("geod_sphere: x not on sphere")
    }
    
    if (is_on_sphere(y)) {
      y = y / norm(y, "2")   # for numerical stability
    } else {
      stop("geod_sphere: y not on sphere")
    }
    
    temp = c(x %*% y)
    res = acos(temp)
  } else if (is.matrix(x) && is.vector(y)) {
    
    if (prod(apply(x, 1, is_on_sphere)) == 1) {
      x = x / apply(x, 1, is_on_sphere)
    } else {
      stop("geod_sphere: x not on sphere")
    }
    
    if (is_on_sphere(y)) {
      y = y / norm(y, "2")   # for numerical stability
    } else {
      stop("geod_sphere: y not on sphere")
    }
    
    temp = c(x %*% y)
    res = acos(temp)
  } else if (is.vector(x) && is.matrix(y)) {
    
    if (is_on_sphere(x)) {
      x = x / norm(x, "2")   # for numerical stability
    } else {
      stop("geod_sphere: x not on sphere")
    }
    
    if (prod(apply(y, 1, is_on_sphere)) == 1) {
      y = y / apply(y, 1, is_on_sphere)
    } else {
      stop("geod_sphere: y not on sphere")
    }
    
    temp = c(y %*% x)
    res = acos(temp)
  } else if (is.matrix(x) && is.matrix(y)) {
    
    if (prod(apply(x, 1, is_on_sphere)) == 1) {
      x = x / apply(x, 1, is_on_sphere)
    } else {
      stop("geod_sphere: x not on sphere")
    }
    
    if (prod(apply(y, 1, is_on_sphere)) == 1) {
      y = y / apply(y, 1, is_on_sphere)
    } else {
      stop("geod_sphere: y not on sphere")
    }
    
    temp = c(diag(x %*% t(y)))
    res = acos(temp)
  } else {
    stop("geod_sphere: x and y must be 2-dim arrays or vectors.")
  }
  
  return (res)
}

#' exponential map on the sphere
#' 
#' @param x an array of tangent vectors
#' @param mu base point 
Exp_sphere = function (x, mu) {

  if (is.matrix(x)) {
    
    if (sum(abs(x)) == 0) {
      return (matrix(mu, nrow = nrow(x), ncol = ncol(x), byrow = T))
    }
    
    x_norm = sqrt(rowSums(x^2))
    x_norm[which(x_norm == 0)] = 1
    std_x = x / x_norm
    res = outer(cos(x_norm), mu) + sin(x_norm) * std_x
    
    res = res / sqrt(rowSums(res^2)) # normalize again to avoid numerical instability
  } else {
    
    if (sum(abs(x)) == 0) {
      return (mu)
    }
    
    x_norm = sqrt(sum(x^2))
    res = cos(x_norm) * mu + sin(x_norm) * x / x_norm
    
    res = res / sqrt(sum(res^2))
  }
  
  return (res)
}

#' logarithmic map for the sphere
#' 
#' @param x an array of points on the sphere
#' @param mu a base point
Log_sphere = function (x, mu) {

  if (is.matrix(x)) {
    w = x - matrix(mu, nrow = nrow(x), ncol = ncol(x), byrow = T)
    Proj = w - outer(c(w %*% mu), mu)
    Proj = Proj / sqrt(rowSums(Proj^2))
    
    res = acos(c(x %*% mu)) * Proj
    
    if (any(rowSums(abs(w)) < 1e-7)) {
      res[which(rowSums(abs(w)) < 1e-7),] = 0
    }
    
    
  } else {
    w = x - mu
    Proj = w - mu * c(t(mu) %*% w)
    
    if (sum(abs(w)) < 1e-7) {
      res = rep(0, length(mu))
    } else {
      res = acos(c(x %*% mu)) * Proj / sqrt(sum(Proj^2))
    }
  }
  
  return (res)
}

#' parallel transport for the sphere
#' 
#' @param x starting point
#' @param y end point
#' @param v tangent vector at x
#' 
pt_sphere = function (x, y, v) {
  
  if (abs(c(x %*% v)) > 1e-6) {
    stop("pt_sphere: v is not tangent at x")
  }
  
  e1 = x
  e2 = y - c(x %*% y) * x
  e2 = e2 / norm(e2, "2")
  
  v_perp = v - c(v %*% e2) * e2
  
  a = c(v %*% e2)
  theta = acos(c(y %*% x))
  
  res = a * (cos(theta) * e2 - sin(theta) * e1) + v_perp

  return(res)
}

#' computes the Frechet mean on the sphere
#' 
#' @param x (n by q) array of data
#' 
mean_on_sphere = function (x, tau = 0.1, tol = 1e-8, max.iter = 1000, verbose = FALSE) {

  n = nrow(x)
  mu = x[sample(n, 1),]
  for (i in 1:max.iter) {
    grad = colMeans(Log_sphere(x, mu))
    mu_new = Exp_sphere(mu + tau * grad, mu)
    
    temp = c(x %*% mu_new / sqrt(rowSums(x^2)))
    if (any(temp > 1.01) || any(temp < -1.01)) {
      stop("mean_on_sphere: something must be wrong")
    } else {
      temp = pmin(pmax(temp, -1), 1)
    }
    loss = mean(acos(temp))
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
    }
    if (verbose) {
      cat("mean_on_sphere: iter", i, "; loss", round(loss, 4), "\n")
    }
    mu = mu_new
    loss_old = loss
  }
  
  return (mu)
}

LYB_fm = function (x, r, h, demean = TRUE) {
  # x: (n by d) observation matrix
  # r: number of factors
  # h: number of lags to use
  # Estimate the factor model of Lam, Yao, and Bathia
  
  if (demean) {
    mean = colMeans(x)
    x = t(t(x) - colMeans(x))
  }
  
  # Compute the auxiliary positive definite matrix
  pd = 0
  H = h + 1
  n = nrow(x)
  for (i in 1:h) {
    temp = t(x[H:n,]) %*% x[(H - i):(n - i),] / n
    pd = pd + temp %*% t(temp)
  }
  
  # Eigenanalysis
  model = eigen(pd)
  Evec = model$vectors
  V = Evec[,1:r] # (d by r)
  evals = model$values
  
  # Extract factors and residuals
  f_hat = t(t(V) %*% t(x)) # (n by r)
  e_hat = x - f_hat %*% t(V)
  
  # Estimate the number of factors
  ratios = evals[2:r] / evals[1:(r - 1)]
  r_hat = which.min(ratios)
  
  return (list("V" = V, "f_hat" = f_hat, "e_hat" = e_hat, 
               "fitted.val" = f_hat %*% t(V), "r_hat" = r_hat,
               "mean" = mean))
}

predict_fm = function (V, mu, new_x) {
  
  if (is.matrix(new_x)) {
    z_temp = t(t(new_x) - mu)
    Factor = z_temp %*% V
    z_hat = Factor %*% t(V)
    x_hat = t(t(z_hat) + mu)
  } else if (is.vector(new_x)) {
    z_temp = new_x - mu
    Factor = c(t(z_temp) %*% V)
    z_hat = V %*% Factor
    x_hat = mu + z_hat
  } else {
    stop("predict_fm: new_x must be matrix or vector")
  }
  
  return (x_hat)
}

tan_basis_sphere = function (mu) {
  n = length(mu)
  U = svd(mu, nu = n, nv = n)$u
  
  return(U[,-1])
}

#' Helper that slices the correct indices
q_index = function (qs, j, type = "ambient") {
  if (type == "ambient") {
    NULL
  } else if (type == "intrinsic") {
    qs = qs - 1
  } else {
    stop("q_index: unsupported slicing type")
  }
  
  if (j == 1) {
    res = 1:qs[j]
  } else {
    res = sum(qs[1:(j - 1)]) + c(1:qs[j])
  }
  
  return (res)
}

#' estimate the RFM for the (product of spheres)
#' 
#' @param x a list of (n by q_j) arrays, j = 1,2,...,d
#' 
#' @export
rfm_sphere = function (x, r, h = 6, tau = 0.5, max.iter = 100) {
  d = length(x)
  n = dim(x[[1]])[1]
  qs = rep(NA, d)
  for (j in 1:d) {
    qs[j] = dim(x[[j]])[2]
  }
  
  # estimate the Frechet mean
  mus = vector("list", length = d)
  for (j in 1:d) {
    mus[[j]] = mean_on_sphere(x[[j]], tau = tau, max.iter = max.iter,
                              verbose = FALSE) 
  }
  
  # construct log-mapped data
  log_x_vec = NULL
  bases = vector("list", d)
  for (j in 1:d) {
    bases[[j]] = tan_basis_sphere(mus[[j]])
    
    log_x = Log_sphere(x[[j]], mus[[j]])
    temp = log_x %*% bases[[j]]
    
    log_x_vec = cbind(log_x_vec, temp)
  }
  if (ncol(log_x_vec) != (sum(qs) - d)) {
    stop("rfm_sphere: incorrect log_x_vec dimension") # for prototyping
  }
  
  # estimate factor model
  model = LYB_fm(log_x_vec, r = r, h = h)
  
  return (list("A" = model$V, "f_hat" = model$f_hat, "E" = bases,
               "mu_hat" = mus, "factor_model" = model,
               "r_hat" = model$r_hat))  

}

#' This function is used internally to evaluate the factor models
#' 
#' Computes evaluation metric for RFM on (product-)sphere
#' 
#' @param x_test raw test data
#' @param factor_model output from rfm_sphere
#' @param evaluation_type to compute geodesic or Euclidean distance
#' @param  fraction whether to return fraction of variance unexplained or squared prediction errors
#' 
Frac_Var_sphere = function (x_test, RFM_model, Euclidean_mean,
                            evaluation_type = "Sphere", fraction = TRUE) {
  
  mu_hat = RFM_model$mu_hat
  E = RFM_model$E  # list of q_j by (q_j-1)
  factor_model = RFM_model$factor_model
  
  V = factor_model$V
  z_mean = factor_model$mean
  r = dim(V)[2]
  
  n = nrow(x_test[[1]])
  d = length(x_test)
  qs = rep(NA, d)
  for (j in 1:d) {
    qs[j] = ncol(x_test[[j]])
  }
  
  # construct log-mapped data
  log_x_vec = NULL
  for (j in 1:d) {
    log_x = Log_sphere(x_test[[j]], mu_hat[[j]])
    temp = log_x %*% E[[j]]
    
    log_x_vec = cbind(log_x_vec, temp)
  }
  if (ncol(log_x_vec) != (sum(qs) - d)) {
    stop(paste("Frac_Var_sphere: incorrect log_x_vec dimension:", ncol(log_x_vec)))
  }
  
  # make predictions and evaluate
  res = rep(0, r)
  total_var = 0
  for (i in 1:r) {
    z_hat = predict_fm(V[,1:i], z_mean, log_x_vec)
    
    # predict
    x_hat = vector("list", d)
    for (j in 1:d) {
      idx = q_index(qs, j, "intrinsic")
      
      temp = z_hat[,idx] # n by (q_j - 1)
      x_hat[[j]] = Exp_sphere(temp %*% t(E[[j]]), mu_hat[[j]])
    }
    
    # evaluate
    if (evaluation_type == "Sphere") {
      for (j in 1:d) {
        res[i] = res[i] + sum(geod_sphere(x_hat[[j]], x_test[[j]])^2)
        if (i == 1) {
          total_var = total_var + sum(geod_sphere(mu_hat[[j]], x_test[[j]])^2)
        }
      }
    } else if (evaluation_type == "Euclidean") {
      for (j in 1:d) {
        res[i] = res[i] + sum(norm(x_hat[[j]] - x_test[[j]], "F")^2)
        if (i == 1) {
          for (m in 1:n) {
            total_var = total_var + sum((Euclidean_mean[[j]] - x_test[[j]][m,])^2)
          }
        }
      }
    } else {
      stop("Frac_Var_sphere: unsupported evaluation type")
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
#' @param mu_hat (a list of) sphere Frechet mean, externally given
#' @param evaluation_type to compute geodesic or Euclidean distances
#' @param fraction whether to return fraction of variance unexplained or squared prediction errors
#' 
Frac_Var_LYB = function (x_test, factor_model, mu_hat, Euclidean_mean,
                         evaluation_type = "Sphere",
                         fraction = TRUE) {
  
  V = factor_model$V     # (q_1 + q_2 + ... + q_d) by r
  z_mean = factor_model$mean
  r = dim(V)[2]
  
  n = nrow(x_test[[1]])
  d = length(x_test)
  qs = rep(NA, d)
  for (j in 1:d) {
    qs[j] = ncol(x_test[[j]])
  }
  
  # Naively encode the product sphere
  z = NULL
  for (j in 1:d) {
    z = cbind(z, x_test[[j]]) # n by (q_1 + ... + q_d)
  }
  
  # make predictions and evaluate
  res = rep(0, r)
  total_var = 0
  for (i in 1:r) {
    z_hat = predict_fm(V[,1:i], z_mean, z)
    
    # predict
    x_hat = vector("list", d)
    for (j in 1:d) {
      idx = q_index(qs, j, "ambient")
      temp = z_hat[,idx] # n by q_j
      
      if (evaluation_type == "Sphere") {
        x_hat[[j]] = temp / apply(temp, 1, norm, "2") # project to the sphere
      } else if (evaluation_type == "Euclidean") {
        x_hat[[j]] = temp
      } else {
        stop("Frac_Var_LYB: unsupported evaluation type")
      }
    }
    
    # evaluate
    if (evaluation_type == "Sphere") {
      for (j in 1:d) {
        res[i] = res[i] + sum(geod_sphere(x_hat[[j]], x_test[[j]])^2)
        if (i == 1) {
          total_var = total_var + sum(geod_sphere(mu_hat[[j]], x_test[[j]])^2)
        }
      }
    } else if (evaluation_type == "Euclidean") {
      for (j in 1:d) {
        res[i] = res[i] + sum(norm(x_hat[[j]] - x_test[[j]], "F")^2)
        if (i == 1) {
          for (m in 1:n) {
            total_var = total_var + sum((Euclidean_mean[[j]] - x_test[[j]][m,])^2)
          }
        }
      }
    }
  }
  
  if (fraction) {
    res = res / total_var
  }
  
  return (res)
}
















