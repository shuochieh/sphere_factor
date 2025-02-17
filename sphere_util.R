# Essential functions for spherical factor analysis
library(Matrix)

mean_on_sphere = function (x, tau = 0.1, tol = 1e-8, max.iter = 1000, verbose = FALSE) {
  # x: (n by d) matrix of data
  # tau: stepsize
  # estimate intrinsic mean on the sphere
  
  n = nrow(x)
  mu = x[sample(n, 1),]
  for (i in 1:max.iter) {
    grad = colMeans(Log_sphere(x, mu))
    mu_new = Exp_sphere(mu + tau * grad, mu)
    
    loss = mean(acos(x %*% mu_new / sqrt(rowSums(x^2))))
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

Exp_sphere = function (x, mu) {
  # x: an (m by d) matrix or a vector on the tangent space
  # mu: an d-dimensional vector of the reference point
  # tangent space --> sphere
  
  if (is.matrix(x)) {
    
    if (sum(abs(x)) == 0) {
      return (matrix(mu, nrow = nrow(x), ncol = ncol(x), byrow = T))
    }
    
    x_norm = sqrt(rowSums(x^2))
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

Log_sphere = function (x, mu) {
  # x: an (m by d) matrix or a vector on the sphere
  # mu: an d-dimensional vector of the reference point
  # sphere --> tangent space
  
  if (is.matrix(x)) {
    w = x - matrix(mu, nrow = nrow(x), ncol = ncol(x), byrow = T)
    Proj = w - outer(c(w %*% mu), mu)
    Proj = Proj / sqrt(rowSums(Proj^2))
    
    res = acos(c(x %*% mu)) * Proj
    if (any(rowSums(abs(w)) == 0)) {
      res[which(rowSums(abs(w)) == 0),] = 0
    }
    
  } else {
    w = x - mu
    Proj = w - mu * c(t(mu) %*% w)
    
    if (sum(abs(w)) == 0) {
      res = rep(0, length(mu))
    } else {
      res = acos(c(x %*% mu)) * Proj / sqrt(sum(Proj^2))
    }
  }
  
  return (res)
}

geod_sphere = function (x, y) {
  # x, y: d-dimensional unit vectors
  # x, y can also be n by d matrices (each row is a unit vector)
  # compute the geodesic distance between x, y
  if (is.vector(x) && is.vector(y)) {
    temp = t(x) %*% y
    if (any(temp > 1.001)) {
      stop("geod_sphere: clear violation of unit vectors")
    } else if (any(temp < -1.001)) {
      stop("geod_sphere: clear violation of unit vectors")
    }
    temp = max(min(temp, 1), -1)
    res = acos(temp)
  } else if (is.matrix(x) && is.vector(y)) {
    temp = c(x %*% y)
    if (any(temp > 1.001)) {
      stop("geod_sphere: clear violation of unit vectors")
    } else if (any(temp < -1.001)) {
      stop("geod_sphere: clear violation of unit vectors")
    }
    temp = pmax(pmin(temp, 1), -1)
    res = acos(temp)
  } else if (is.vector(x) && is.matrix(y)) {
    temp = c(y %*% x)
    if (any(temp > 1.001)) {
      stop("geod_sphere: clear violation of unit vectors")
    } else if (any(temp < -1.001)) {
      stop("geod_sphere: clear violation of unit vectors")
    }
    temp = pmax(pmin(temp, 1), -1)
    res = acos(temp)
  } else if (is.matrix(x) && is.matrix(y)) {
    temp = c(diag(x %*% t(y)))
    if (any(temp > 1.001)) {
      stop("geod_sphere: clear violation of unit vectors")
    } else if (any(temp < -1.001)) {
      stop("geod_sphere: clear violation of unit vectors")
    }
    temp = pmax(pmin(temp, 1), -1)
    res = acos(temp)
  } else {
    stop("geod_sphere: x, y has to be vectors or matrices")
  }
  
  return (res)
}

perp_proj = function (x, mu) {
  # x: d-dimensional vector
  # mu: d-dimensional vector
  # Projects x to the perpendicular subspace of mu
  
  res = mu * c(t(mu) %*% x) / (norm(mu, type = "2")^2)
  res = x - res
  
  return (res)
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
  Evec = eigen(pd)$vectors
  V = Evec[,1:r] # (d by r)
  
  # Extract factors and residuals
  f_hat = t(t(V) %*% t(x)) # (n by r)
  e_hat = x - f_hat %*% t(V)
  
  return (list("V" = V, "f_hat" = f_hat, "e_hat" = e_hat, 
               "fitted.val" = f_hat %*% t(V),
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