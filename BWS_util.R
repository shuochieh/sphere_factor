# Essential functions for Bures-Wasserstein factor analysis
library(maotai)
library(expm)

geod_BWS_core = function (X, Y) {
    temp = sqrtm(X)
    temp = sqrtm(temp %*% Y %*% temp)
    res = sum(diag(X)) + sum(diag(Y)) - 2 * sum(diag(temp))
    
    return (sqrt(res))
}

geod_BWS = function (X, Y) {
  
  dimX = length(dim(X))
  dimY = length(dim(Y))
  
  if (dimX == 3 && dimY == 3) {
    if (dim(X)[1] != dim(Y)[1]) {
      stop("geod_BWS: X and Y have different sample sizes")
    }
    n = dim(X)[1]
    res = rep(NA, n)
    for (i in 1:n) {
      res[i] = geod_BWS_core(X[i,,], Y[i,,])
    }
  } else if (dimX == 3 && dimY == 2) {
    n = dim(X)[1]
    res = rep(NA, n)
    for (i in 1:n) {
      res[i] = geod_BWS_core(X[i,,], Y)
    }
  } else if (dimX == 2 && dimY == 3) {
    n = dim(Y)[1]
    res = rep(NA, n)
    for (i in 1:n) {
      res[i] = geod_BWS_core(X, Y[i,,])
    }
  } else if (dimX == 2 && dimY == 2) {
    res = geod_BWS_core(X, Y)
  } else {
    stop(paste("geod_BWS: incorrect dimensions X:",
               dim(X),
               "Y:", dim(Y)))
  }
  
  return (res)
}

Exp_BWS_core = function (X, M) {
  L = lyapunov(M ,X)
  p = dim(X)[1]
  diag(L) = diag(L) + 1
  output = L %*% M %*% L
  
  return (output)
}

Exp_BWS = function (X, M) {
  # X: (n by p by p) or (p by p) array of data on the tangent space
  # M: (p by p) reference point
  # tangent space --> sphere
  
  output = array(NA, dim = dim(X))
  
  if (length(dim(X)) == 3) {
    for (i in 1:dim(X)[1]) {
      output[i,,] = Exp_BWS_core(X[i,,], M)
    }
  } else if (length(dim(X)) == 2) {
    output = Exp_BWS_core(X, M)
  } else {
    stop(paste("Exp_BWS: Incorrect dimensions X", dim(X)))
  }
  
  return (output)
}

Log_BWS_core = function (X, M) {
  
  # check if X in the injectivity radius
  test_X = lyapunov(M, X)
  diag(test_X) = diag(test_X) + 1
  eig = eigen(test_X, symmetric = T, only.values = T)$values
  if (!all(eig > 0)) {
    stop("Log_BWS: X outside injectivity radius")
  }
  
  output = sqrtm(X %*% M)
  output = output + sqrtm(M %*% X) - 2 * M
  
  return (output)
}

Log_BWS = function (X, M) {
  if (length(dim(X)) == 3) {
    output = array(NA, dim = dim(X))
    for (i in 1:dim(X)[1]) {
      output[i,,] = Log_BWS_core(X[i,,], M)
    }
  } else if (length(dim(X)) == 2) {
    output = Log_BWS_core(X, M)
  } else {
    stop(paste("Log_BWS: Incorrect dimensions X", dim(X)))
  }
  
  return (output)
}

symmetric_to_vector = function (X) {
  return (X[upper.tri(X, diag = TRUE)])
}

vector_to_symmetric = function (v, n) {
  res = matrix(NA, n, n)
  res[upper.tri(res, diag = TRUE)] = v
  res[lower.tri(res)] = t(res)[lower.tri(res)]
  
  return (res)
}

is.spd = function (x, tol = 1e-8) {
  if (!isSymmetric.matrix(x)) {
    return (FALSE)
  } 
  lambda = eigen(x, symmetric = T, only.values = T)$values
  if (all(lambda > tol)) {
    return (TRUE)
  } else {
    return (FALSE)
  }
}

project_to_SPD <- function(A, epsilon = 0) {
  # epsilon is to ensure positive-definiteness (ocassionally used)
  
  eig <- eigen(A)
  eig$values[eig$values < epsilon] <- epsilon
  A_SPD <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
  return(A_SPD)
}

mean_on_BWS = function (X, tau = 0.1, tol = 1e-4, max.iter = 1000,
                             verbose = FALSE) {
  # X: (n by ...) data
  # tau: step size
  
  n = dim(X)[1]
  mu = X[1,,]
  for (i in 1:max.iter) {
    grad = apply(Log_BWS(X, mu), c(2, 3), mean)
    mu_new = Exp_BWS(tau * grad, mu)
    
    loss = mean(geod_BWS(X, mu_new))
    if (verbose) {
      cat("mean_on_BWS: iter", i, "loss", round(loss, 4), "\n")
    }
    if (i > 1 && (loss_old - loss < tol)) {
      mu = mu_new
      break
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

#' Computes a canonical basis for the tangent space at a point in BWS
#' 
#' @param Sigma Base point
#' @return An array of orthonormal basis
#' @export
tan_basis_bws = function (Sigma) {
  p = dim(Sigma)[1]
  
  model = eigen(Sigma)
  lambdas = model$values
  P = model$vectors
  
  # construct canonical basis
  E = array(NA, dim = c(p * (p + 1) / 2, p, p))
  E_lyapunov = array(NA, dim = c(p * (p + 1) / 2, p, p))
  counter = 0
  for (i in 1:p) {
    for (j in i:p) {
      counter = counter + 1
      S = matrix(0, ncol = p, nrow = p)
      S_tilde = matrix(0, ncol = p, nrow = p)
      if (i == j) {
        S[i,j] = sqrt(2 * (lambdas[i] + lambdas[j]))
        S_tilde[i,j] = 1 / sqrt(lambdas[i])
      } else {
        S[i,j] = sqrt(lambdas[i] + lambdas[j])
        S[j,i] = sqrt(lambdas[i] + lambdas[j])
        
        S_tilde[i,j] = 1 / sqrt(lambdas[i] + lambdas[j])
        S_tilde[j,i] = 1 / sqrt(lambdas[i] + lambdas[j])
      }
      E[counter,,] = P %*% S %*% t(P)
      E_lyapunov[counter,,] = P %*% S_tilde %*% t(P) # L_{mu_hat}(E)
    }
  }
  
  return (list("E" = E, "E_lyapunov" = E_lyapunov))
}


#' Estimate the RFM in Bures-Wasserstein manifold
#' 
#' This function estimates the Riemannian factor model (RFM) with data being
#' symmetric positive definite matrices equipped with the Bures-Wasserstein
#' metric.
#' 
#' @param x an $n \times p \times p$ array of data where sample is a $p \times p$ SPD matrix
#' @param r number of factors to extract
#' @param h number of lags used in estimating the factor model
#' @param mu_tol optimization tolerance level for estimating the Fr\'{e}chet mean
#' @return 
#' \describe{
#'  \item{A}{Estimated loading matrix}
#'  \item{f_hat}{Estimated factor process}
#'  \item{E}{A set of orthonormal basis for the tangent space}
#'  \item{mu_hat}{Estimated Fr\'{e}chet mean}
#' } 
#' @export
rfm_bws = function (x, r, h = 6, mu_tol) {
  n = dim(x)[1]
  p = dim(x)[2]
  
  # Estimate mu
  mu_hat = mean_on_BWS(x, tol = mu_tol, verbose = FALSE)
  
  # Construct a set of orthonormal basis
  coord = tan_basis_bws(mu_hat)
  E = coord$E
  E_lyapunov = coord$E_lyapunov
  

  # construct log-mapped data
  log_x = Log_BWS(x, mu_hat)
  log_x_vec = array(NA, dim = c(n, p * (p + 1) / 2))
  counter = 0
  for (m in 1:n) {
    for (i in 1:p) {
      for (j in i:p) {
        counter = counter + 1
        log_x_vec[m, counter] = 0.5 * sum(diag(E_lyapunov[counter,,] %*% log_x))
      }
    }
  }
  
  # Estimate the factor model
  model = LYB_fm(log_x_vec, r = r, h = h)
  
  return(list("A" = model$V, "f_hat" = model$f_hat, "E" = E, "mu_hat" = mu_hat))
}












