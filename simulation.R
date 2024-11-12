### Simulations for factor model for spherical time series

library(Matrix)

Exp_sphere = function (x, mu) {
  # x: an (m by d) matrix or a vector on the tangent space
  # mu: an d-dimensional vector of the reference point
  # tangent space --> sphere
  
  if (is.matrix(x)) {
    x_norm = sqrt(rowSums(x^2))
    std_x = x / x_norm
    res = outer(cos(x_norm), mu) + sin(x_norm) * std_x
    
    res = res / sqrt(rowSums(res^2)) # normalize again to avoid numerical instability
  } else {
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

LYB_fm = function (x, r, h, demean = TRUE) {
  # x: (n by d) observation matrix
  # r: number of factors
  # h: number of lags to use
  # Estimate the factor model of Lam, Yao, and Bathia
  
  if (demean) {
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
               "fitted.val" = f_hat %*% t(V)))
}

perp_proj = function (x, mu) {
  # x: d-dimensional vector
  # mu: d-dimensional vector
  # Projects x to the perpendicular subspace of mu
  
  res = mu * c(t(mu) %*% x) / (norm(mu, type = "2")^2)
  res = x - res
  
  return (res)
}

dta_gen = function (n, ds, r, b1 = 0.1, b2 = 0.1, spec = 1) {
  # n: sample size
  # ds: a vector of the dimensions of the ambient Euclidean spaces for each spherical
  #     component
  # r: number of (true) factors
  # Generate synthetic data
  
  if (spec == 1) {
    m = length(ds)
    
    # generate reference points
    mus = vector("list", length = m)
    for (i in 1:m) {
      temp = rnorm(ds[i])
      mus[[i]] = temp / norm(temp, type = "2")
    }
    
    # generate latent factor processes
    Fs = matrix(0, nrow = r, ncol = n + 100)
    for (t in 2:(n + 100)) {
      Fs[,t] = 0.9 * Fs[,t - 1] + runif(r, min = -b1, max = b1)
    }
    Fs = Fs[,-c(1:100)] # r by n
    
    # generate factor loadings
    A = matrix(runif((sum(ds) - m) * r, min = -2, max = 2), nrow = sum(ds) - m, ncol = r)
    # A is (d1 + d2 + ... + dm - m) by r
    
    # generate latent observations (embeded in the ambient space)
    Z = matrix(NA, nrow = n, ncol = sum(ds))
    Vs = vector("list", length = m)
    for (j in 1:m) {
      Vs[[j]] = svd(mus[[j]], nu = ds[j], nv = ds[j])$u[,-1] # Coordinate system
    }
    V = bdiag(Vs) # (d1 + ... + dm) by (d1 + ... + dm - m)
    
    A_tilde = as.matrix(V %*% A) # (d1 + ... + dm) by r
    for (t in 1:n) {
      Z[t,] = A_tilde %*% c(Fs[,t])
    }
    
    # generate observations
    X = matrix(NA, nrow = n, ncol = sum(ds))
    for (j in 1:m) {
      if (j == 1) {
        lower_idx = 0
      } else {
        lower_idx = sum(ds[1:(j - 1)])
      }
      idx_range = (lower_idx + 1):(lower_idx + ds[j])
      temp = Exp_sphere(t(mus[[j]] + t(Z[,idx_range])), mus[[j]]) # n by dj
      for (t in 1:n) {
        X[t,idx_range] = Exp_sphere(perp_proj(temp[t,] + runif(ds[j], min = -b2, max = b2), 
                                              temp[t,]), 
                                    temp[t,])
      }
    }
  }
  
  if (spec == 2) {
    m = length(ds)
    
    # generate reference points
    mus = vector("list", length = m)
    for (i in 1:m) {
      temp = rnorm(ds[i])
      mus[[i]] = temp / norm(temp, type = "2")
    }
    
    # generate latent factor processes
    Fs = matrix(0, nrow = r, ncol = n + 100)
    for (t in 2:(n + 100)) {
      Fs[,t] = 0.9 * Fs[,t - 1] + runif(r, min = -b1, max = b1)
    }
    Fs = Fs[,-c(1:100)] # r by n
    
    # generate factor loadings
    A = matrix(runif((sum(ds) - m) * r, min = -2, max = 2), nrow = sum(ds) - m, ncol = r)
    A[sample(sum(ds) - m, floor(0.3 * (sum(ds) - m))),] = 0
    # A is (d1 + d2 + ... + dm - m) by r
    
    # generate latent observations (embeded in the ambient space)
    Z = matrix(NA, nrow = n, ncol = sum(ds))
    Vs = vector("list", length = m)
    for (j in 1:m) {
      Vs[[j]] = svd(mus[[j]], nu = ds[j], nv = ds[j])$u[,-1] # Coordinate system
    }
    V = bdiag(Vs) # (d1 + ... + dm) by (d1 + ... + dm - m)
    
    A_tilde = as.matrix(V %*% A) # (d1 + ... + dm) by r
    for (t in 1:n) {
      Z[t,] = A_tilde %*% c(Fs[,t])
    }
    
    # generate observations
    X = matrix(NA, nrow = n, ncol = sum(ds))
    for (j in 1:m) {
      if (j == 1) {
        lower_idx = 0
      } else {
        lower_idx = sum(ds[1:(j - 1)])
      }
      idx_range = (lower_idx + 1):(lower_idx + ds[j])
      temp = Exp_sphere(t(mus[[j]] + t(Z[,idx_range])), mus[[j]]) # n by dj
      for (t in 1:n) {
        X[t,idx_range] = Exp_sphere(perp_proj(temp[t,] + runif(ds[j], min = -b2, max = b2), 
                                              temp[t,]), 
                                    temp[t,])
      }
    }
  }
  
  if (spec == 3) {
    m = length(ds)
    
    # generate reference points
    mus = vector("list", length = m)
    for (i in 1:m) {
      temp = rnorm(ds[i])
      mus[[i]] = temp / norm(temp, type = "2")
    }
    
    # generate latent factor processes
    Fs = matrix(0, nrow = r, ncol = n + 100)
    for (t in 2:(n + 100)) {
      Fs[,t] = 0.9 * Fs[,t - 1] + runif(r, min = -b1, max = b1)
    }
    Fs = Fs[,-c(1:100)] # r by n
    
    # generate factor loadings
    A = matrix(runif((sum(ds) - m) * r, min = -0.5, max = 0.5), nrow = sum(ds) - m, ncol = r)
    # A is (d1 + d2 + ... + dm - m) by r
    
    # generate latent observations (embeded in the ambient space)
    Z = matrix(NA, nrow = n, ncol = sum(ds))
    Vs = vector("list", length = m)
    for (j in 1:m) {
      Vs[[j]] = svd(mus[[j]], nu = ds[j], nv = ds[j])$u[,-1] # Coordinate system
    }
    V = bdiag(Vs) # (d1 + ... + dm) by (d1 + ... + dm - m)
    
    A_tilde = as.matrix(V %*% A) # (d1 + ... + dm) by r
    for (t in 1:n) {
      Z[t,] = A_tilde %*% c(Fs[,t])
    }
    
    # generate observations
    X = matrix(NA, nrow = n, ncol = sum(ds))
    for (j in 1:m) {
      if (j == 1) {
        lower_idx = 0
      } else {
        lower_idx = sum(ds[1:(j - 1)])
      }
      idx_range = (lower_idx + 1):(lower_idx + ds[j])
      temp = Exp_sphere(t(mus[[j]] + t(Z[,idx_range])), mus[[j]]) # n by dj
      for (t in 1:n) {
        X[t,idx_range] = Exp_sphere(perp_proj(temp[t,] + runif(ds[j], min = -b2, max = b2), 
                                              temp[t,]), 
                                    temp[t,])
      }
    }
  }
  
  
  return (list("X" = X, "A_tilde" = A_tilde, "A" = A, "Fs" = Fs, "Z" = Z, "mus" = mus, "Vs" = Vs))
}

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

subspace_loss = function (A_hat , A) {
  
  Q1 = qr.Q(qr(A))
  Q2 = qr.Q(qr(A_hat))
  
  if (ncol(Q1) != ncol(Q2)) {
    stop("subspace_loss: Q1 and Q2 is of different dimensions")
  }
  q = ncol(Q1)
  temp = 1 - (1 / q) * sum(diag(Q1 %*% t(Q1) %*% Q2 %*% t(Q2)))
  
  return (sqrt(temp))
}

simulation = function (n, ds, r, sd1, sd2, spec, mu_tol = 1e-3) {
  dta = dta_gen (n, ds, r, sd1, sd2, spec)
  
  m = length(ds)
  mu_hat = vector("list", length = m)
  trans_x = matrix(NA, nrow = n, ncol = sum(ds))
  trans_x_oracle = matrix(NA, nrow = n, ncol = sum(ds))
  for (j in 1:m) {
    lower_idx = ifelse(j == 1, 0, sum(ds[1:(j - 1)]))
    idx_range = (lower_idx + 1):(lower_idx + ds[j])
    mu_hat[[j]] = mean_on_sphere(dta$X[,idx_range], tol = mu_tol)
    trans_x[,idx_range] = Log_sphere(dta$X[,idx_range], mu_hat[[j]])
    trans_x_oracle[,idx_range] = Log_sphere(dta$X[,idx_range], dta$mus[[j]])
  }
  
  # Quality of mu estimate
  mu_loss = rep(NA, m)
  for (j in 1:m) {
    mu_loss[j] = norm(mu_hat[[j]] - dta$mu[[j]], "2") # acos(c(mu_hat[[j]] %*% dta$mu[[j]]))
  }
  
  # joint modeling
  model = LYB_fm(trans_x, r, floor(sum(ds)^(1/3)))
  loss1 = subspace_loss(model$V, dta$A_tilde)
  RSS1 = (norm(trans_x_oracle - model$fitted.val, "2")^2) / n
  
  # separate modeling
  RSS2 = 0
  for (j in 1:m) {
    lower_idx = ifelse(j == 1, 0, sum(ds[1:(j - 1)]))
    idx_range = (lower_idx + 1):(lower_idx + ds[j])
    model = LYB_fm(trans_x[,idx_range], r, floor(sum(ds)^(1/3)))
    RSS2 = RSS2 + (norm(trans_x_oracle[,idx_range] - model$fitted.val, "2")^2) / n
  }
  
  return (list("subspace_loss" = loss1, "RSS1" = RSS1, "RSS2" = RSS2,
               "mu_loss" = mu_loss, "dta" = dta, "mu_hat" = mu_hat))
}

# Simulation
max.sim = 500
n = 600
ds = rep(10, 5)
r = 3
snr = 0.6
b1 = 0.5 * snr * pi / (20 * r * sqrt(max(ds)))
b2 = 0.5 * (1 - snr) * pi / (20 * r * sqrt(max(ds)))
spec = 3

A_loss = rep(0, max.sim)
RSS1 = rep(0, max.sim)
RSS2 = rep(0, max.sim)
mu_loss = rep(0, max.sim)

for (sim in 1:max.sim) {
  res = simulation(n, ds, r, b1, b2, spec)
  
  A_loss[sim] = res$subspace_loss
  RSS1[sim] = sqrt(res$RSS1)
  RSS2[sim] = sqrt(res$RSS2)
  mu_loss[sim] = mean(res$mu_loss)
  
  # cat("iteration", sim, "\n")
}

mean(A_loss); sd(A_loss)
# mean(RSS1); sd(RSS1)
# mean(RSS2); sd(RSS2)
mean(RSS1 / RSS2); sd(RSS1 / RSS2)

mean(mu_loss); sd(mu_loss)
















 


