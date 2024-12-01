### Simulations for factor model for spherical time series

library(Matrix)
library(vMF)

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

dta_gen = function (n, ds, r, lambda = 0.1, b2 = 0.1, spec = 1) {
  # n: sample size
  # ds: a vector of the dimensions of the ambient Euclidean spaces for each 
  #     spherical component
  # r: number of true factors
  # spec: 1. strong factors
  #       2. weak factors (all loaded)
  # lambda: scaling of the latent variable (controls signal strength)
  # Generate synthetic data
  
  m = length(ds)
  
  # generate reference points
  mus = vector("list", length = m)
  for (i in 1:m) {
    temp = rnorm(ds[i])
    mus[[i]] = temp / norm(temp, type = "2")
  }
  
  # generate latent factor processes
  if (spec %in% c(1:3)) {
    Fs = matrix(0, nrow = r, ncol = n + 100)
    for (t in 2:(n + 100)) {
      Fs[,t] = 0.8 * Fs[,t - 1] + rnorm(r)
    }
    Fs = Fs[,-c(1:100)] # r by n
  }
  
  # generate factor loadings
  if (spec == 1) {
    A = matrix(runif((sum(ds) - m) * r, min = -2, max = 2), 
               nrow = sum(ds) - m, ncol = r)
  } else if (spec == 2) {
    A = matrix(runif((sum(ds) - m) * r, min = -0.5, max = 0.5), 
               nrow = sum(ds) - m, ncol = r)
  } else if (spec == 3) {
    A = matrix(runif((sum(ds) - m) * r, min = -0.1, max = 0.1), 
               nrow = sum(ds) - m, ncol = r)
  }
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
    Z[t,] = lambda * A_tilde %*% c(Fs[,t])
  }
  
  # generate observations
  X = matrix(NA, nrow = n, ncol = sum(ds))
  X_noiseless = matrix(NA, nrow = n, ncol = sum(ds))
  for (j in 1:m) {
    if (j == 1) {
      lower_idx = 0
    } else {
      lower_idx = sum(ds[1:(j - 1)])
    }
    idx_range = (lower_idx + 1):(lower_idx + ds[j])
    X_noiseless[,idx_range] = Exp_sphere(t(mus[[j]] + t(Z[,idx_range])), mus[[j]]) # n by dj
    temp = X_noiseless[,idx_range]
    for (t in 1:n) {
      X[t, idx_range] = rvMF(1, temp[t,])
      ####
      # Remember to adjust noise level
      ####
      
      # X[t,idx_range] = Exp_sphere(perp_proj(temp[t,] + runif(ds[j], min = -b2, max = b2), 
      #                                       temp[t,]), 
      #                             temp[t,])
    }
  }
  
  return (list("X" = X, "A_tilde" = A_tilde, "X_nless" = X_noiseless, "Fs" = Fs,
               "Z" = Z, "mus" = mus))
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
  temp = max(temp, 0) # avoid numerical issues
  
  return (sqrt(temp))
}

parallel_transport = function (a, x1, x2) {
  # a: a d-dimensional vector or a (d by r) matrix
  # x1, x2: d-dimensional unit vectors
  # parallel transport a from the tangent space of x1 to x2
  
  d = length(x1)
  theta = acos(t(x1) %*% x2)[1]
  e1 = x1
  temp = svd(cbind(x1, x2), nu = d, nv = d)$u[,-c(1,2)]
  
  H = qr(cbind(e1, temp))
  Q = qr.Q(H)
  en = x2 - Q %*% t(Q) %*% x2
  en = c(en / norm(en, "2"))
  
  R = diag(1, d) + sin(theta) * (outer(en, e1) - outer(e1, en)) + (cos(theta) - 1) * (outer(e1, e1) + outer(en, en))
  
  if (is.vector(a)) {
    res = R %*% a
    return (c(res))
  } else {
    return (R %*% a)
  }
}

oos_eval = function (x, A_hat, mu_hat, noiseless_x) {
  # x: n by p observation matrix
  # A_hat: p by r loading matrix
  # noiseless_x: n by p target to compare against
  z = Log_sphere(x, mu_hat)
  z_hat = z %*% A_hat %*% t(A_hat)
  x_hat = Exp_sphere(z_hat, mu_hat)
  
  res = acos(c(diag(noiseless_x %*% t(x_hat))))
  
  return (sum(res))
}

simulation = function (n, ds, r, b1, b2, spec, h = 5, mu_tol = 1e-3) {
  dta = dta_gen(2 * n, ds, r, b1, b2, spec)
  X_train = dta$X[1:n,]
  X_test = dta$X[-c(1:n),]
  X_nless = dta$X_nless[-c(1:n),]
  
  m = length(ds)
  mu_hat = vector("list", length = m)
  trans_x = matrix(NA, nrow = n, ncol = sum(ds))
  for (j in 1:m) {
    lower_idx = ifelse(j == 1, 0, sum(ds[1:(j - 1)]))
    idx_range = (lower_idx + 1):(lower_idx + ds[j])
    mu_hat[[j]] = mean_on_sphere(X_train[,idx_range], tol = mu_tol)
    trans_x[,idx_range] = Log_sphere(X_train[,idx_range], mu_hat[[j]])
  }
  
  # Quality of mu estimate
  mu_loss = rep(NA, m)
  for (j in 1:m) {
    mu_loss[j] = acos(c(mu_hat[[j]] %*% dta$mu[[j]]))
  }
  
  # joint modeling
  model = LYB_fm(trans_x, r, h)
  loading_est = model$V
  loading_est_transport = matrix(NA, ncol = ncol(model$V), nrow = nrow(model$V))
  dist_space_each_joint = rep(NA, m)
  for (j in 1:m) {
    lower_idx = ifelse(j == 1, 0, sum(ds[1:(j - 1)]))
    idx_range = (lower_idx + 1):(lower_idx + ds[j])
    temp = parallel_transport(loading_est[idx_range,], mu_hat[[j]], dta$mu[[j]])
    loading_est_transport[idx_range,] = temp

    dist_space_each_joint[j] = subspace_loss(temp, dta$A_tilde[idx_range,])
  }
  
  dist_space_joint = subspace_loss(loading_est_transport, dta$A_tilde)

  x_temp = matrix(NA, nrow = n, ncol = ncol(X_test))
  x_temp_oracle = matrix(NA, nrow = n, ncol = ncol(X_test))
  for (j in 1:m) {
    lower_idx = ifelse(j == 1, 0, sum(ds[1:(j - 1)]))
    idx_range = (lower_idx + 1):(lower_idx + ds[j])
    x_temp[,idx_range] = Log_sphere(X_test[,idx_range], mu_hat[[j]])
    x_temp_oracle[,idx_range] = Log_sphere(X_test[,idx_range], dta$mus[[j]])
  }
  z_hat = x_temp %*% loading_est %*% t(loading_est)
  temp = qr(dta$A_tilde)
  AR_oracle = qr.Q(temp)
  z_hat_oracle = x_temp_oracle %*% AR_oracle %*% t(AR_oracle)
  
  pred_error_joint = 0
  pred_error_oracle = 0
  for (j in 1:m) {
    lower_idx = ifelse(j == 1, 0, sum(ds[1:(j - 1)]))
    idx_range = (lower_idx + 1):(lower_idx + ds[j])
    temp = Exp_sphere(z_hat[,idx_range], mu_hat[[j]])
    pred_error_joint = pred_error_joint + 
      sum(acos(c(diag(X_nless[,idx_range] %*% t(temp)))))
    temp = Exp_sphere(z_hat_oracle[,idx_range], dta$mus[[j]])
    pred_error_oracle = pred_error_oracle + 
      sum(acos(c(diag(X_nless[,idx_range] %*% t(temp)))))
  }
  pred_error_joint = pred_error_joint / (n * m)
  pred_error_oracle = pred_error_oracle / (n * m)
  
  # separate modeling
  dist_space_each_separate = rep(NA, m)
  x_pred_oos = matrix(NA, nrow = n, ncol = sum(ds))
  pred_error_separate = 0
  for (j in 1:m) {
    lower_idx = ifelse(j == 1, 0, sum(ds[1:(j - 1)]))
    idx_range = (lower_idx + 1):(lower_idx + ds[j])
    model = LYB_fm(trans_x[,idx_range], r, h)
    temp = parallel_transport(model$V, mu_hat[[j]], dta$mu[[j]])
    dist_space_each_separate[j] = subspace_loss(temp, dta$A_tilde[idx_range,])
    
    pred_error_separate = pred_error_separate +
      oos_eval(X_test[,idx_range], model$V, mu_hat[[j]], X_nless[,idx_range])
  }
  
  pred_error_separate = pred_error_separate / (n * m)

  
  return (list("mu_loss" = mu_loss, "D_joint" = dist_space_joint,
               "D_each_joint" = dist_space_each_joint,
               "D_each_separate" = dist_space_each_separate,
               "pred_err_joint" = pred_error_joint,
               "pred_err_separate" = pred_error_separate,
               "pred_err_oracle" = pred_error_oracle))
}

# Simulation
max.sim = 1000
n = 300
m = 5
ds = rep(5, m)
r = 3
snr = 0.7
spec = 1
if (spec == 1) {
  A = 2
} else if (spec == 2) {
  A = 0.5
} else if (spec == 3) {
  A = 0.1
}
b1 = snr * (1 - 0.8) * pi / (A * r)
b2 = (1 - snr) * (1 - 0.8) * pi / (A * r)

dist_loss = rep(0, max.sim)
sub_loss_J = matrix(0, nrow = max.sim, ncol = m)
sub_loss_S = matrix(0, nrow = max.sim, ncol = m)
PE_J = rep(0, max.sim)
PE_O = rep(0, max.sim)
PE_S = rep(0, max.sim)

for (sim in 1:max.sim) {
  res = simulation(n, ds, r, b1, b2, spec, mu_tol = 1e-4)
  
  dist_loss[sim] = res$D_joint
  sub_loss_J[sim,] = res$D_each_joint
  sub_loss_S[sim,] = res$D_each_separate
  PE_J[sim] = res$pred_err_joint
  PE_S[sim] = res$pred_err_separate
  PE_O[sim] = res$pred_err_oracle
  
  if (sim > 1 && sim %% 100 == 0) {
    cat("iteration", sim, "\n")
  }
}

mean(dist_loss); sd(dist_loss)
mean(PE_J); sd(PE_J)
mean(PE_S); sd(PE_S)
mean(PE_O); sd(PE_O)

colMeans(sub_loss_J); apply(sub_loss_J, 2, sd)
colMeans(sub_loss_S); apply(sub_loss_S, 2, sd)






