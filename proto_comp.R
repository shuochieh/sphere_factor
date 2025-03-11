source("./main_func.R")

z_max = rep(0, 300 - 1)
z_mean = rep(0, 300 - 1)
z_min = rep(0, 300 - 1)
for (i in 2:300) {
  x = rep(1, i) / sqrt(i)
  y = matrix(abs(rnorm(i * 3000)), ncol = i)
  y = y / apply(y, 1, norm, "2")
  z_max[i-1] = max(geod_sphere(x, y))
  z_mean[i-1] = mean(geod_sphere(x, y))
  z_min[i-1] = min(geod_sphere(x, y))
  cat("iter", i, "\n")
}

plot(z_max, type = "b", ylim = c(0, pi / 2), xlab = "dimension of the sphere",
     ylab = "geodesic distance")
lines(z_mean, type = "b", col = 2)
lines(z_min, type = "b", col = 3)
legend("topright", legend = c("Maximum", "Mean", "Minimum"), 
       col = c(1, 2, 3), lty = 1, pch = 1)

z = rep(0, 300 - 1)
for (i in 2:300) {
  z[i - 1] = geod_sphere(rep(1, i) / sqrt(i),
                         c(1, rep(0, i - 1)))
  cat("iter", i, "\n")
}
plot(z, type = "b", ylim = c(0, pi / 2), xlab = "dimension of the sphere",
     ylab = "geodesic distance")
abline(h = pi / 2, col = 2)


dta_gen_composition = function (n, p, type, geo_rng = 0.3, d = 1, r = 1) {

  factors = array(NA, dim = c(n, r))
  for (i in 1:r) {
    if (type %in% c(5, 6, 7)) {
      factors[,i] = arima.sim(list(ar = c(0.8)), n = n, sd = geo_rng / 2)
    }
  }
  
  if (type %in% c(5, 6, 7)) {
    A = array(rnorm((p * d - d) * r), dim = c(p * d - d, r))
    for (i in 1:r) {
      if (i > 1) {
        for (j in 1:(i - 1)) {
          A[,i] = perp_proj(A[,i], A[,j])
        }
      }
      A[,i] = A[,i] / norm(A[,i], "2")
    }
  }
  
  z = array(NA, dim = c(n, p * d))
  x = vector("list", d)
  V = array(NA, dim = c(d * p, r))
  for (j in 1:d) {
    v = matrix(rnorm(p * (p - 1)), nrow = p)
    v = t(t(v) / apply(v, 2, norm, "2"))
    for (i in 1:(p - 1)) {
      v[,i] = perp_proj(v[,i], rep(1, p) / sqrt(p))
    }
    temp_z = factors %*% t(A[((j - 1) * (p - 1) + 1):(j * (p - 1)),]) %*% t(v) / sqrt(r)
    z[,((j - 1) * p + 1):(j * p)] = temp_z
    x[[j]] = Exp_sphere(temp_z, rep(1, p) / sqrt(p))
    V[((j - 1) * p + 1):(j * p),] = v %*% A[((j - 1) * (p - 1) + 1):(j * (p - 1)),]
  }
  
  return (list("x" = x, "z" = z, "factors" = factors, "v" = V))
}

# dta_gen_composition = function (n, p, type, geo_rng = 0.3, d = NULL) {
#   if (type == 1) {
#     mu = rep(1, p) / sqrt(p)
#     factors = seq(from = -geo_rng, to = geo_rng, length.out = n)
#     
#     v = rnorm(p)
#     v = perp_proj(v, mu)
#     v = v / norm(v, "2")
#     
#     z = factors %o% v 
#     
#     x = Exp_sphere(z, mu)
#   }
#   
#   if (type == 2) {
#     if (is.null(d)) {
#       stop("d should be provided")
#     }
#     
#     mu = vector("list", length = d)
#     x = vector("list", length = d)
#     z = vector("list", length = d)
#     v = vector("list", length = d)
#     
#     factors = seq(from = -geo_rng, to = geo_rng, length.out = n)
#     for (i in 1:d) {
#       mu[[i]] = rep(1, p) / sqrt(p)
# 
#       v[[i]] = rnorm(p)
#       v[[i]] = perp_proj(v[[i]], mu[[i]])
#       v[[i]] = v[[i]] / norm(v[[i]], "2")
#       
#       z[[i]] = factors %o% v[[i]] 
#       
#       x[[i]] = Exp_sphere(z[[i]], mu[[i]])
#     }
#   }
#   
#   if (type == 3) {
#     mu = rep(1, p) / sqrt(p)
#     factors = matrix(0, nrow = n, ncol = 2)
#     factors[,1] = seq(from = -geo_rng, to = geo_rng, length.out = n)
#     factors[,2] = geo_rng * 2 * ((1/(1+exp(-seq(from = -6, to = 6, length.out = n)))) - 0.5)
#       # sin(seq(from = 0, to = 1, length.out = n) * 2 * pi)
#       # c(arima.sim(list(order = c(1,0,0), ar = 0.5), n = n, sd = 0.2))
#     
#     v = matrix(rnorm(2 * p), ncol = 2)
#     v[,1] = perp_proj(v[,1], mu)
#     v[,1] = v[,1] / norm(v[,1], "2")
#     v[,2] = perp_proj(v[,2], mu)
#     v[,2] = v[,2] - c(v[,2] %*% v[,1]) * v[,1] / norm(v[,1], "2")
#     v = t(t(v) / apply(v, 2, norm, "2"))
#     
#     z = factors[,1] %o% v[,1] + factors[,2] %o% v[,2]
#     
#     x = Exp_sphere(z, mu)
#   }
#   
#   if (type == 4) {
#     if (is.null(d)) {
#       stop("d should be provided")
#     }
#     
#     mu = vector("list", length = d)
#     x = vector("list", length = d)
#     z = vector("list", length = d)
#     v = vector("list", length = d)
#     
#     for (i in 1:d) {
#       mu[[i]] = rep(1, p) / sqrt(p)
#       factors = matrix(0, nrow = n, ncol = 2)
#       factors[,1] = seq(from = -geo_rng, to = geo_rng, length.out = n)
#       factors[,2] = geo_rng * 2 * ((1/(1+exp(-seq(from = -6, to = 6, length.out = n)))) - 0.5)
#       
#       v[[i]] = matrix(rnorm(2 * p), ncol = 2)
#       v[[i]][,1] = perp_proj(v[[i]][,1], mu[[i]])
#       v[[i]][,1] = v[[i]][,1] / norm(v[[i]][,1], "2")
#       v[[i]][,2] = perp_proj(v[[i]][,2], mu[[i]])
#       v[[i]][,2] = v[[i]][,2] - c(v[[i]][,2] %*% v[[i]][,1]) * v[[i]][,1] / norm(v[[i]][,1], "2")
#       v[[i]] = t(t(v[[i]]) / apply(v[[i]], 2, norm, "2"))
#       
#       z[[i]] = factors[,1] %o% v[[i]][,1] + factors[,2] %o% v[[i]][,2]
# 
#       x[[i]] = Exp_sphere(z[[i]], mu[[i]])
#     }
#   }
#   
#   if (type == 5) {
#     mu = rep(1, p) / sqrt(p)
#     factors = arima.sim(list(ar = c(0.8)), n = n, sd = geo_rng / 2)
#     
#     v = rnorm(p)
#     v = perp_proj(v, mu)
#     v = v / norm(v, "2")
#     
#     z = factors %o% v 
#     
#     x = Exp_sphere(z, mu)
#   }
#   
#   if (type == 6) {
#     if (is.null(d)) {
#       stop("d should be provided")
#     }
#     
#     mu = vector("list", length = d)
#     x = vector("list", length = d)
#     z = vector("list", length = d)
#     v = vector("list", length = d)
# 
#     factors = arima.sim(list(ar = c(0.8)), n = n, sd = geo_rng / 2)
#     
#     for (i in 1:d) {
#       mu[[i]] = rep(1, p) / sqrt(p)
# 
#       v[[i]] = rnorm(p)
#       v[[i]] = perp_proj(v[[i]], mu)
#       v[[i]] = v[[i]] / norm(v[[i]], "2")
#       
#       z[[i]] = factors %o% v[[i]] 
#       
#       x[[i]] = Exp_sphere(z[[i]], mu[[i]])
#     }
#   }
#   
#   if (type == 7) {
#     if (is.null(d)) {
#       stop("d should be provided")
#     }
#     
#     mu = vector("list", length = d)
#     x = vector("list", length = d)
#     z = vector("list", length = d)
#     v = vector("list", length = d)
#     
#     factors = array(NA, dim = c(n, d)) 
#     for (j in 1:d) {
#       factors[,j] = arima.sim(list(ar = c(0.8)), n = n, sd = geo_rng / 2)
#     }
#     
#     for (i in 1:d) {
#       mu[[i]] = rep(1, p) / sqrt(p)
#       
#       v[[i]] = array(rnorm(d * p), dim = c(d, p))
#       for (j in 1:d) {
#         v[[i]][j,] = perp_proj(v[[i]][j,], mu[[i]])
#         if (j > 1) {
#           for (k in 1:(j - 1)) {
#             v[[i]][j,] = perp_proj(v[[i]][j,], v[[i]][k,])
#           }
#         }
#       }
#       v[[i]] = v[[i]] / apply(v[[i]], 1, norm, "2")
#       
#       for (j in 1:d) {
#         if (j == 1) {
#           z[[i]] = (1 / sqrt(d)) * factors[,j] %o% v[[i]][j,]
#         } else{
#           z[[i]] = z[[i]] + (1 / d) * factors[,j] %o% v[[i]][j,]
#         }
#       }
#       
#       x[[i]] = Exp_sphere(z[[i]], mu[[i]])
#     }
#   }
#   
#   
#   return (list("x" = x, "z" = z, "factors" = factors, "v" = v, "mu" = mu))
# }

TV_calculator = function (x) {
  x = t(t(x) - colMeans(x))
  
  return (sum(x^2))
}

matrix_dist = function (A, B) {
  r1 = ncol(A)
  k1 = nrow(A)
  r2 = ncol(B)
  k2 = nrow(B)
  
  if (!is.matrix(A)) {
    r1 = 1
  }
  
  temp = sum(diag(A %*% t(A) %*% B %*% t(B)))
  
  return (sqrt(1 - (1 / r1) * temp))
}

geo_rng = 0.3

# One sphere - one AR factor
par(mfcol = c(3, 6))
for (p in 4 * 2^c(0:5)) {
  sim_number = 30
  
  collect_FVU = array(NA, dim = c(sim_number, 4))
  collect_PE = array(NA, dim = c(sim_number, 4))
  collect_FVU_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, 4))
  collect_TV = array(NA, dim = c(sim_number, 1))
  
  for (zz in 1:sim_number) {
    
    temp = dta_gen_composition(200, p, 5, geo_rng = geo_rng, d = 1, r = 1)
    x = temp$x
    res = main_uneven_sphere(x, r = 4, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    collect_PE[zz,] = c(res$pe_e * 100)
    collect_PE_linear[zz,] = c(res$pe_e_linear * 100)
    collect_PE_linear_direct[zz,] = c(res$pe_e_linear_direct * 100)
    collect_TV[zz,] = sqrt(TV_calculator(x[[1]][151:200,]^2)) * 100
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", ylim = c(0, 0.05),
       main = paste("FVU; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
  plot(x = colMeans(collect_PE), type = "b", ylim = c(0, 0.05 * 100),
       main = paste("PE; p =", p, "; TV = ",
                    round(mean(c(collect_TV)), 2), "%"),
       ylab = "Prediction error (%)")
  lines(colMeans(collect_PE_linear), type = "b", col = 2)
  lines(colMeans(collect_PE_linear_direct), type = "b", col = 3)
  
  plot(x = colMeans(collect_PE) / mean(c(collect_TV)), type = "b", ylim = c(0, 0.1),
       main = paste("PE; p =", p, "; as a fraction of TV"),
       ylab = "")
  lines(colMeans(collect_PE_linear) / mean(c(collect_TV)), type = "b", col = 2)
  lines(colMeans(collect_PE_linear_direct) / mean(c(collect_TV)), type = "b", col = 3)
  
}

par(mfrow = c(2, 3))
for (p in 4 * 2^c(0:5)) {
  temp = dta_gen_composition(200, p, 5, geo_rng = geo_rng, d = 1, r = 1)
  x = temp$x[[1]]
  mu_hat = mean_on_sphere(x)
  plot(geod_sphere(x, mu_hat), type = "l",
       xlab = "Time", ylab = "Geo Dist to mean",
       main = paste("p =", p))
}

par(mfrow = c(2, 3))
for (p in 4 * 2^c(0:5)) {
  sim_number = 30
  
  collect_Vdist = array(NA, dim = c(sim_number, min(p, 8)))

  for (zz in 1:sim_number) {
    
    temp = dta_gen_composition(200, p, 5, geo_rng = geo_rng, d = 1, r = 1)
    x = temp$x
    res = main_uneven_sphere(x, r = min(p, 8), test_size = 50, h = 3)
    
    for (j in 1:min(p, 8)) {
      collect_Vdist[zz,j] = matrix_dist(res$V[,1:j], res$V_linear[,1:j])
    }

    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_Vdist), type = "b", ylim = c(0, 1),
       main = paste("Dist. btw Vs ; p =", p))
}

# Product spheres - one AR factor
par(mfrow = c(1, 3))
p = 4
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_FVU = array(NA, dim = c(sim_number, 15))
  collect_FVU_linear = array(NA, dim = c(sim_number, 15))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 6, d = d, geo_rng = geo_rng, r = 1)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", 
       ylim = c(0, 0.008),
       main = paste("d = ", d, "; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
}

for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 15))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 6, d = d, geo_rng = geo_rng, r = 1)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }
  
  # par(mfrow = c(d / 5, 5))  
  # for (j in 1:d) {
  #   mu_hat = mean_on_sphere(temp$x[[j]])
  #   plot(geod_sphere(temp$x[[j]], mu_hat), type = "l",
  #        xlab = "Time", ylab = "Geo Dist to mean",
  #        main = paste(j, "/", d))
  # }  
  
  par(mfrow = c(d / 5, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]), type = "b", 
         ylab = "prediction error (%)",
         ylim = c(0, 0.5),
         main = paste("d = ", d, "; TV = ",
                      round(mean(collect_TV[,j]), 2), "%")
    )
    lines(colMeans(collect_PE_linear[,j,]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]), type = "b", col = 3)
  }
}

for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 15))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 6, d = d, geo_rng = geo_rng, r = 1)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }
  
  par(mfrow = c(d / 5, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]) / mean(collect_TV[,j]), type = "b", 
         ylab = "",
         ylim = c(0, 0.02),
         main = paste("d = ", d, "; as a fraction of TV")
    )
    lines(colMeans(collect_PE_linear[,j,]) / mean(collect_TV[,j]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]) / mean(collect_TV[,j]), type = "b", col = 3)
  }
}

# One sphere - two AR factors
par(mfcol = c(3, 6))
for (p in 4 * 2^c(0:5)) {
  sim_number = 30
  
  collect_FVU = array(NA, dim = c(sim_number, 4))
  collect_PE = array(NA, dim = c(sim_number, 4))
  collect_FVU_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, 4))
  collect_TV = array(NA, dim = c(sim_number, 1))
  
  for (zz in 1:sim_number) {
    
    temp = dta_gen_composition(200, p, 5, geo_rng = geo_rng, d = 1, r = 2)
    x = temp$x
    res = main_uneven_sphere(x, r = 4, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    collect_PE[zz,] = c(res$pe_e * 100)
    collect_PE_linear[zz,] = c(res$pe_e_linear * 100)
    collect_PE_linear_direct[zz,] = c(res$pe_e_linear_direct * 100)
    collect_TV[zz,] = sqrt(TV_calculator(x[[1]][151:200,]^2)) * 100
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", ylim = c(0, 0.5),
       main = paste("FVU; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
  plot(x = colMeans(collect_PE), type = "b", ylim = c(0, 0.1 * 100),
       main = paste("PE; p =", p, "; TV = ",
                    round(mean(c(collect_TV)), 2), "%"),
       ylab = "Prediction error (%)")
  lines(colMeans(collect_PE_linear), type = "b", col = 2)
  lines(colMeans(collect_PE_linear_direct), type = "b", col = 3)
  
  plot(x = colMeans(collect_PE) / mean(c(collect_TV)), type = "b", ylim = c(0, 0.12),
       main = paste("PE; p =", p, "; as a fraction of TV"),
       ylab = "")
  lines(colMeans(collect_PE_linear) / mean(c(collect_TV)), type = "b", col = 2)
  lines(colMeans(collect_PE_linear_direct) / mean(c(collect_TV)), type = "b", col = 3)
  
}

par(mfrow = c(2, 3))
for (p in 4 * 2^c(0:5)) {
  sim_number = 30
  
  collect_Vdist = array(NA, dim = c(sim_number, min(p, 8)))
  
  for (zz in 1:sim_number) {
    
    temp = dta_gen_composition(200, p, 5, geo_rng = geo_rng, d = 1, r = 2)
    x = temp$x
    res = main_uneven_sphere(x, r = min(p, 8), test_size = 50, h = 3)
    
    for (j in 1:min(p, 8)) {
      collect_Vdist[zz,j] = matrix_dist(res$V[,1:j], res$V_linear[,1:j])
    }
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_Vdist), type = "b", ylim = c(0, 1),
       main = paste("Dist. btw Vs ; p =", p))
}


# Product spheres - two AR factors
par(mfrow = c(1, 3))
p = 4
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_FVU = array(NA, dim = c(sim_number, 15))
  collect_FVU_linear = array(NA, dim = c(sim_number, 15))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 7, d = d, geo_rng = geo_rng, r = 2)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", 
       ylim = c(0, 0.5),
       main = paste("d = ", d, "; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
}

par(mfrow = c(2, 4))
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 15))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 7, d = d, geo_rng = geo_rng, r = 2)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }
  
  # par(mfrow = c(d / 5, 5))  
  # for (j in 1:d) {
  #   mu_hat = mean_on_sphere(temp$x[[j]])
  #   plot(geod_sphere(temp$x[[j]], mu_hat), type = "l",
  #        xlab = "Time", ylab = "Geo Dist to mean",
  #        main = paste(j, "/", d))
  # }  
  
  par(mfrow = c(d / 5, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]), type = "b", 
         ylab = "prediction error (%)",
         ylim = c(0, 6),
         main = paste("d = ", d, "; TV = ",
                      round(mean(collect_TV[,j]), 2), "%")
    )
    lines(colMeans(collect_PE_linear[,j,]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]), type = "b", col = 3)
  }
}

for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 15))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 6, d = d, geo_rng = geo_rng, r = 2)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }
  
  par(mfrow = c(d / 5, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]) / mean(collect_TV[,j]), type = "b", 
         ylab = "",
         ylim = c(0, 0.15),
         main = paste("d = ", d, "; as a fraction of TV")
    )
    lines(colMeans(collect_PE_linear[,j,]) / mean(collect_TV[,j]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]) / mean(collect_TV[,j]), type = "b", col = 3)
  }
}

# One spheres - 8 AR factors
par(mfcol = c(3, 5))
for (p in 4 * 2^c(1:5)) {
  sim_number = 30
  
  collect_FVU = array(NA, dim = c(sim_number, min(15, p)))
  collect_PE = array(NA, dim = c(sim_number, min(15, p)))
  collect_FVU_linear = array(NA, dim = c(sim_number, min(15, p)))
  collect_PE_linear = array(NA, dim = c(sim_number, min(15, p)))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, min(15, p)))
  collect_TV = array(NA, dim = c(sim_number, 1))
  
  for (zz in 1:sim_number) {
    
    temp = dta_gen_composition(200, p, 5, geo_rng = geo_rng, d = 1, r = 8)
    x = temp$x
    res = main_uneven_sphere(x, r = min(15, p), test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    collect_PE[zz,] = c(res$pe_e * 100)
    collect_PE_linear[zz,] = c(res$pe_e_linear * 100)
    collect_PE_linear_direct[zz,] = c(res$pe_e_linear_direct * 100)
    collect_TV[zz,] = sqrt(TV_calculator(x[[1]][151:200,]^2)) * 100
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", ylim = c(0, 0.8),
       main = paste("FVU; p =", p, "; r =", 8))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
  plot(x = colMeans(collect_PE), type = "b", ylim = c(0, 0.1 * 100),
       main = paste("PE; p =", p, "; TV = ",
                    round(mean(c(collect_TV)), 2), "%"),
       ylab = "Prediction error (%)")
  lines(colMeans(collect_PE_linear), type = "b", col = 2)
  lines(colMeans(collect_PE_linear_direct), type = "b", col = 3)
  
  plot(x = colMeans(collect_PE) / mean(c(collect_TV)), type = "b", ylim = c(0, 0.2),
       main = paste("PE; p =", p, "; as a fraction of TV"),
       ylab = "")
  lines(colMeans(collect_PE_linear) / mean(c(collect_TV)), type = "b", col = 2)
  lines(colMeans(collect_PE_linear_direct) / mean(c(collect_TV)), type = "b", col = 3)
  
}

par(mfrow = c(2, 3))
for (p in 4 * 2^c(1:5)) {
  sim_number = 30
  
  collect_Vdist = array(NA, dim = c(sim_number, min(p, 10)))
  
  for (zz in 1:sim_number) {
    
    temp = dta_gen_composition(200, p, 5, geo_rng = geo_rng, d = 1, r = 8)
    x = temp$x
    res = main_uneven_sphere(x, r = min(p, 10), test_size = 50, h = 3)
    
    for (j in 1:min(p, 10)) {
      collect_Vdist[zz,j] = matrix_dist(res$V[,1:j], res$V_linear[,1:j])
    }
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_Vdist), type = "b", ylim = c(0, 1),
       main = paste("Dist. btw Vs ; p =", p))
}

# Product spheres - 8 AR factors
par(mfrow = c(1, 3))
p = 4
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_FVU = array(NA, dim = c(sim_number, 15))
  collect_FVU_linear = array(NA, dim = c(sim_number, 15))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 7, d = d, geo_rng = geo_rng, r = 8)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", 
       ylim = c(0, 0.8),
       main = paste("d = ", d, "; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
}

par(mfrow = c(2, 4))
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 15))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 7, d = d, geo_rng = geo_rng, r = 8)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }
  
  par(mfrow = c(d / 5, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]), type = "b", 
         ylab = "prediction error (%)",
         ylim = c(0, 6),
         main = paste("d = ", d, "; TV = ",
                      round(mean(collect_TV[,j]), 2), "%")
    )
    lines(colMeans(collect_PE_linear[,j,]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]), type = "b", col = 3)
  }
}

for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 15))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 6, d = d, geo_rng = geo_rng, r = 8)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }
  
  par(mfrow = c(1, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]) / mean(collect_TV[,j]), type = "b", 
         ylab = "",
         ylim = c(0, 0.15),
         main = paste("d = ", d, "; as a fraction of TV")
    )
    lines(colMeans(collect_PE_linear[,j,]) / mean(collect_TV[,j]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]) / mean(collect_TV[,j]), type = "b", col = 3)
  }
}

















# type = 1
par(mfcol = c(2, 6))
for (p in 4 * 2^c(0:5)) {
  sim_number = 30
  
  collect_FVU = array(NA, dim = c(sim_number, 4))
  collect_PE = array(NA, dim = c(sim_number, 4))
  collect_FVU_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, 4))
  collect_TV = array(NA, dim = c(sim_number, 1))
  
  for (zz in 1:sim_number) {

    temp = dta_gen_composition(200, p, 1, geo_rng = geo_rng)
    x = temp$x
    res = main_uneven_sphere(list("x" = x), r = 4, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    collect_PE[zz,] = c(res$pe_e * 100)
    collect_PE_linear[zz,] = c(res$pe_e_linear * 100)
    collect_PE_linear_direct[zz,] = c(res$pe_e_linear_direct * 100)
    collect_TV[zz,] = sqrt(TV_calculator(x[151:200,]^2)) * 100
  }
  plot(colMeans(collect_FVU), type = "b", ylim = c(0, 0.004),
       main = paste("FVU; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
  plot(x = colMeans(collect_PE), type = "b", ylim = c(0, 0.025 * 100),
       main = paste("PE; p =", p, "; TV = ",
                    round(mean(c(collect_TV)), 2), "%"),
       ylab = "Prediction error (%)")
  lines(colMeans(collect_PE_linear), type = "b", col = 2)
  lines(colMeans(collect_PE_linear_direct), type = "b", col = 3)
}

par(mfrow = c(2, 3))
for (p in 4 * 2^c(0:5)) {
  temp = dta_gen_composition(200, p, 1, geo_rng = geo_rng)
  x = temp$x
  mu_hat = mean_on_sphere(x)
  plot(geod_sphere(x, mu_hat), type = "l",
       xlab = "Time", ylab = "Geo Dist to mean",
       main = paste("p =", p))
}

# type = 3
par(mfcol = c(2, 6))
for (p in 4 * 2^c(0:5)) {
  sim_number = 100
  
  collect_FVU = array(NA, dim = c(sim_number, 4))
  collect_PE = array(NA, dim = c(sim_number, 4))
  collect_FVU_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, 4))
  collect_TV = array(NA, dim = c(sim_number, 1))
  
  for (zz in 1:sim_number) {
    
    temp = dta_gen_composition(200, p, 3, geo_rng = geo_rng)
    x = temp$x
    res = main_uneven_sphere(list("x" = x), r = 4, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    collect_PE[zz,] = c(res$pe_e * 100)
    collect_PE_linear[zz,] = c(res$pe_e_linear * 100)
    collect_PE_linear_direct[zz,] = c(res$pe_e_linear_direct * 100)
    collect_TV[zz,] = sqrt(TV_calculator(x[151:200,]^2)) * 100
  }
  plot(colMeans(collect_FVU), type = "b", ylim = c(0, 0.05),
       main = paste("FVU; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
  plot(x = colMeans(collect_PE), type = "b", ylim = c(0, 0.05 * 100),
       main = paste("PE; p =", p, "; TV = ",
                    round(mean(c(collect_TV)), 2), "%"),
       ylab = "Prediction error (%)")
  lines(colMeans(collect_PE_linear), type = "b", col = 2)
  lines(colMeans(collect_PE_linear_direct), type = "b", col = 3)
}

par(mfrow = c(2, 3))
for (p in 4 * 2^c(0:5)) {
  temp = dta_gen_composition(200, p, 3, geo_rng = geo_rng)
  x = temp$x
  mu_hat = mean_on_sphere(x)
  plot(geod_sphere(x, mu_hat), type = "l",
       xlab = "Time", ylab = "Geo Dist to mean",
       main = paste("p =", p))
}


# type = 2
par(mfrow = c(1, 3))
p = 4
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_FVU = array(NA, dim = c(sim_number, 5))
  collect_FVU_linear = array(NA, dim = c(sim_number, 5))

  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 2, d = d, geo_rng = geo_rng)
    res = main_uneven_sphere(temp$x, 5, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", 
       ylim = c(0, 0.01),
       main = paste("d = ", d, "; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
}

par(mfrow = c(2, 4))
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 5))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 5))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 5))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 2, d = d, geo_rng = geo_rng)
    res = main_uneven_sphere(temp$x, 5, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }

  # par(mfrow = c(d / 5, 5))  
  # for (j in 1:d) {
  #   mu_hat = mean_on_sphere(temp$x[[j]])
  #   plot(geod_sphere(temp$x[[j]], mu_hat), type = "l",
  #        xlab = "Time", ylab = "Geo Dist to mean",
  #        main = paste(j, "/", d))
  # }  
  
  par(mfrow = c(d / 5, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]), type = "b", 
         ylab = "prediction error (%)",
         ylim = c(0, 4),
         main = paste("d = ", d, "; TV = ",
                      round(mean(collect_TV[,j]), 2), "%")
        )
    lines(colMeans(collect_PE_linear[,j,]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]), type = "b", col = 3)
  }
}

# type = 4
par(mfrow = c(1, 3))
p = 4
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_FVU = array(NA, dim = c(sim_number, 5))
  collect_FVU_linear = array(NA, dim = c(sim_number, 5))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 4, d = d, geo_rng = geo_rng)
    res = main_uneven_sphere(temp$x, 5, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", 
       ylim = c(0, 0.03),
       main = paste("d = ", d, "; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
}

par(mfrow = c(2, 4))
p = 4
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 5))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 5))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 5))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 4, d = d, geo_rng = geo_rng)
    res = main_uneven_sphere(temp$x, 5, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }
  
  # par(mfrow = c(d / 5, 5))  
  # for (j in 1:d) {
  #   mu_hat = mean_on_sphere(temp$x[[j]])
  #   plot(geod_sphere(temp$x[[j]], mu_hat), type = "l",
  #        xlab = "Time", ylab = "Geo Dist to mean",
  #        main = paste(j, "/", d))
  # }  
  
  par(mfrow = c(d / 5, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]), type = "b", 
         ylab = "prediction error (%)",
         ylim = c(0, 10),
         main = paste("d = ", d, "; TV = ",
                      round(mean(collect_TV[,j]), 2), "%")
    )
    lines(colMeans(collect_PE_linear[,j,]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]), type = "b", col = 3)
  }
}





# type = 5
par(mfcol = c(2, 6))
for (p in 4 * 2^c(0:5)) {
  sim_number = 30
  
  collect_FVU = array(NA, dim = c(sim_number, 4))
  collect_PE = array(NA, dim = c(sim_number, 4))
  collect_FVU_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear = array(NA, dim = c(sim_number, 4))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, 4))
  collect_TV = array(NA, dim = c(sim_number, 1))
  
  for (zz in 1:sim_number) {
    
    temp = dta_gen_composition(200, p, 5, geo_rng = geo_rng)
    x = temp$x
    res = main_uneven_sphere(list("x" = x), r = 4, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    collect_PE[zz,] = c(res$pe_e * 100)
    collect_PE_linear[zz,] = c(res$pe_e_linear * 100)
    collect_PE_linear_direct[zz,] = c(res$pe_e_linear_direct * 100)
    collect_TV[zz,] = sqrt(TV_calculator(x[151:200,]^2)) * 100
  }
  plot(colMeans(collect_FVU), type = "b", ylim = c(0, 0.004),
       main = paste("FVU; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
  plot(x = colMeans(collect_PE), type = "b", ylim = c(0, 0.025 * 100),
       main = paste("PE; p =", p, "; TV = ",
                    round(mean(c(collect_TV)), 2), "%"),
       ylab = "Prediction error (%)")
  lines(colMeans(collect_PE_linear), type = "b", col = 2)
  lines(colMeans(collect_PE_linear_direct), type = "b", col = 3)
}

par(mfrow = c(2, 3))
for (p in 4 * 2^c(0:5)) {
  temp = dta_gen_composition(200, p, 5, geo_rng = geo_rng)
  x = temp$x
  mu_hat = mean_on_sphere(x)
  plot(geod_sphere(x, mu_hat), type = "l",
       xlab = "Time", ylab = "Geo Dist to mean",
       main = paste("p =", p))
}

# type = 6
par(mfrow = c(1, 3))
p = 4
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_FVU = array(NA, dim = c(sim_number, 15))
  collect_FVU_linear = array(NA, dim = c(sim_number, 15))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 6, d = d, geo_rng = geo_rng)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", 
       ylim = c(0, 0.8),
       main = paste("d = ", d, "; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
}

par(mfrow = c(2, 4))
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 15))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 6, d = d, geo_rng = geo_rng)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }
  
  # par(mfrow = c(d / 5, 5))  
  # for (j in 1:d) {
  #   mu_hat = mean_on_sphere(temp$x[[j]])
  #   plot(geod_sphere(temp$x[[j]], mu_hat), type = "l",
  #        xlab = "Time", ylab = "Geo Dist to mean",
  #        main = paste(j, "/", d))
  # }  
  
  par(mfrow = c(d / 5, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]), type = "b", 
         ylab = "prediction error (%)",
         ylim = c(0, 20),
         main = paste("d = ", d, "; TV = ",
                      round(mean(collect_TV[,j]), 2), "%")
    )
    lines(colMeans(collect_PE_linear[,j,]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]), type = "b", col = 3)
  }
}

# type = 7
par(mfrow = c(1, 3))
p = 4
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_FVU = array(NA, dim = c(sim_number, 15))
  collect_FVU_linear = array(NA, dim = c(sim_number, 15))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 7, d = d, geo_rng = geo_rng)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_FVU[zz,] = res$FVU_e
    collect_FVU_linear[zz,] = res$FVU_e_linear
    
    cat("iteration", zz, "\n")
  }
  plot(colMeans(collect_FVU), type = "b", 
       ylim = c(0, 0.5),
       main = paste("d = ", d, "; p =", p))
  lines(colMeans(collect_FVU_linear), type = "b", col = 2)
}

par(mfrow = c(2, 4))
for (d in c(5 * c(1:3))) {
  sim_number = 30
  collect_PE = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear = array(NA, dim = c(sim_number, d, 15))
  collect_PE_linear_direct = array(NA, dim = c(sim_number, d, 15))
  collect_TV = array(NA, dim = c(sim_number, d))
  
  for (zz in 1:sim_number) {
    temp = dta_gen_composition(200, p, 7, d = d, geo_rng = geo_rng)
    res = main_uneven_sphere(temp$x, 15, test_size = 50, h = 3)
    
    collect_PE[zz,,] = res$pe_e * 100
    collect_PE_linear[zz,,] = res$pe_e_linear * 100
    collect_PE_linear_direct[zz,,] = res$pe_e_linear_direct * 100
    for (j in 1:d) {
      collect_TV[zz,j] = sqrt(TV_calculator(temp$x[[j]][151:200,]^2)) * 100
    }
    
    cat("iteration", zz, "\n")
  }
  
  # par(mfrow = c(d / 5, 5))  
  # for (j in 1:d) {
  #   mu_hat = mean_on_sphere(temp$x[[j]])
  #   plot(geod_sphere(temp$x[[j]], mu_hat), type = "l",
  #        xlab = "Time", ylab = "Geo Dist to mean",
  #        main = paste(j, "/", d))
  # }  
  
  par(mfrow = c(d / 5, 5))  
  for (j in 1:d) {
    plot(colMeans(collect_PE[,j,]), type = "b", 
         ylab = "prediction error (%)",
         ylim = c(0, 6),
         main = paste("d = ", d, "; TV = ",
                      round(mean(collect_TV[,j]), 2), "%")
    )
    lines(colMeans(collect_PE_linear[,j,]), type = "b", col = 2)
    lines(colMeans(collect_PE_linear_direct[,j,]), type = "b", col = 3)
  }
}
