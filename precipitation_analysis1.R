read_precip_data <- function(file) {
  df <- read.csv(file, skip = 3, na.strings = "-99")
  
  colnames(df) <- c("Date", "Value")
  
  return(df)
}

convert_to_date <- function(numeric_date) {
  year <- as.integer(numeric_date / 100)
  month <- sprintf("%02d", numeric_date %% 100)
  return(paste0(year, "-", month, "-01"))
}

# state-station list
# West
ss_ls_pacific = list("CA" = 7, "OR" = 9, "WA" = 10) # Pacific
ss_ls_mountain = list("AZ" = 7, "CO" = 5, "ID" = 10, # Mountain
                      "MT" = 7, "NV" = 4, "NM" = 8,
                      "UT" = 7, "WY" = 10)

# South
ss_ls_sAtlantic = list("FL" = 7, "GA" = 9, "NC" = 8, # South Atlantic
                       "SC" = 7, "VA" = 6, "MD" = 8,
                       "DE" = 2, "WV" = 6) 
ss_ls_ESCentral = list("AL" = 8, "KY" = 4, "MS" = 10, # East South Central
                       "TN" = 4)
ss_ls_WSCentral = list("AR" = 9, "LA" = 9, "OK" = 9, # West South Central
                       "TX" = 10)

# Midwest
ss_ls_ENCentral = list("IL" = 9, "IN" = 9, "MI" = 10, # East North Central
                       "OH" = 10, "WI" = 9)
ss_ls_WNCentral = list("IA" = 9, "KS" = 9, "MN" = 9, # West North Central
                       "MO" = 6, "NE" = 8, "ND" = 9,
                       "SD" = 9)

# Northeast
ss_ls_NewEng = list("CT" = 3, "ME" = 3, "MA" = 3, # New England
                    "NH" = 2, "RI" = 1, "VT" = 3)
ss_ls_MidAtl = list("NY" = 10, "NJ" = 3, "PA" = 10) # Middle Atlantic

# Analysis by region: Create multi-compositional data
comp_data_get = function (ss_ls) {
  p = length(ss_ls)
  for (i in 1:p) {
    state_name = names(ss_ls)[i]
    for (j in 1:ss_ls[[i]]) {
      dir_name = paste0("./precipitation/", state_name, j, ".csv")
      temp = read_precip_data(dir_name)
      if (i == 1 && j == 1) {
        n = nrow(temp)
        Dates = convert_to_date(temp[,1])
        res = matrix(0, nrow = n, ncol = p)
      } 
      res[,i] = res[,i] + temp[,2]
    }
    res[,i] = res[,i] / ss_ls[[i]]
  }
  
  res = res / rowSums(res)
  colnames(res) = names(ss_ls)
  return (list("dta" = res, "Dates" = Dates))
}

x = vector("list", length = 1)
x[[1]] = comp_data_get(ss_ls_pacific)$dta
Dates = comp_data_get(ss_ls_pacific)$Dates
x[[2]] = comp_data_get(ss_ls_mountain)$dta
x[[3]] = comp_data_get(ss_ls_sAtlantic)$dta
x[[4]] = comp_data_get(ss_ls_ESCentral)$dta
x[[5]] = comp_data_get(ss_ls_WSCentral)$dta
x[[6]] = comp_data_get(ss_ls_ENCentral)$dta
x[[7]] = comp_data_get(ss_ls_WNCentral)$dta
x[[8]] = comp_data_get(ss_ls_NewEng)$dta
x[[9]] = comp_data_get(ss_ls_MidAtl)$dta

mu1 = mean_on_sphere(sqrt(x[[1]]))
mu2 = mean_on_sphere(sqrt(x[[2]]))
mu3 = mean_on_sphere(sqrt(x[[3]]))
mu4 = mean_on_sphere(sqrt(x[[4]]))
mu5 = mean_on_sphere(sqrt(x[[5]]))
mu6 = mean_on_sphere(sqrt(x[[6]]))
mu7 = mean_on_sphere(sqrt(x[[7]]))
mu8 = mean_on_sphere(sqrt(x[[8]]))
mu9 = mean_on_sphere(sqrt(x[[9]]))

par(mfrow = c(3, 3))
temp = x[[1]] ; mu_temp = mu1
# plot(tail(temp[,1], 120), type = "l", ylim = c(0, 1),
#      xlab = "", ylab = "", main = "Pacific region")
# for(j in 2:dim(temp)[2]) {
#   lines(tail(temp[,j], 120), col = j)
# }
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geo Dist to mean", main = "Pacific Region")

temp = x[[2]] ; mu_temp = mu2
# plot(tail(temp[,1], 120), type = "l", ylim = c(0, 1),
#      xlab = "", ylab = "", main = "Mountain region")
# for(j in 2:dim(temp)[2]) {
#   lines(tail(temp[,j], 120), col = j)
# }
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geo Dist to mean", main = "Mountain Region")

temp = x[[3]] ; mu_temp = mu3
# plot(tail(temp[,1], 120), type = "l", ylim = c(0, 1),
#      xlab = "", ylab = "", main = "South Atlantic region")
# for(j in 2:dim(temp)[2]) {
#   lines(tail(temp[,j], 120), col = j)
# }
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geo Dist to mean", main = "South Atlantic Region")

temp = x[[4]] ; mu_temp = mu4
# plot(tail(temp[,1], 120), type = "l", ylim = c(0, 1),
#      xlab = "", ylab = "", main = "East South Central region")
# for(j in 2:dim(temp)[2]) {
#   lines(tail(temp[,j], 120), col = j)
# }
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geo Dist to mean", main = "East South Central Region")

temp = x[[5]] ; mu_temp = mu5
# plot(tail(temp[,1], 120), type = "l", ylim = c(0, 1),
#      xlab = "", ylab = "", main = "West South Central region")
# for(j in 2:dim(temp)[2]) {
#   lines(tail(temp[,j], 120), col = j)
# }
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geo Dist to mean", main = "West South Central Region")

temp = x[[6]] ; mu_temp = mu6
# plot(tail(temp[,1], 120), type = "l", ylim = c(0, 1),
#      xlab = "", ylab = "", main = "East North Central region")
# for(j in 2:dim(temp)[2]) {
#   lines(tail(temp[,j], 120), col = j)
# }
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geo Dist to mean", main = "East North Central Region")

temp = x[[7]] ; mu_temp = mu7
# plot(tail(temp[,1], 120), type = "l", ylim = c(0, 1),
#      xlab = "", ylab = "", main = "West North Central region")
# for(j in 2:dim(temp)[2]) {
#   lines(tail(temp[,j], 120), col = j)
# }
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geo Dist to mean", main = "West North Central Region")

temp = x[[8]] ; mu_temp = mu8
# plot(tail(temp[,1], 120), type = "l", ylim = c(0, 1),
#      xlab = "", ylab = "", main = "New England region")
# for(j in 2:dim(temp)[2]) {
#   lines(tail(temp[,j], 120), col = j)
# }
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geo Dist to mean", main = "New England Region")

temp = x[[9]] ; mu_temp = mu9
# plot(tail(temp[,1], 120), type = "l", ylim = c(0, 1),
#      xlab = "", ylab = "", main = "Middle Atlantic region")
# for(j in 2:dim(temp)[2]) {
#   lines(tail(temp[,j], 120), col = j)
# }
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geo Dist to mean", main = "Middle Atlantic Region")


par(mfrow = c(3, 3))
temp = x[[1]] ; mu_temp = mu1
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "Pacific Region")

temp = x[[2]] ; mu_temp = mu2
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "Mountain Region")

temp = x[[3]] ; mu_temp = mu3
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "South Atlantic Region")

temp = x[[4]] ; mu_temp = mu4
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "East South Central Region")

temp = x[[5]] ; mu_temp = mu5
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "West South Central Region")

temp = x[[6]] ; mu_temp = mu6
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "East North Central Region")

temp = x[[7]] ; mu_temp = mu7
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "West North Central Region")

temp = x[[8]] ; mu_temp = mu8
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "New England Region")

temp = x[[9]] ; mu_temp = mu9
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "Middle Atlantic Region")

# Analysis
source("./main_func.R")

for (i in 1:9) {
  x[[i]] = sqrt(x[[i]])
}

results = main_uneven_sphere(x, r = 20, h = 25, test_size = 360)

par(mfrow = c(1, 1))
plot(results$FVU_e, type = "b", xlab = "", ylab = "FVEU")
lines(results$FVU_e_linear, type = "b", col = 2)

par(mfrow = c(3, 3))
plot(results$pe_e[1,], type = "b", xlab = "", ylab = "", main = "Pacific Region",
     ylim = c(0, max(results$pe_e[1,])))
lines(results$pe_e_linear[1,], type = "b", col = 2)
# lines(results$pe_e_linear_direct[1,], type = "b", col = 3)
plot(results$pe_e[2,], type = "b", xlab = "", ylab = "", main = "Mountain Region",
     ylim = c(0, max(results$pe_e[2,])))
lines(results$pe_e_linear[2,], type = "b", col = 2)
# lines(results$pe_e_linear_direct[2,], type = "b", col = 3)

plot(results$pe_e[3,], type = "b", xlab = "", ylab = "", main = "South Arlantic Region",
     ylim = c(0, max(results$pe_e[3,])))
lines(results$pe_e_linear[3,], type = "b", col = 2)
# lines(results$pe_e_linear_direct[3,], type = "b", col = 3)
plot(results$pe_e[4,], type = "b", xlab = "", ylab = "", main = "East South Central Region",
     ylim = c(0, max(results$pe_e[4,])))
lines(results$pe_e_linear[4,], type = "b", col = 2)
# lines(results$pe_e_linear_direct[4,], type = "b", col = 3)
plot(results$pe_e[5,], type = "b", xlab = "", ylab = "", main = "West South Central Region",
     ylim = c(0, max(results$pe_e[5,])))
lines(results$pe_e_linear[5,], type = "b", col = 2)
# lines(results$pe_e_linear_direct[5,], type = "b", col = 3)

plot(results$pe_e[6,], type = "b", xlab = "", ylab = "", main = "East North Central Region",
     ylim = c(0, max(results$pe_e[6,])))
lines(results$pe_e_linear[6,], type = "b", col = 2)
# lines(results$pe_e_linear_direct[6,], type = "b", col = 3)
plot(results$pe_e[7,], type = "b", xlab = "", ylab = "", main = "West North Central Region",
     ylim = c(0, max(results$pe_e[7,])))
lines(results$pe_e_linear[7,], type = "b", col = 2)
# lines(results$pe_e_linear_direct[7,], type = "b", col = 3)

plot(results$pe_e[8,], type = "b", xlab = "", ylab = "", main = "New England Region",
     ylim = c(0, max(results$pe_e[8,])))
lines(results$pe_e_linear[8,], type = "b", col = 2)
# lines(results$pe_e_linear_direct[8,], type = "b", col = 3)
plot(results$pe_e[9,], type = "b", xlab = "", ylab = "", main = "Middle Atlantic Region",
     ylim = c(0, max(results$pe_e[9,])))
lines(results$pe_e_linear[9,], type = "b", col = 2)
# lines(results$pe_e_linear_direct[9,], type = "b", col = 3)









# Analysis as a whole:
large_ls = c(ss_ls_ENCentral, ss_ls_ESCentral, ss_ls_MidAtl, ss_ls_mountain,
             ss_ls_NewEng, ss_ls_pacific, ss_ls_sAtlantic, ss_ls_WNCentral,
             ss_ls_WSCentral)
x = comp_data_get(large_ls)$dta

mu = mean_on_sphere(sqrt(x))
temp = x ; mu_temp = mu
par(mfrow = c(1, 1))
plot(tail(temp[,1], 120), type = "l", ylim = c(0, 0.15),
     xlab = "", ylab = "", main = "All states")
for(j in 2:dim(temp)[2]) {
  lines(tail(temp[,j], 120), col = j)
}
plot(x = as.Date(Dates), y = geod_sphere(sqrt(temp), mu_temp), type = "l", 
     xlab = "", ylab = "Geodesic distance to Frechet mean", main = "All states")
hist(geod_sphere(sqrt(temp), mu_temp), xlab = "Geodesic distance to Frechet mean",
     main = "All states")

x = list("All_states" = sqrt(x))
results = main_uneven_sphere(x, r = 20, h = 13, test_size = 360)
par(mfrow = c(1, 1))
plot(results$FVU_e, type = "b", xlab = "number of factors",
     ylab = "Fraction of (Euclidean) variances unexplained")
lines(results$FVU_e_linear, type = "b", col = 2)
par(mfrow = c(1, 1))
plot(results$pe_e[1,] * 100, type = "b", xlab = "number of factors",
     ylab = "Prediction error (percentage)",
     ylim = c(0, 100 * max(results$pe_e[1,])))
lines(results$pe_e_linear[1,] * 100, type = "b", col = 2)
lines(results$pe_e_linear_direct[1,] * 100, type = "b", col = 3)
