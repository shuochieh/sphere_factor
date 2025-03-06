library(dplyr)
library(lubridate)

source("./main_func.R")

read_LT <- function(file) {
  df <- read.csv(file, skip = 2, header = T, sep = "")
  return(df[,c("Year", "Age", "dx")])
}
countries = c("USA", 
              "CAN", # Canada
              "CHE", # Switzerland
              "DNK", # Denmark
              "ESP", # Spain
              "FIN", # Finland
              "FRATNP", # France Total
              "GBR_NP", # United Kingdom
              "NOR" # Norway
              )


male_x = array(NA, dim = c(90, length(countries), 40))
female_x = array(NA, dim = c(90, length(countries), 40))

for (i in 1:length(countries)) {
  country = countries[i]
  dir_name = paste0("./life_tables/", country, ".mltper_1x1.txt")
  temp = read_LT(dir_name)
  idx = which(is.na(as.numeric(temp$Age)))
  temp = temp[-idx,]
  
  temp2 = temp %>%
    mutate(Age = as.numeric(Age)) %>%
    group_by(Age) %>%
    arrange(Age) %>%
    ungroup()
  
  temp2 = temp2[which(temp2$Year > 1932),]
  temp2 = temp2[which(temp2$Year < 2023),]
  
  for (j in 61:100) {
    male_x[,i,(j - 60)] = c(unlist(temp2[which(temp2$Age == j),"dx"]))
  }
  
  country = countries[i]
  dir_name = paste0("./life_tables/", country, ".fltper_1x1.txt")
  temp = read_LT(dir_name)
  idx = which(is.na(as.numeric(temp$Age)))
  temp = temp[-idx,]
  
  temp2 = temp %>%
    mutate(Age = as.numeric(Age)) %>%
    group_by(Age) %>%
    arrange(Age) %>%
    ungroup()
  
  temp2 = temp2[which(temp2$Year > 1932),]
  temp2 = temp2[which(temp2$Year < 2023),]
  
  for (j in 60:100) {
    female_x[,i,(j - 60)] = c(unlist(temp2[which(temp2$Age == j),"dx"]))
  }
}
        
male_x_temp = array(NA, dim = c(90, length(countries), 8))
female_x_temp = array(NA, dim = c(90, length(countries), 8))

for (i in 1:length(countries)) {
  for (j in 1:8) {
    male_x_temp[,i,j] = rowSums(male_x[,i,((j - 1) * 5 + 1):(j * 5)])
    female_x_temp[,i,j] = rowSums(female_x[,i,((j - 1) * 5 + 1):(j * 5)])
  }
  
  male_x_temp[,i,] = male_x_temp[,i,] / rowSums(male_x_temp[,i,])
  female_x_temp[,i,] = female_x_temp[,i,] / rowSums(female_x_temp[,i,])
}

male_x = vector("list", length = length(countries))
female_x = vector("list", length = length(countries))

for (i in 1:length(countries)) {
  male_x[[i]] = sqrt(male_x_temp[,i,])
  female_x[[i]] = sqrt(female_x_temp[,i,])
}

dates = seq(from = as.Date("1933-01-01"), to = as.Date("2022-01-01"), by = "year")
par(mfrow = c(3, 10))
for (i in 1:length(countries)) {
  country = countries[i]
  for (j in 1:10) {
    barplot(male_x_temp[(j - 1) * 9 + 3,i,], main = paste(country, year(dates[(j - 1) * 9 + 3])))
  }
}
for (i in 1:length(countries)) {
  country = countries[i]
  for (j in 1:10) {
    barplot(female_x_temp[(j - 1) * 9 + 3,i,], main = paste(country, year(dates[(j - 1) * 9 + 3])))
  }
}

par(mfrow = c(3, 3))
for (i in 1:length(countries)) {
  temp = geod_sphere(male_x[[i]], mean_on_sphere(male_x[[i]]))
  plot(x = year(dates),
       y = temp, type = "l", main = countries[i],
       xlab = "", ylab = "Geo Dist to Mean")
}

for (i in 1:length(countries)) {
  temp = geod_sphere(female_x[[i]], mean_on_sphere(female_x[[i]]))
  plot(x = year(dates),
       y = temp, type = "l", main = countries[i],
       xlab = "", ylab = "Geo Dist to Mean")
}

results = main_uneven_sphere(male_x, r = 10, test_size = 20, h = 3)

par(mfrow = c(1, 1))
plot(results$FVU_e, type = "b", xlab = "number of factors", ylab = "", ylim = c(0, 0.2))
lines(results$FVU_e_linear, type = "b", col = 2)

par(mfrow = c(3, 3))
for (i in 1:length(countries)) {
  plot(results$pe_e[i,] * 100, type = "b", 
       main = countries[i],
       ylim = c(0, max(c(results$pe_e, results$pe_e_linear) * 100)),
       xlab = "number of factors",
       ylab = "Prediction error (percentage)")
  lines(results$pe_e_linear[i,] * 100, type = "b",
        col = 2)
  # lines(results$pe_e_linear_direct[i,] * 100, type = "b",
  #       col = 3)
}


results = main_uneven_sphere(female_x, r = 10, test_size = 20, h = 3)

par(mfrow = c(1, 1))
plot(results$FVU_e, type = "b", xlab = "number of factors", ylab = "", ylim = c(0, 0.2))
lines(results$FVU_e_linear, type = "b", col = 2)

par(mfrow = c(3, 3))
for (i in 1:length(countries)) {
  plot(results$pe_e[i,] * 100, type = "b", 
       main = countries[i],
       ylim = c(0, max(c(results$pe_e, results$pe_e_linear) * 100)),
       xlab = "number of factors",
       ylab = "Prediction error (percentage)")
  lines(results$pe_e_linear[i,] * 100, type = "b",
        col = 2)
  # lines(results$pe_e_linear_direct[i,] * 100, type = "b",
  #       col = 3)
}


