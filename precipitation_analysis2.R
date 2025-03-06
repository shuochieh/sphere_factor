source("./main_func.R")

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

station_data_get = function (ss_ls) {
  p = length(ss_ls)
  for (i in 1:p) {
    state_name = names(ss_ls)[i]
    for (j in 1:ss_ls[[i]]) {
      dir_name = paste0("./precipitation/", state_name, j, ".csv")
      temp = read_precip_data(dir_name)
      if (i == 1 && j == 1) {
        n = nrow(temp)
        Dates = convert_to_date(temp[,1])
        res = matrix(0, nrow = n, ncol = ss_ls[[1]])
      } 
      res[,j] = temp[,2]
    }
  }
  res[which(rowSums(res) == 0),] = 1 / ss_ls[[1]]
  res = res / rowSums(res)
  
  return (list("dta" = res, "Dates" = Dates))
}

x = list("CA" = station_data_get(list("CA" = 7))$dta)
# x[["OR"]] = station_data_get(list("OR" = 9))$dta
# x[["WA"]] = station_data_get(list("WA" = 10))$dta

x[["AZ"]] = station_data_get(list("AZ" = 7))$dta
# x = list("AZ" = station_data_get(list("AZ" = 7))$dta)
# x[["CO"]] = station_data_get(list("CO" = 5))$dta
# x[["ID"]] = station_data_get(list("ID" = 10))$dta
# x[["MT"]] = station_data_get(list("MT" = 7))$dta
x[["NV"]] = station_data_get(list("NV" = 4))$dta
# x[["NM"]] = station_data_get(list("NM" = 8))$dta
# x[["UT"]] = station_data_get(list("UT" = 7))$dta
# x[["WY"]] = station_data_get(list("WY" = 10))$dta

# x = list("FL" = station_data_get(list("FL" = 7))$dta)
# x[["GA"]] = station_data_get(list("GA" = 9))$dta
# x[["NC"]] = station_data_get(list("NC" = 8))$dta
# x[["SC"]] = station_data_get(list("SC" = 7))$dta
# x[["VA"]] = station_data_get(list("VA" = 6))$dta
# x[["MD"]] = station_data_get(list("MD" = 8))$dta
# x[["WV"]] = station_data_get(list("WV" = 6))$dta

# x = list("AL" = station_data_get(list("AL" = 8))$dta)
# x[["KY"]] = station_data_get(list("KY" = 4))$dta
# x[["MS"]] = station_data_get(list("MS" = 10))$dta
# x[["TN"]] = station_data_get(list("TN" = 4))$dta

# x = list("AR" = station_data_get(list("AR" = 9))$dta)
# x[["LA"]] = station_data_get(list("LA" = 9))$dta
# x[["OK"]] = station_data_get(list("OK" = 9))$dta
# x[["TX"]] = station_data_get(list("TX" = 10))$dta

# x = list("IL" = station_data_get(list("IL" = 9))$dta)
# x[["IN"]] = station_data_get(list("IN" = 9))$dta
# x[["MI"]] = station_data_get(list("MI" = 10))$dta
# x[["OH"]] = station_data_get(list("OH" = 10))$dta
# x[["WI"]] = station_data_get(list("WI" = 9))$dta

# x = list("IA" = station_data_get(list("IA" = 9))$dta)
# x[["KS"]] = station_data_get(list("KS" = 9))$dta
# x[["MN"]] = station_data_get(list("MN" = 9))$dta
# x[["MO"]] = station_data_get(list("MO" = 6))$dta
# x[["NE"]] = station_data_get(list("NE" = 8))$dta
# x[["ND"]] = station_data_get(list("ND" = 9))$dta
# x[["SD"]] = station_data_get(list("SD" = 9))$dta

# x = list("CT" = station_data_get(list("CT" = 3))$dta)
# x[["ME"]] = station_data_get(list("ME" = 3))$dta
# x[["MA"]] = station_data_get(list("MA" = 3))$dta
# x[["VT"]] = station_data_get(list("VT" = 3))$dta
# x[["NY"]] = station_data_get(list("NY" = 10))$dta
# x[["NJ"]] = station_data_get(list("NJ" = 3))$dta
# x[["PA"]] = station_data_get(list("PA" = 10))$dta


dates = station_data_get(list("CA" = 7))$Date

for (state in names(x)) {
  cat(state, "total variation:", round(100 * sqrt(sum(apply(x[[state]], 2, var))), 4), "% \n")
  x[[state]] = sqrt(x[[state]])
}

par(mfrow = c(3, 4))
for (state in names(x)) {
  plot(x = as.Date(dates),
       y = geod_sphere(x[[state]], mean_on_sphere(x[[state]])),
       type = "l",
       xlab = "",
       ylab = "Geo Dist to mean",
       main = state)
}

results = main_uneven_sphere(x, r = 15, h = 25, test_size = 360)

par(mfrow = c(1, 1))
plot(results$FVU_e, type = "b", xlim = c(0, 35), ylim = c(0,1),
     xlab = "number of factors",
     ylab = "")
lines(results$FVU_e_linear, type = "b", col = 2)
# lines(x = c(1:length(results$FVU_e_sprt)) * length(x), y = results$FVU_e_sprt, type = "b", col = 4)
abline(h = 0.1)

par(mfrow = c(1, 3))
for (i in 1:length(names(x))) {
  plot(results$pe_e[i,] * 100, type = "b", xlab = "number of factors",
       ylab = "percentage", main = names(x)[i],
       ylim = c(0, 100 * max(results$pe_e[i,])))
  lines(results$pe_e_linear[i,] * 100, type = "b", col = 2)
  lines(results$pe_e_linear_direct[i,] * 100, type = "b", col = 3)
}

