source("main_func.r")

vol_dol = read.csv("./sp500_covariance/volume_dollar_data")

# Categorize
temp = read.csv("./sp500_covariance/sp500.csv")[c("Symbol", "GICS.Sector")]
sectors = unique(temp[,2])
sectors = sectors[c(1,2,3,6,8,11)]

for (sector in sectors) {
  cat(sector, ":", length(which(temp[,2] == sector)), "\n")
}

# treat missing as zero
vol_dol[,which(colSums(is.na(vol_dol)) == 0)]
vol_dol[,-1][is.na(vol_dol[,-1])] = 0

sector_comp = matrix(NA, ncol = length(sectors), nrow = nrow(vol_dol))
for (i in 1:length(sectors)) {
  sector = sectors[i]
  companies = temp[which(temp[,2] %in% sector),1]
  sub_m = vol_dol[,which(colnames(vol_dol) %in% companies)]
  sector_comp[,i] = rowSums(sub_m)
}
colnames(sector_comp) = c(sectors)
dates = vol_dol[,1]

sector_comp = sector_comp / rowSums(sector_comp) # turn into composition

par(mfrow = c(2, 4))
for (i in 1:length(sectors)) {
  plot(x = as.Date(dates), y = sector_comp[,i] * 100, type = "l",
       xlab = "Time", ylab = "percentage",
       main = sectors[i])
}

sector_comp = sqrt(sector_comp)

mu_hat = mean_on_sphere(sector_comp, verbose = T)

par(mfrow = c(1, 1))
plot(x=as.Date(dates), y=geod_sphere(sector_comp, mu_hat), type="l")

x = array(sector_comp, dim = c(nrow(sector_comp), 1, ncol(sector_comp)))

results = main(x, 3, test_size = 752)

results$FVU_e
results$FVU_e_linear
results$pe_e
results$pe_e_linear











