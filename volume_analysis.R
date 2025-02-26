source("main_func.r")

vol_dol = read.csv("./sp500_covariance/volume_dollar_data")
dates = as.Date(vol_dol[,1])

# Categorize
temp = read.csv("./sp500_covariance/sp500.csv")
temp = temp[c("Symbol", "GICS.Sector", "GICS.Sub.Industry")]

sectors = unlist(unique(temp["GICS.Sector"]))
sub_industries = vector("list", length = length(sectors))
for (i in 1:length(sectors)) {
  sub_industries[[i]] = unique(temp[which(temp[,"GICS.Sector"] == sectors[i]),"GICS.Sub.Industry"])
}

# treat missing as zero
vol_dol[,-1][is.na(vol_dol[,-1])] = 0

par(mfrow = c(1 , 1))
x = vector("list", length = 11)
for (zz in 1:11){
  sub_ind_composition = matrix(NA, nrow = nrow(vol_dol), ncol = length(sub_industries[[zz]]))
  for (i in 1:length(sub_industries[[zz]])) {
    sub_industry = sub_industries[[zz]][i]
    companies = temp[which(temp[,"GICS.Sub.Industry"] == sub_industry), "Symbol"]
    
    sub_matrix = vol_dol[,which(colnames(vol_dol) %in% companies)]
    
    if (!is.data.frame(sub_matrix)) {
      sub_ind_composition[,i] = sub_matrix
    } else {
      sub_ind_composition[,i] = rowSums(sub_matrix)
    }
  }
  
  sub_ind_composition = sub_ind_composition / rowSums(sub_ind_composition)
  

  temp_x = sqrt(sub_ind_composition)
  
  mu_hat = mean_on_sphere(temp_x, verbose = T)
  plot(x=as.Date(dates), y=geod_sphere(temp_x, mu_hat), type="l",
       main = sectors[zz],
       xlab = "",
       ylab = "Geo. Dist. to Frechet Mean")
  for (j in 1:ncol(sub_ind_composition)) {
    if (j == 1) {
      pal = rainbow(ncol(sub_ind_composition))
      plot(x = dates, y = sub_ind_composition[,j], type = "l", ylim = range(sub_ind_composition),
           xlab = "",
           ylab = "composition",
           main = sectors[zz],
           col = pal[1])
    } else {
      lines(x = dates, y = sub_ind_composition[,j], type = "l", col = pal[j])
    }
    legend("topleft",
           lty = 1,
           col = pal,
           legend = sub_industries[[zz]],
           cex = 0.8)
  }
  
  x[[zz]] = temp_x
}

# Sector-specific
results = main(array(x[[6]], c(nrow(x[[6]]), 1, ncol(x[[6]]))), 15, test_size = 1004)
plot(results$FVU_e, type = "b", 
     xlab = "Number of factors",
     ylab = "Fraction of (Euclidean) variations unexplained")
lines(results$FVU_e_linear, col = 2, type = "b")
legend("topright",
       col = c(1,2),
       lty = c(1, 1),
       pch = c(16, 16),
       legend = c("sphere factor", "linear factor"))
plot(results$pe_e, type = "b",
     xlab = "Number of factors",
     ylab = "(Euclidean) Prediction errors", ylim = c(0, 0.25))
lines(results$pe_e_linear, col = 2, type = "b")
legend("topright",
       col = c(1,2),
       lty = c(1, 1),
       pch = c(16, 16),
       legend = c("sphere factor", "linear factor"))

par(mfrow = c(1, 1))
plot(results$FVU_e, type = "b")
lines(results$FVU_e_linear, type = "b", col = 2)
plot(results$pe_e, type = "b", ylim = c(0,max(results$pe_e)))
lines(results$pe_e_linear, type = "b", col = 2)

# Multiple sectors
results = main_uneven_sphere(x[c(3,6)], 15, test_size = 1004)
plot(results$FVU_e, type = "b")
lines(results$FVU_e_linear, col = 2, type = "b")
plot(results$pe_e, type = "b")
lines(results$pe_e_linear, col = 2, type = "b")

par(mfrow = c(1, 1))
plot(results$FVU_e, type = "b")
lines(results$FVU_e_linear, type = "b", col = 2)
plot(results$pe_e, type = "b")
lines(results$pe_e_linear, type = "b", col = 2)


















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

results = main(x, 9, test_size = 752)

plot(results$FVU_e, type="b")
lines(x=1:9, y=results$FVU_e_linear, type="b", col=2)
plot(results$pe_e,type="b")
lines(x=1:9, y=results$pe_e_linear, type="b", col=2)











