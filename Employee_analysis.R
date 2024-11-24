### Analysis of employment data

# Need res from previous script

states = c("AL", "AK", "AZ", "AR", 
           "CA", "CO", "CT", 
           "DE", 
           "FL", 
           "GA", 
           "HI", 
           "ID", "IL", "IA", "IN", 
           "KS", "KY", 
           "LA", 
           "ME", "MA", "MI", "MN", "MS", "MO", "MT", 
           "MD",
           "NV", "NH", "NJ", "NM", "NY", "NC", "ND", 
           "NE",
           "OH", "OK", "OR", 
           "PA", "RI", 
           "SC", "SD", 
           "TN", "TX", 
           "UT", 
           "VT", "VA", 
           "WA", "WI", "WV", "WY")

industries = c("TRADN", "INFON", "FIREN", "PBSVN", "EDUHN", "LEIHN", "SRVON")


NewEng_idx = which(states %in% c("ME", "NH", "VT", "MA", "RI", "CT"))
MidAtl_idx = which(states %in% c("NY", "NJ", "PA"))
Northeast = c(NewEng_idx, MidAtl_idx)

ENCen_idx = which(states %in% c("OH", "IN", "IL", "MI", "WI"))
WNCen_idx = which(states %in% c("MN", "IA", "MO", "ND", "SD", "NE", "KS"))
Midwest = c(ENCen_idx, WNCen_idx)

SouAtl_idx = which(states %in% c("DE", "MD", "VA", "WV", "NC", "SC", "GA", "FL"))
ESCen_idx = which(states %in% c("KY", "TN", "AL", "MS"))
WSCen_idx = which(states %in% c("AR", "LA", "OK", "TX"))
South = c(SouAtl_idx, ESCen_idx, WSCen_idx)

Montn_idx = which(states %in% c("MT", "ID", "WY", "NV", "UT", "CO", "AZ", "NM"))
PCF_idx = which(states %in% c("WA", "OR", "CA", "AK", "HI"))
West = c(Montn_idx, PCF_idx)

x = array(NA, dim = c(418, 7, 50))

for (t in 1:418) {
  x[t,,] = t(res[t,,])
}

x = sqrt(x) # square root transform

# Estimate mu
mu_hat = matrix(NA, nrow = 7, ncol = 50)
trans_x = matrix(NA, nrow = 418, ncol = 7 * 50)
for (j in 1:7) {
  mu_hat[j,] = mean_on_sphere(x[,j,])
  trans_x[,((j - 1) * 50 + 1):(j * 50)] = Log_sphere(x[,j,], mu_hat[j,])
}

colnames(mu_hat) = states
for (j in 1:length(industries)) {
  barplot(mu_hat[j,], cex.names = 0.9, las = 3, ylim = c(0, ifelse(j == 2, 0.41, 0.35)), 
          main = paste("Estimated mean composition of", industries[j]))
}

r = 4
model = LYB_fm(trans_x, r = r, h = 5)
V = model$V
Fac = model$f_hat

Northeast_V = matrix(NA, nrow = length(Northeast) * 7, ncol = r)
Midwest_V = matrix(NA, nrow = length(Midwest) * 7, ncol = r)
South_V = matrix(NA, nrow = length(South) * 7, ncol = r)
West_V = matrix(NA, nrow = length(West) * 7, ncol = r)
Vs = vector("list", 7)
for (j in 1:7) {
  Vs[[j]] = V[((j - 1) * 50 + 1):(j * 50),]
  Northeast_V[((j - 1) * length(Northeast) + 1):(j * length(Northeast)),] = Vs[[j]][Northeast,]
  Midwest_V[((j - 1) * length(Midwest) + 1):(j * length(Midwest)),] = Vs[[j]][Midwest,]
  South_V[((j - 1) * length(South) + 1):(j * length(South)),] = Vs[[j]][South,]
  West_V[((j - 1) * length(West) + 1):(j * length(West)),] = Vs[[j]][West,]
}

par(mfrow = c(2, 2))
for (i in 1:7) {
  for (j in 1:1) {
    plot(x = Vs[[i]][,j], y = Vs[[i]][,j+1], pch = 20, 
         xlab = paste("Loading for factor", j),
         ylab = paste("Loading for factor", j + 1),
         main = paste("Estimated factor loading for", industries[i]))
    
    points(x = Vs[[i]][5,j], y = Vs[[i]][5,j+1], pch = 19, col = 2, cex = 2)
    text(x = Vs[[i]][5,j], y = Vs[[i]][5,j+1], states[5], col = 2, pos = 2)
    
    points(x = Vs[[i]][31,j], y = Vs[[i]][31,j+1], pch = 19, col = 3,
           xlim = c(-1.5, 1), ylim = c(-1.5, 1), cex = 2)
    text(x = Vs[[i]][31,j], y = Vs[[i]][31,j+1], states[31], col = 3, pos = 2)
    
    points(x = Vs[[i]][43,j], y = Vs[[i]][43,j+1], pch = 19, col = 4, cex = 2)
    text(x = Vs[[i]][43,j], y = Vs[[i]][43,j+1], states[43], col = 4, pos = 2)
    
    points(x = Vs[[i]][9,j], y = Vs[[i]][9,j+1], pch = 19, col = 5, cex = 2)
    text(x = Vs[[i]][9,j], y = Vs[[i]][9,j+1], states[9], col = 5, pos = 2)
    
    points(x = Vs[[i]][38,j], y = Vs[[i]][38,j+1], pch = 19, col = 6, cex = 2)
    text(x = Vs[[i]][38,j], y = Vs[[i]][38,j+1], states[38], col = 6, pos = 2)
    
    points(x = Vs[[i]][13,j], y = Vs[[i]][13,j+1], pch = 19, col = 7, cex = 2)
    text(x = Vs[[i]][13,j], y = Vs[[i]][13,j+1], states[13], col = 7, pos = 2)
  }
}

par(mfrow = c(2, 2))
plot(x = seq(from = as.Date("1990-01-01"), to = as.Date("2024-10-01"), by = "month"),
     y = Fac[,1], type = "l",
     xlab = "")
plot(x = seq(from = as.Date("1990-01-01"), to = as.Date("2024-10-01"), by = "month"),
     y = Fac[,2], type = "l",
     xlab = "")
plot(x = seq(from = as.Date("1990-01-01"), to = as.Date("2024-10-01"), by = "month"),
     y = Fac[,3], type = "l",
     xlab = "")
plot(x = seq(from = as.Date("1990-01-01"), to = as.Date("2024-10-01"), by = "month"),
     y = Fac[,4], type = "l",
     xlab = "")

par(mfrow = c(1,1))
Dates = seq(from = as.Date("1990-01-01"), to = as.Date("2024-10-01"), by = "month")
for (i in 1:7) {
  par(mfrow = c(3, 1))
  for (j in 1:3) {
    y_max = max(abs(Vs[[i]][,j] %o% Fac[,j]))
    if (j == 1) {
      temp = "First"
    } else if (j == 2) {
      temp = "Second"
    } else if (j == 3) {
      temp = "Third"
    } else {
      stop("Too many factors")
    }
    plot(x = Dates, y = Vs[[i]][5,j] * Fac[,j], type = "l", lwd = 1.5, 
         col = 2,
         xlab = "Time",
         ylab = "",
         main = paste(temp, "factor-driven dynamic for", industries[i]),
         ylim = c(-y_max, y_max))
    text(x = tail(Dates, 1), y = tail(Vs[[i]][5,j] * Fac[,j], 1), 
         states[5], col = 2, pos = 1)
    
    lines(x = Dates, y = Vs[[i]][31,j] * Fac[,j], col = 3, lwd = 1.5)
    text(x = tail(Dates, 1), y = tail(Vs[[i]][31,j] * Fac[,j], 1), 
         states[31], col = 3, pos = 1)

    lines(x = Dates, y = Vs[[i]][43,j] * Fac[,j], col = 4, lwd = 1.5)
    text(x = tail(Dates, 1), y = tail(Vs[[i]][43,j] * Fac[,j], 1), 
         states[43], col = 4, pos = 1)

    lines(x = Dates, y = Vs[[i]][9,j] * Fac[,j], col = 5, lwd = 1.5)
    text(x = tail(Dates, 1), y = tail(Vs[[i]][9,j] * Fac[,j], 1), 
         states[9], col = 5, pos = 1)

    lines(x = Dates, y = Vs[[i]][38,j] * Fac[,j], col = 6, lwd = 1.5)
    text(x = tail(Dates, 1), y = tail(Vs[[i]][38,j] * Fac[,j], 1), 
         states[38], col = 6, pos = 1)

    lines(x = Dates, y = Vs[[i]][13,j] * Fac[,j], col = 7, lwd = 1.5)
    text(x = tail(Dates, 1), y = tail(Vs[[i]][13,j] * Fac[,j], 1), 
         states[13], col = 7, pos = 1)
  }
}










