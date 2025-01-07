library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

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

# No NRMNN data for DE, HI, MD, NE
# 144 Missing data for FL

sectors = c(# "NRMNN", "CONSN", "MFGN", 
               "TRADN", "INFON", "FIREN", "PBSVN", "EDUHN", "LEIHN", "SRVON")

data_mat = array(NA, dim = c(418, length(states), length(sectors)))

for (i in 1:length(sectors)) {
  industry = sectors[i]
  for (j in 1:length(states)) {
    state = states[j]
    temp = read.csv(paste0("./Employee_data/", state, industry, ".csv"), header = T)
    
    data_mat[,j,i] = temp$value
    
    if (sum(is.na(temp$value)) > 0) {
      cat(state, industry, sum(is.na(temp$value)), "\n")
    }
  }
}

res = data_mat
if (by_cat == "states") {
  for (t in 1:418) {
    temp = res[t,,] / rowSums(res[t,,])
    res[t,,] = temp
  }
  x = res
} else if (by_cat == "sectors") {
  for (t in 1:418) {
    temp = t(t(res[t,,]) / colSums(res[t,,]))
    res[t,,] = temp
  }
  x = array(NA, dim = dim(res)[c(1,3,2)])
  for (t in 1:418) {
    x[t,,] = t(res[t,,])
  }
} else {
  stop("by_cat must be either states or sectors")
}

rm(list=setdiff(ls(), c("x", "states", "sectors", "by_cat")))


# Analysis by states: Product of 50 6-dimensional spheres
# res = data_mat
# for (t in 1:418) {
#   temp = res[t,,] / rowSums(res[t,,])
#   res[t,,] = temp
# }

# plots = vector("list", length(states))
# for (j in 1:length(states)) {
#   temp = res[,j,]
#   colnames(temp) = sectors
#   ts_temp = data.frame(
#     Date = seq(as.Date("1990-01-01"), as.Date("2024-10-01"), by = "month"),
#     as.data.frame(temp)
#   )
#   
#   long_data <- ts_temp %>%
#     pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value")
#   
#   # Create the plot
#   if (j == 1) {
#     plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
#       geom_line(size = 0.5) + 
#       scale_x_date(
#         date_breaks = "5 years", 
#         date_labels = "%y"      
#       ) +
#       labs(
#         title = paste0(states[j], " Compositions"),
#         x = "Year",
#         y = "Composition",
#         color = "Sector"
#       ) +
#       theme_minimal()
#   } else {
#     plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
#       geom_line(size = 0.5) + 
#       scale_x_date(
#         date_breaks = "5 years", 
#         date_labels = "%y"      
#       ) +
#       labs(
#         title = paste0(states[j], " Compositions"),
#         x = "Year",
#         y = "Composition",
#         color = "Sector"
#       ) +
#     theme_minimal() +
#     theme(legend.position = "none")        
#   }
#   # plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
#   #   geom_line(size = 0.5) + 
#   #   scale_x_date(
#   #     date_breaks = "5 years", 
#   #     date_labels = "%Y"      
#   #   ) +
#   #   labs(
#   #     title = paste0(states[j], " Compositions"),
#   #     x = "Year",
#   #     y = "Composition",
#   #     color = "Sector"
#   #   ) #+
#     #theme_minimal() #+ 
#     #theme(legend.position = "none") #+        
#     #theme(
#     #  text = element_text(size = 14),  
#     #  axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels for readability
#     #)
# }
# for (i in 1:5) {
#   grid.arrange(grobs = plots[((i - 1) * 10 + 1):(i * 10)], ncol = 2, nrow = 5)
# }


# Analysis by sectors: Product of 7 49-dimensional spheres
# res = data_mat
# for (t in 1:418) {
#   temp = t(t(res[t,,]) / colSums(res[t,,]))
#   res[t,,] = temp
# }

# plots = vector("list", length(sectors))
# for (j in 1:length(sectors)) {
#   temp = res[,,j]
#   colnames(temp) = states
#   ts_temp = data.frame(
#     Date = seq(as.Date("1990-01-01"), as.Date("2024-10-01"), by = "month"),
#     as.data.frame(temp)
#   )
#   
#   long_data <- ts_temp %>%
#     pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value")
#   
#   # Create the plot
#   if (j == 1) {
#     plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
#       geom_line(size = 0.5) + 
#       scale_x_date(
#         date_breaks = "5 years", 
#         date_labels = "%Y"      
#       ) +
#       labs(
#         title = paste0(sectors[j], " Compositions"),
#         x = "Year",
#         y = "Composition",
#         color = "Sector"
#       ) +
#       theme_minimal() 
#   } else {
#     plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
#       geom_line(size = 0.5) + 
#       scale_x_date(
#         date_breaks = "5 years", 
#         date_labels = "%Y"      
#       ) +
#       labs(
#         title = paste0(sectors[j], " Compositions"),
#         x = "Year",
#         y = "Composition",
#         color = "Sector"
#       ) +
#       theme_minimal() + 
#       theme(legend.position = "none")
#   }
# }
# for (i in 1:7) {
#   grid.arrange(grobs = plots[((i - 1) * 1 + 1):(i * 1)], ncol = 1, nrow = 1)
# }
# # Trimmed time plots
# plots = vector("list", length(sectors))
# for (j in 1:length(sectors)) {
#   temp = res[,,j]
#   colnames(temp) = states
#   
#   state_idx = which(states %in% c("CA", "NY", "TX", "FL", "PA", "IL")) #which(temp[1,] > sort(temp[1,], decreasing = T)[6])
#   
#   ts_temp = data.frame(
#     Date = seq(as.Date("1990-01-01"), as.Date("2024-10-01"), by = "month"),
#     as.data.frame(temp[,state_idx])
#   )
#   
#   long_data <- ts_temp %>%
#     pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value")
#   
#   # Create the plot
#   plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
#     geom_line(size = 0.5) + 
#     scale_x_date(
#       date_breaks = "5 years", 
#       date_labels = "%Y"      
#     ) +
#     labs(
#       title = paste0(sectors[j], " Compositions"),
#       x = "Year",
#       y = "Composition",
#       color = ""
#     ) +
#     theme_minimal() 
# }
# for (i in 1:7) {
#   grid.arrange(grobs = plots[((i - 1) * 1 + 1):(i * 1)], ncol = 1, nrow = 1)
# }







