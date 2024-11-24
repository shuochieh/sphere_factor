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

industries = c(# "NRMNN", "CONSN", "MFGN", 
               "TRADN", "INFON", "FIREN", "PBSVN", "EDUHN", "LEIHN", "SRVON")

data_mat = array(NA, dim = c(418, length(states), length(industries)))

for (i in 1:length(industries)) {
  industry = industries[i]
  for (j in 1:length(states)) {
    state = states[j]
    temp = read.csv(paste0("./Employee_data/", state, industry, ".csv"), header = T)
    
    data_mat[,j,i] = temp$value
    
    if (sum(is.na(temp$value)) > 0) {
      cat(state, industry, sum(is.na(temp$value)), "\n")
    }
  }
}

# Prelim analysis: by states
res = data_mat
for (t in 1:418) {
  temp = res[t,,] / rowSums(res[t,,])
  res[t,,] = temp
}

plots = vector("list", length(states))
for (j in 1:length(states)) {
  temp = res[,j,]
  colnames(temp) = industries
  ts_temp = data.frame(
    Date = seq(as.Date("1990-01-01"), as.Date("2024-10-01"), by = "month"),
    as.data.frame(temp)
  )
  
  long_data <- ts_temp %>%
    pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value")
  
  # Create the plot
  if (j == 1) {
    plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
      geom_line(size = 0.5) + 
      scale_x_date(
        date_breaks = "5 years", 
        date_labels = "%y"      
      ) +
      labs(
        title = paste0(states[j], " Compositions"),
        x = "Year",
        y = "Composition",
        color = "Sector"
      ) +
      theme_minimal()
  } else {
    plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
      geom_line(size = 0.5) + 
      scale_x_date(
        date_breaks = "5 years", 
        date_labels = "%y"      
      ) +
      labs(
        title = paste0(states[j], " Compositions"),
        x = "Year",
        y = "Composition",
        color = "Sector"
      ) +
    theme_minimal() +
    theme(legend.position = "none")        
  }
  # plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
  #   geom_line(size = 0.5) + 
  #   scale_x_date(
  #     date_breaks = "5 years", 
  #     date_labels = "%Y"      
  #   ) +
  #   labs(
  #     title = paste0(states[j], " Compositions"),
  #     x = "Year",
  #     y = "Composition",
  #     color = "Sector"
  #   ) #+
    #theme_minimal() #+ 
    #theme(legend.position = "none") #+        
    #theme(
    #  text = element_text(size = 14),  
    #  axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels for readability
    #)
}

for (i in 1:5) {
  grid.arrange(grobs = plots[((i - 1) * 10 + 1):(i * 10)], ncol = 2, nrow = 5)
}


# Prelim analysis: by sectors
res = data_mat
for (t in 1:418) {
  temp = t(t(res[t,,]) / colSums(res[t,,]))
  res[t,,] = temp
}

plots = vector("list", length(industries))
for (j in 1:length(industries)) {
  temp = res[,,j]
  colnames(temp) = states
  ts_temp = data.frame(
    Date = seq(as.Date("1990-01-01"), as.Date("2024-10-01"), by = "month"),
    as.data.frame(temp)
  )
  
  long_data <- ts_temp %>%
    pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value")
  
  # Create the plot
  if (j == 1) {
    plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
      geom_line(size = 0.5) + 
      scale_x_date(
        date_breaks = "5 years", 
        date_labels = "%Y"      
      ) +
      labs(
        title = paste0(industries[j], " Compositions"),
        x = "Year",
        y = "Composition",
        color = "Sector"
      ) +
      theme_minimal() 
  } else {
    plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
      geom_line(size = 0.5) + 
      scale_x_date(
        date_breaks = "5 years", 
        date_labels = "%Y"      
      ) +
      labs(
        title = paste0(industries[j], " Compositions"),
        x = "Year",
        y = "Composition",
        color = "Sector"
      ) +
      theme_minimal() + 
      theme(legend.position = "none")
  }
}

for (i in 1:7) {
  grid.arrange(grobs = plots[((i - 1) * 1 + 1):(i * 1)], ncol = 1, nrow = 1)
}

# Trimmed time plots

plots = vector("list", length(industries))
for (j in 1:length(industries)) {
  temp = res[,,j]
  colnames(temp) = states
  
  state_idx = which(states %in% c("CA", "NY", "TX", "FL", "PA", "IL")) #which(temp[1,] > sort(temp[1,], decreasing = T)[6])
  
  ts_temp = data.frame(
    Date = seq(as.Date("1990-01-01"), as.Date("2024-10-01"), by = "month"),
    as.data.frame(temp[,state_idx])
  )
  
  long_data <- ts_temp %>%
    pivot_longer(cols = -Date, names_to = "Variable", values_to = "Value")
  
  # Create the plot
  plots[[j]] = ggplot(long_data, aes(x = Date, y = Value, color = Variable)) +
    geom_line(size = 0.5) + 
    scale_x_date(
      date_breaks = "5 years", 
      date_labels = "%Y"      
    ) +
    labs(
      title = paste0(industries[j], " Compositions"),
      x = "Year",
      y = "Composition",
      color = ""
    ) +
    theme_minimal() 
}

for (i in 1:7) {
  grid.arrange(grobs = plots[((i - 1) * 1 + 1):(i * 1)], ncol = 1, nrow = 1)
}







