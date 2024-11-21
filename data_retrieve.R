library("httr")
library("jsonlite")

states = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", 
           "IL", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", 
           "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", 
           "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WI", "WY")

industries = c("NRMNN", "CONSN", "MFGN", 
               "TRADN", "INFON", "FIREN", "PBSVN", "EDUHN", "LEIHN", "SRVON")

Sys.sleep(30)
for (industry in industries) {
  for (state in states) {
    series_id = paste0(state, industry)
    dir = paste0("https://api.stlouisfed.org/fred/series/observations?series_id=", series_id, "&")
    dir = paste0(dir, "api_key=", api_key, "&file_type=json")
    x = GET(dir)
    
    if (status_code(x) == 200) {
      data <- content(x, as = "text", encoding = "UTF-8")
      parsed_data <- fromJSON(data, flatten = TRUE)  # Flatten for easy data frame conversion
      
      res = parsed_data$observations[,c("date", "value")]
      res[,c("date")] = as.Date(res[,c("date")])
      res[,c("value")] = as.numeric(res[,c("value")])
      
      write.csv(res, paste0("./Employee_data/", state, industry, ".csv"), row.names = FALSE)
    } else {
      cat(c(state, industry), "\n")
      cat(paste("Error:", status_code(x), content(x, "text")), "\n")
    }
  }
  
  Sys.sleep(65)
}

