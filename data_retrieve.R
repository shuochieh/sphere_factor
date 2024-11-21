library("httr")
library("jsonlite")


states = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", 
           "IL", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", 
           "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", 
           "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WI", "WY")

industries = c("NRMNN", "CONSN", "MFGN", 
               "TRADN", "INFON", "FIREN", "PBSVN", "EDUHN", "LEIHN", "SRVON")

for (state in states) {
  for (industry in industries) {
    series_id = paste0(state, industry)
    dir = paste0("https://api.stlouisfed.org/fred/series/observations?series_id=", series_id, "&")
    dir = paste0(dir, "api_key=", api_key, "&file_type=json")
    x = GET(dir)
    
    data <- content(x, as = "text", encoding = "UTF-8")
    parsed_data <- fromJSON(data, flatten = TRUE)  # Flatten for easy data frame conversion
    
    
  }
}

series_id = paste0(state, "NRMNN")
dir = paste0("https://api.stlouisfed.org/fred/series/observations?series_id=", series_id, "&")
dir = paste0(dir, "api_key=", api_key, "&file_type=json")
x = GET(dir)

data <- content(x, as = "text", encoding = "UTF-8")
parsed_data <- fromJSON(data, flatten = TRUE)  # Flatten for easy data frame conversion
print(parsed_data)









x = GET("https://api.stlouisfed.org/fred/series/observations?series_id=KSNRMNN&api_key=57f6a31436a3039bdd7e1b23106506d9&file_type=json")

# Step 1: Define API details
base_url <- "https://api.example.com/data"
query_params <- list(
  type = "income",
  year = "2023"
)
headers <- add_headers(
  Authorization = paste("Bearer", "YOUR_API_KEY")
)

# Step 2: Make the request
response <- GET(base_url, headers, query = query_params)

# Step 3: Check the response
if (status_code(response) == 200) {
  # Step 4: Parse the JSON data
  data <- content(x, as = "text", encoding = "UTF-8")
  parsed_data <- fromJSON(data, flatten = TRUE)  # Flatten for easy data frame conversion
  print(parsed_data)
  
  # Step 5: Save the data
  write.csv(parsed_data, "income_data.csv", row.names = FALSE)
  print("Data saved to 'income_data.csv'")
} else {
  # Handle errors
  print(paste("Error:", status_code(response), content(response, "text")))
}