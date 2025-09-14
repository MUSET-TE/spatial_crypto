library(crypto2)
library(dplyr)
library(purrr)
library(writexl)

# Hitorical listing of the top 100 cryptos
top100 <- crypto_listings(which = "latest", sort = "cmc_rank") %>%
  filter(cmc_rank <= 100) %>%
  arrange(cmc_rank)

# In a next step we download time series data for these coins.

i = 1

for (i in 1:nrow(top100)) {
  
  symbol_ <- top100$symbol[i]
  message(symbol_)
  
  quot_ <- crypto_history(coins[i,], 
                          start_date = "20130501",
                          end_date = "20250901")
  
  if(!exists("quotations_all")) quotations_all <- quot_ else 
    quotations_all <- rbind(quotations_all, quot_)
  
  rm(symbol_, quot_)
  
}
write_xlsx(quotations_all, "top100_20130501_20250901.xlsx")