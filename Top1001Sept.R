library(crypto2)
library(dplyr)
library(purrr)

# preia listarea istorică
coins <- crypto_list(only_active = TRUE)
top100 <- crypto_listings(
  which      = "latest",
  limit      = 500,        # API poate întoarce până la 500
  sort       = "cmc_rank"  # ordonare după rank oficial
)

# selectăm primele 100
top100 <- top100 %>%
  arrange(cmc_rank) %>%
  slice(1:100)

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
write_xlsx(quotations_all, "toatecoinurile.xlsx")
#save(list = "quotations_all",
 #    file = "quotations_top100.RData")

#quotations_top100 <- quotations_all %>%
 # filter(symbol %in% top100$symbol) %>%
  #arrange(symbol)

#quotations_top100 <- quotations_top100 %>%
 # left_join(
  #  top100 %>% select(symbol, cmc_rank),
   # by = "symbol"
  #)

#quotations_top100 <- quotations_top100 %>%
 # arrange(cmc_rank, time_open)  # sau date, depinde ce coloana ai pentru timp

# verificare
#unique(quotations_top100$symbol)  # ar trebui să fie 100

####inca odata\This one stops at ethereum idk why

#quotations_all <- NULL

#for(i in 1:nrow(top100)) {
 # symbol_ <- top100$symbol[i]
  #message("Downloading: ", symbol_)
  
# quot_ <- tryCatch(
 #   crypto_history(top100[i,], 
  #                 start_date = top100$start_date[i], 
   #              end_date = "20250901"),
    #error = function(e) {
    #  message("No data for ", symbol_, ", skipping.")
     # return(NULL)
    }
    
  )
  
#  if(!is.null(quot_)) {
 #   if(is.null(quotations_all)) quotations_all <- quot_ else quotations_all <- rbind(quotations_all, quot_)
  #}
  
#rm(symbol_, quot_)
}
#write_xlsx(top100, "top100.xlsx")
#write_xlsx(quotations_all, "toatecoinurile.xlsx")

#coins <- crypto_list(only_active = FALSE)

# Add id + date_added to your top100
#top100 <- top100 %>%
 # left_join(coins %>% select(id, symbol, name, first_historical_data),
  #          by = "id.x")

# check matches
#sum(!is.na(top100$id))   # should be close to 100



#top100history = top100 %>%
  #left_join(coins, by= "symbol")

#top100_with_first_history_clean <- top100history %>%
  #distinct(symbol, .keep_all = TRUE)


#top100_selected <- top100_with_first_history_clean %>%
  #select(id.x, name.x, symbol, cmc_rank, first_historical_data)

#top100_selected <- top100_selected %>%
 # mutate(first_historical_data = as.Date(first_historical_data))

#quotations_secondtry <- NULL

#for (i in 1:nrow(top100_selected)) {
  
 # id_     <- top100_selected$id.x[i]
  #symbol_ <- top100_selected$symbol[i]
  #start_  <- format(top100_selected$first_historical_data[i], "%Y%m%d")  # YYYYMMDD
  #end_    <- "20250901"
  
  #message("Downloading: ", symbol_, " (id=", id_, ") from ", start_)
  
  # try to download history
 # quot_ <- tryCatch(
  #  crypto_history(coin_list = data.frame(id = id_, symbol = symbol_),
   #                start_date = start_,
    #               end_date   = end_),
    #error = function(e) {
     # message("No data for ", symbol_, " → skipping.")
      #return(NULL)
  #  }
  #)
  
  # store results
  #if(!is.null(quot_)) {
   # quot_$symbol <- symbol_  # tag with symbol
    #if(is.null(quotations_secondtry)) {
    #  quotations_secondtry <- quot_
   # } else {
    #  quotations_secondtry <- bind_rows(quotations_secondtry, quot_)
   # }
 # }
#}
