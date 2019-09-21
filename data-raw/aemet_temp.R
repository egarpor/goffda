
library(fda.usc)

# Load the results from running aemet.local.R used for creating fda.usc's 
# aemet. The script reads the data downloaded from AEMET's FTP in 2014. 
# Unfortunately, the data retrieval approach does not work anymore due to a
# change in AEMET's data access API, so one must have locally the old raw files
# to reconstruct aemet.raw.RData
load("aemet.raw.RData")

# Fix duplicated ID in "A CORUÑA AEROPUERTO" station
aemet.raw$df$ind[aemet.raw$df$name == "A CORUÑA AEROPUERTO"] <- "1387A"

# Retain only the temperature
aemet_temp <- list("temp" = aemet.raw$temp, "df" = aemet.raw$df)

# Purge to use only the 73 stations in fda.usc::aemet
data(aemet)
aemet$df$ind[2] <- "1387A" # Fix duplicated ID for consistency
ind_purge <- aemet_temp$df$ind %in% aemet$df$ind
aemet_temp$df <- aemet_temp$df[ind_purge, ]
aemet_temp$temp <- aemet_temp$temp[ind_purge]

# Use the same names from fda.usc::aemet
ind_names <- match(x = aemet_temp$df$ind, table = aemet$df$ind)
aemet_temp$df$name <- aemet$df$name[ind_names]

# Remove useless info
aemet_temp$df$province <- NULL
aemet_temp$df$longitude <- NULL
aemet_temp$df$latitude <- NULL
aemet_temp$df$altitude <- NULL

# Add metainfo
aemet_temp$temp$rangeval <- c(0, 365)
aemet_temp$temp$names <- list(main = "Temperature", xlab = "Day", 
                              ylab = "Temperature (ºC)")

# Save dataset
save(file = "aemet_temp", file = "aemet_temp.rda", compress = "bzip2")
