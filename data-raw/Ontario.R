
library(R.matlab)
library(fda.usc)
library(lubridate)

# Read data downloaded from David Benatia's personal website
# (https://davidbenatia.com/research/) as companion for the paper
# "Functional linear regression with functional response" on the Journal of
# Econometrics (Benatia et al., 2017)
data <- readMat("DATA.mat")

# From Benatia et al. (2017) 's description, the .mat file contains two main
# arrays:
# * DATA_Y (25 x 368 x 35). The functional response, electricity consumption in
#   Ontario, is DATA_Y(:, :, 1)
# * DATA_Z (73 x 368 x 35). The functional predictor, temperature in
#   Ontario, is DATA_Z(:, :, 1)
# * TIME_Y (25 x 368). The time variable for functional data.
# The remaining arrays of DATA_Y and DATA_Z contain the seasonal dummies,
# excep for the last array that contains a variable for year to separate the
# samples

# Store the temperature as fdata
temp <- fda.usc::fdata(mdata = t(data$DATA.Z[, , 1]), argvals = (-24):48,
                       names = list(main = "Temperature 2010-2014",
                                    xlab = "Hour = [Day - 1, Day, Day + 1]",
                                    ylab = "Temperature (ºC)"))

# Note that there are 25 argvals recording times in data$DATA.Y, accounting
# for the sequence of hours 0:24. That means that the last observation of the
# i-th day is the same as the first observation of the (i + 1)-th day
# Electricity as fdata
elec <- fda.usc::fdata(mdata = t(data$DATA.Y[, , 1]), argvals = 0:24,
                       names = list(main = "Electricity consumption 2010-2014",
                                    xlab = "Hour", ylab = "Power (GW)"))

# In order to add the time information to the data, we run in MATLAB:
# % Export TIME_Y (not readable through R.matlab)
# writetable(cell2table(cellstr(TIME_Y)), "time_y.txt")

# Then read the time_y.txt dataset
time_y <- read.table(file = "time_y.txt", header = TRUE, sep = ",")

# Extract only the first row, as it contains the days of the observations
# (the rows are the hours for each day)
x <- parse_date_time(x = as.matrix(time_y)[1, ], orders = "%d-%b!-%y %H:%M:%S")
tz(x) <- "America/Toronto"

# Add time information
df <- data.frame("date" = x, "weekday" = weekdays(x, abbreviate = TRUE))

# Create dataset
ontario <- list("temp" = temp, "elec" = elec, "df" = df)

# Save data
save(list = "ontario", file = "ontario.rda", compress = "bzip2")
