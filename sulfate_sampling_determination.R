#determing sulfate windows for sample analysis 2024

source("data_availability_function.R")

#published 2025
current_df <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/272/9/f246b36c591a888cc70ebc87a5abbcb7")

BVRphytos <- current_df %>% 
  filter(Reservoir == "BVR", Site == 50)%>%
  mutate(Date  = as_date(DateTime)) |> 
  filter((hour(DateTime) >= 8), (hour(DateTime) <= 18))|>
  filter(!(CastID == 592))|> #filter out weird drop in 2017
  filter(!(CastID == 395))|> #weird drop in 2016
  #filter(Flag_TotalConc_ugL != 2,Flag_TotalConc_ugL != 3)|> #2 is instrument malfunction and #3 is low transmission value
  mutate(Week = week(Date))|>
  mutate(Year = year(Date))|>
  mutate(DOY = yday(Date))|>
  filter(Year == 2024, TotalConc_ugL>25)

variables <- c("TotalConc_ugL")

data_availability(BVRphytos, variables)

FCRphytos <- current_df %>% 
  filter(Reservoir == "FCR", Site == 50)%>%
  mutate(Date  = as_date(DateTime)) |> 
  filter((hour(DateTime) >= 8), (hour(DateTime) <= 18))|>
  filter(!(CastID == 592))|> #filter out weird drop in 2017
  filter(!(CastID == 395))|> #weird drop in 2016
  #filter(Flag_TotalConc_ugL != 2,Flag_TotalConc_ugL != 3)|> #2 is instrument malfunction and #3 is low transmission value
  mutate(Week = week(Date))|>
  mutate(Year = year(Date))|>
  mutate(DOY = yday(Date))|>
  filter(Year == 2016, TotalConc_ugL>25)

variables <- c("TotalConc_ugL")

data_availability(FCRphytos, variables)

