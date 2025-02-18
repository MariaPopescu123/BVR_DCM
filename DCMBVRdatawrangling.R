#### section to start editing: 2/7 Adding PAR, DO, DOsat_percent, cond, ORP, pH, temp ####

# Maria DCM BVR
#data includes chlorophyll maxima, PAR, secchi, light attenuation, metals, ghgs, nutrients


pacman::p_load(tidyverse, lubridate, akima, reshape2, pracma,
               gridExtra, grid, colorRamps, RColorBrewer, rLakeAnalyzer,
               reader, cowplot, dplyr, tidyr, ggplot2, zoo, purrr, beepr, forecast, ggthemes, splines)

source("interpolate_variable.R")
source("data_availability_function.R")


#need to update links for all the data now

#### Loading Data  ####

#ctd data https://portal.edirepository.org/nis/metadataviewer?packageid=edi.200.14
#not updated for 2024
CTD <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/200/14/0432a298a90b2b662f26c46071f66b8a")

#flora data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.272.9
#nonupdated
current_df <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/272/8/0359840d24028e6522f8998bd41b544e")

# metals data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.455.8
metalsdf <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/455/8/9c8c61b003923f4f03ebfe55cea8bbfd")
#removed flags for 68 as per Cece's advice

#ghgs data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.551.8
ghgs <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/551/8/454c11035c491710243cae0423efbe7b")
#not sure whether or not to remove those with flags 3 and 4
#3 = The difference between the reps are above the limit of quantification and >30% and <50% different from each other. Both replicates were retained but flagged
#4 = The difference between the reps are above the limit of quantification and >50% different from each other. Both replicates were retained but flagged

#secchi data https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=198&revision=13
#updated 2025
secchiframe <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/198/13/3ee0ddb9f2183ad4d8c955d50d1b8fba")

#ysi https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=198&revision=13
#updated 2025
ysi_profiles <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/198/13/3ee0ddb9f2183ad4d8c955d50d1b8fba")

#data from here https://portal.edirepository.org/nis/mapbrowse?packageid=edi.199.12
chemistry <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/199/12/a33a5283120c56e90ea414e76d5b7ddb")

#meteorological data from FCR https://portal.edirepository.org/nis/mapbrowse?packageid=edi.389.8
options(timeout = 300)
metdata <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/389/8/d4c74bbb3b86ea293e5c52136347fbb0")

#bathymetry data for BVR https://portal.edirepository.org/nis/metadataviewer?packageid=edi.1254.1
bath <- read.csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.1254.1&entityid=f7fa2a06e1229ee75ea39eb586577184")

#waterlevel data
wtrlvl <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/725/4/43476abff348c81ef37f5803986ee6e1") 

#waterlevel data using the pressure sensor (platform data) https://portal.edirepository.org/nis/metadataviewer?packageid=edi.725.4
#for past 2020
BVRplatform <- read.csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.725.4&entityid=9adadd2a7c2319e54227ab31a161ea12")

####weekly dataframe for interpolation####
start_date <- as.Date("2014-01-01")
end_date <- as.Date("2024-12-31")

weekly_dates <- data.frame(
  Date_fake = seq.Date(from = start_date, to = end_date, by = "week")
) %>%
  mutate(Year = year(Date_fake),
         Week = week(Date_fake))|>
  mutate(Depth_m = NA)

Depth_fake = seq(0, 13, by = 0.1)

# Expand grid to get each date with each depth
expanded_dates <- expand_grid(Date_fake = weekly_dates$Date_fake, Depth_m = Depth_fake)

# Add year and week info to the expanded data
expanded_dates <- expanded_dates %>%
  mutate(Year = year(Date_fake),
         Week = week(Date_fake),
         DOY = yday(Date_fake), 
         Date = Date_fake)|>
  select(-Date_fake)

#adding columns with total_conc max and the depth at which it occurs
phytos <- current_df %>% 
  filter(Reservoir == "BVR", Site == 50)%>%
  mutate(Date  = as_date(DateTime)) |> 
  filter((hour(DateTime) >= 8), (hour(DateTime) <= 18))|>
  filter(!(CastID == 592))|> #filter out weird drop in 2017
  filter(!(CastID == 395))|> #weird drop in 2016
  #filter(Flag_TotalConc_ugL != 2,Flag_TotalConc_ugL != 3)|> #2 is instrument malfunction and #3 is low transmission value
  mutate(Week = week(Date))|>
  mutate(Year = year(Date))|>
  mutate(DOY = yday(Date))

phytos2018 <- current_df %>% 
  filter(Reservoir == "BVR", Site == 50)%>%
  mutate(Date  = as_date(DateTime)) |> 
  filter((hour(DateTime) >= 8), (hour(DateTime) <= 18))|>
  filter(!(CastID == 592))|> #filter out weird drop in 2017
  filter(!(CastID == 395))

filtered<- 

#write.csv(phytos, "phytos.csv", row.names = FALSE)



####flora instrument data availability####

#days on the x axis, years on the y axis
plot_dat <- phytos %>%
  filter(!is.na(TotalConc_ugL)) %>%
  mutate(Year = year(Date), 
         DayOfYear = yday(Date))|> # Extract year and day of the year
  select(Date, Year, DayOfYear, TotalConc_ugL, Depth_m)

# Find the maximum TotalConc_ugL value for each year
max_totals_per_year <- plot_dat %>%
  group_by(year(Date)) %>%
  slice(which.max(TotalConc_ugL)) %>%
  ungroup()

# Plot: x-axis is DayOfYear, y-axis is Year, with a line and highlighted points
ggplot(plot_dat, aes(x = DayOfYear, y = as.factor(Year), group = Year)) +
  geom_line() +  # Line for each year
  geom_point() +  # Data points
  geom_point(data = max_totals_per_year, aes(x = DayOfYear, y = as.factor(Year)), 
             color = "red", size = 3) +  # Highlight max points in red
  geom_text(data = max_totals_per_year, 
            aes(x = DayOfYear, y = as.factor(Year), 
                label = paste0("Max: ", round(TotalConc_ugL, 2), " µg/L\nDepth: ", Depth_m, " m")), 
            vjust = 1.5, hjust = 0.5, color = "black", size = 3) +  # Smaller text and place below the point
  theme_bw() +
  labs(x = "Day of Year", y = "Year", title = "Fluoroprobe Data Availability") +
  scale_x_continuous(breaks = seq(1, 365, by = 30), limits = c(1, 365)) +  # Set x-axis limits and breaks
  theme(panel.grid.minor = element_blank())+  # Optional: remove minor grid lines
  geom_vline(xintercept = 133, linetype = "dashed", color = "red") +  # Vertical dashed line at DayOfYear 133
  geom_vline(xintercept = 286, linetype = "dashed", color = "red")  # Vertical dashed line at DayOfYear 286

#see that the data is dispersed at random intervals (and also different depths)
#because of this going to interpolate and join this to a weekly dataframe that has all the same depths


#list of DOY for interpolation purpose
DOY_list <- 32:334  # DOYs from February 1 to November 30
years <- unique(year(wtrlvl$Date))
DOY_year_ref <- expand.grid(Year = years, DOY = DOY_list)|>
  arrange(Year, DOY)

#add water level to data frame to use as the max depth for creating sequence of depths to interpolate each cast to
####Waterlevel####
wtrlvl <- wtrlvl |> 
  mutate(Date = as.POSIXct(DateTime))

#Add DOY and Year columns to wtrlvl2, then join with DOY_year_ref
wtrlvl2 <- wtrlvl |>
  mutate(Year = year(Date), DOY = yday(Date))

#join and interpolate WaterLevel_m for each DOY in each year
wtrlvl2_interpolated <- DOY_year_ref %>%
  left_join(wtrlvl2, by = c("Year", "DOY")) %>%
  group_by(Year) %>%
  mutate(
    WaterLevel_m = zoo::na.spline(WaterLevel_m, x = DOY, na.rm = FALSE)
  ) %>%
  filter(Year > 2013) %>%
  arrange(Year, DOY)|>
  select(Year, DOY, WaterLevel_m)

#now for past 2020 
#Add DOY and Year columns to wtrlvl2, then join with DOY_year_ref
BVRplatform2 <- BVRplatform |>
  filter(Flag_LvlPressure_psi_13 != 5)|>#filter flags, questionable value but left in the dataset
  mutate(Date = as.Date(DateTime))|>
  mutate(Year = year(Date), DOY = yday(Date))

#join and interpolate WaterLevel_m for each DOY in each year
BVRplatform2_interpolated <- DOY_year_ref |>
  left_join(BVRplatform2, by = c("Year" = "Year", "DOY" = "DOY")) |>
  group_by(Year) |>
  mutate(
    LvlDepth_m_13 = zoo::na.spline(LvlDepth_m_13, x = DOY, na.rm = FALSE)
  )|>
  filter(Year > 2019, Site == 50)|>
  arrange(Year, DOY)|>
  select(Year, DOY, DateTime, LvlDepth_m_13)

water_levelsjoined <- expanded_dates|>
  left_join(BVRplatform2_interpolated, by = c("Year", "DOY"), relationship = "many-to-many")|>
  left_join(wtrlvl2_interpolated, by = c("Year", "DOY"), relationship = "many-to-many")|>
  filter(Year>2013)

water_levelscoalesced<- water_levelsjoined|>
  mutate(WaterLevel_m = coalesce(LvlDepth_m_13,WaterLevel_m))|>
  select(Year, DOY, WaterLevel_m)|>
  group_by(Year, DOY)|>
  summarise(WaterLevel_m = mean(WaterLevel_m, na.rm = TRUE), .groups = "drop")


#this has all the depths and all the days
water_levels <- expanded_dates|>
  left_join(water_levelscoalesced, by = c("Year", "DOY"))

ggplot(water_levels, aes(x = Date, y = WaterLevel_m)) +
  geom_line(color = "#2C3E50", size = .8) +  # Use a sophisticated dark blue-gray color
  theme_minimal(base_size = 14) +  # Increase base font size for readability
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.3),  # Subtle grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines for a cleaner look
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(color = "black"),  # Dark axis text for contrast
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)  # Centered bold title
  ) +
  labs(
    title = "Water Level Over Time",
    x = "Date",
    y = "Water Level (m)"
  )

####phyto interp####
#1. first bind it to expanded_dates (this has all the depths and all the weeks I want to interpolate to)
#2. then add it to the last completed dataframe (in this case water_levels)

#testing interpolation funciton on phytos
variables <- ("TotalConc_ugL")
phytos_interpolated <- interpolate_variable(phytos, variables, expanded_dates)

# looking<- phytos_interpolated|>
#  filter(!is.na(TotalConc_ugL))

####interp_phyto data availability####

#days on the x axis, years on the y axis
plot_dat <- phytos_interpolated %>%
  filter(!is.na(TotalConc_ugL)) %>%
  mutate(Year = year(Date), 
         DayOfYear = yday(Date))|> # Extract year and day of the year
  select(Date, Year, DayOfYear, TotalConc_ugL, Depth_m)

# Find the maximum TotalConc_ugL value for each year
max_Totals_per_year <- plot_dat %>%
  group_by(year(Date)) %>%
  slice(which.max(TotalConc_ugL)) %>%
  ungroup()

# Plot: x-axis is DayOfYear, y-axis is Year, with a line and highlighted points
ggplot(plot_dat, aes(x = DayOfYear, y = as.factor(Year), group = Year)) +
  geom_line() +  # Line for each year
  geom_point() +  # Data points
  geom_point(data = max_Totals_per_year, aes(x = DayOfYear, y = as.factor(Year)), 
             color = "red", size = 3) +  # Highlight max points in red
  geom_text(data = max_Totals_per_year, 
            aes(x = DayOfYear, y = as.factor(Year), 
                label = paste0("Max: ", round(TotalConc_ugL, 2), " µg/L\nDepth: ", Depth_m, " m")), 
            vjust = 1.5, hjust = 0.5, color = "black", size = 3) +  # Smaller text and place below the point
  theme_bw() +
  labs(x = "Day of Year", y = "Year", title = "Interpolated Data Distribution") +
  scale_x_continuous(breaks = seq(1, 365, by = 30), limits = c(1, 365)) +  # Set x-axis limits and breaks
  theme(panel.grid.minor = element_blank())+  # Optional: remove minor grid lines
  geom_vline(xintercept = 133, linetype = "dashed", color = "red") +  # Vertical dashed line at DayOfYear 133
  geom_vline(xintercept = 286, linetype = "dashed", color = "red")  # Vertical dashed line at DayOfYear 286

#add it to the water_levels 
phytos_waterlevel<- water_levels|>
  left_join(phytos_interpolated, by = c("DOY", "Year", "Depth_m", "Week", "Date"))|>
  group_by(Year, Week) %>%
  mutate(Totals_DCM_conc = max(TotalConc_ugL, na.rm = TRUE))|> #concentration of totals at totals DCM
  mutate(Totals_DCM_depth = ifelse(TotalConc_ugL == Totals_DCM_conc, Depth_m, NA_real_))|>
  fill(Totals_DCM_conc, .direction = "downup")|>
  fill(Totals_DCM_depth, .direction = "downup")|>
  ungroup()


#### metals  ####
metalsdf_filtered <- metalsdf |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  filter(!if_any(starts_with("Flag"), ~. == 68))

#how much raw metal data do we have?

variables <- c("SFe_mgL", "TFe_mgL", "SMn_mgL", "SCa_mgL",
               "TCa_mgL", "TCu_mgL", "SCu_mgL", "SBa_mgL", "TBa_mgL")

#data_availability(metalsdf_filtered, variables)

metals_interpolated <- interpolate_variable(metalsdf_filtered, variables, expanded_dates)

#data_availability(metals_interpolated, variables)


phytos_wtrlvl_metals <- phytos_waterlevel|>
  left_join(metals_interpolated, by = c("DOY", "Year", "Depth_m", "Week", "Date"))

#### ghgs  ####
ghgs_filtered <- ghgs |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))

#see how much raw ghg data is available
variables <- c("CO2_umolL", "CH4_umolL")

data_availability(ghgs_filtered, variables)

ghgs_interpolated <- interpolate_variable(ghgs_filtered, variables, expanded_dates)

phytos_wtrlvl_metals_ghgs <- phytos_wtrlvl_metals|>
  left_join(ghgs_interpolated, by = c("DOY", "Year", "Depth_m", "Week", "Date"))

#### secchi PZ  ####
{
  
  secchi_df <- secchiframe |>
    mutate(Date = as_date(DateTime)) |>
    group_by(Date, Reservoir, Site) |>
    summarise(Secchi_m = mean(Secchi_m, na.rm = TRUE), .groups = "drop") |>
    filter(Reservoir == "BVR" & Site == 50) |>
    mutate(Year = year(Date), DOY = yday(Date))
  
  looking <- secchi_df|>
    filter(Year == 2018)
  
  variables <- c("Secchi_m")
  
  data_availability(secchi_df, variables) #see how much raw secchi data is available
  
  # Ensure DOY_year_ref contains all DOY values for each Year
  DOY_year_ref <- expand.grid(Year = unique(secchi_df$Year), DOY = 1:366) |> # Handle leap years
    filter(!(Year %% 4 != 0 & DOY == 366)) |>  # Remove DOY 366 for non-leap years
    mutate(Date = as.Date(DOY - 1, origin = paste0(Year, "-01-01")))  # Ensure Date matches DOY
  
  # Perform interpolation
  secchi_interpolated <- DOY_year_ref %>%
    left_join(secchi_df, by = c("Year", "DOY")) %>%
    filter(Year > 2013, DOY != 207) %>% #filter out incorrect secchi obseration (doesn't make sense) 
    group_by(Year) %>%
    mutate(
      first_valid_DOY = min(DOY[!is.na(Secchi_m)], na.rm = TRUE),
      last_valid_DOY = max(DOY[!is.na(Secchi_m)], na.rm = TRUE),
      Secchi_m = ifelse(
        DOY >= first_valid_DOY & DOY <= last_valid_DOY,
        na.approx(Secchi_m, x = DOY, na.rm = FALSE),
        NA_real_
      )
    ) |>
    arrange(Year, DOY) |>
    select(Year, DOY, Secchi_m)
           

# Adding Secchi
pwmgs <- phytos_wtrlvl_metals_ghgs|> #first letter of each dataframe for traceability
  left_join(secchi_interpolated, by = c("Year", "DOY"))|>
  group_by(Date)|>
  fill(Secchi_m, .direction = "updown")|>
  ungroup()

data_availability(pwmgs, variables)

# Calculating K_d and light availability from secchi

pwmgsl <- pwmgs |> #add light
  mutate(sec_K_d = 1.7/Secchi_m) |>
  mutate(light_availability_fraction = exp(-sec_K_d * Depth_m)) |>
  mutate(sec_LAP = light_availability_fraction * 100)|> #light availability percentage calculated from secchi
  mutate(PZ = 4.065 /sec_K_d)|>
  group_by(Date)|>
  mutate(PZ = if_else(PZ > WaterLevel_m, WaterLevel_m, PZ))
}

####Adding PAR, DO, DOsat_percent, cond, ORP, pH, temp ####
#####YSI#####
ysi_profiles <- ysi_profiles|>
  mutate(Date = as_date(DateTime))

#come back to remove flags
variables <- c("DO_mgL", "PAR_umolm2s", "DOsat_percent", "Cond_uScm", "ORP_mV", "pH", "Temp_C")

data_availability(ysi_profiles, variables)

# Generate the plot
plot <- data_availability(ysi_profiles, variables)  
# Save the plot with specific dimensions
ggsave("raw_ysi_availability.png", plot = plot, width = 20, height = 15, dpi = 300)


#removing PAR, ORP, cond, and pH due to limited data availability
#keeping temp because YSI has the most temp

variables <- c("DO_mgL","DOsat_percent", "Temp_C")
ysi <- ysi_profiles|>
  select(-PAR_umolm2s, -ORP_mV, -Cond_uScm, -pH)|>
  filter(Reservoir == "BVR", Site == 50)
  
ysi_interpolated <- interpolate_variable(ysi, variables, expanded_dates)

data_availability(ysi_interpolated, variables)

pwmgsly <- pwmgsl|>
  left_join(ysi_interpolated, by = c("DOY", "Year", "Depth_m", "Week", "Date"))

#look at new interpolated values for DO and DOsat
data_availability(ysi_interpolated, variables)
# Generate the plot
plot <- data_availability(ysi_interpolated, variables)  
# Save the plot with specific dimensions
ggsave("data_availability_plot.png", plot = plot, width = 8, height = 6, dpi = 300)


#####CTD#####
CTDfiltered <- CTD|> #flag 2, instrument malfunction. haven't removed flags yet
  filter(Reservoir == "BVR", Site == 50)|>
  filter(!if_any(starts_with("Flag"), ~. == 68))|>
  mutate(Date = as_date(DateTime))|>
  filter(Reservoir == "BVR", Site == 50)
  

variables <- c("DO_mgL", "PAR_umolm2s", "DOsat_percent", "Cond_uScm", "ORP_mV", 
               "pH", "Temp_C")
plot <- data_availability(CTDfiltered, variables)
#ggsave("raw_CTD_availability.png", plot = plot, width = 20, height = 15, dpi = 300)

#not going to use CTD
#variables <- c("Cond_uScm")
#CTDinterpolated <- interpolate_variable(CTDfiltered, variables, expanded_dates)
#plot <- data_availability(CTDinterpolated, variables)
#ggsave("interpolated_CTD_availability.png", plot = plot, width = 20, height = 15, dpi = 300)

#join to dataframe with everything
#pwmgslyc <- pwmgsly|>
# left_join(CTDinterpolated, by = c("DOY", "Year", "Depth_m", "Week", "Date"))

#need 2014 temp data from CTD
CTDfiltered2014 <- CTDfiltered|>
  filter(year(Date) == 2014)
variables <- c("Temp_C")
CTDinterpolated <- interpolate_variable(CTDfiltered2014, variables, expanded_dates)
pwmgslyc <- pwmgsly|>
  left_join(CTDinterpolated, by = c("DOY", "Year", "Depth_m", "Week", "Date"))|>
  mutate(Temp_C = coalesce(Temp_C.x, Temp_C.y))|>
  select(-Temp_C.x, -Temp_C.y)

data_availability(pwmgslyc, variables)


#### Nutrients  ####

chemistry_filtered <- chemistry |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))


variables <- c("TN_ugL", "TP_ugL", "NH4_ugL", "NO3NO2_ugL", "SRP_ugL", 
               "DOC_mgL", "DIC_mgL", "DC_mgL", "DN_mgL")
#raw data availability 
plot <- data_availability(chemistry_filtered, variables)
ggsave("raw_chem_availability.png", plot = plot, width = 20, height = 15, dpi = 300)

#interpolated data availability 
chemistry_interpolated <- interpolate_variable(chemistry_filtered, variables, expanded_dates)
#plot <- data_availability(chemistry_interpolated, variables)
#ggsave("interpolated_chem_availability.png", plot = plot, width = 20, height = 15, dpi = 300)

pwmgslycc <- pwmgslyc|>
  left_join(chemistry_interpolated, by = c("DOY", "Year", "Depth_m", "Week", "Date"))


#### NP ratio  ####

calculate_np_ratio <- function(tn, tp) {
  # Convert concentrations from µg/L to mg/L
  tn_mgL <- tn / 1000
  tp_mgL <- tp / 1000
  
  # Convert mg/L to moles/L
  tn_molL <- tn_mgL / 14.01
  tp_molL <- tp_mgL / 30.97
  
  # Calculate N:P ratio
  calcnp_ratio <- ifelse(is.na(tn_molL) | is.na(tp_molL), NA, tn_molL / tp_molL)
  
  return(calcnp_ratio)
}

# added np ratio to dataframe
pwmgslyccn <- pwmgslycc %>%
  mutate(np_ratio = calculate_np_ratio(TN_ugL,TP_ugL))|>
  relocate(np_ratio, .before = TN_ugL)

#### Visualizing metdata  ####

metdata0 <- metdata|>
  mutate(Date = as_date(DateTime))|>
  mutate(DOY = yday(Date))|>
  relocate(DOY, .before = DateTime)|>
  relocate(Date, .before = DateTime)

#### Function for plotting meteorological variables #### 
metplots <- function(yearz, variable, maxx = NULL){
  
  metviz <- metdata0|>
    filter(year(DateTime) == yearz) #filtering for the year specified
  
  
  ggplot(metviz, aes(x = DateTime, y = {{variable}}))+
    geom_path()+
    ggtitle(paste(deparse(substitute(variable)), yearz))+
    theme_minimal()+
    scale_y_continuous(limits = c(0, maxx))  # setting consistent y-axis limits
  
}

#### Precipitation ####
#not including all of this for now
metdataprecip <- metdata0 |> 
  group_by(Date, year(DateTime))|> 
  mutate(precip_daily = sum(Rain_Total_mm, na.rm = TRUE))|>
  ungroup()|>
  relocate(Date, .before = DateTime)|>
  relocate(precip_daily, .before = Rain_Total_mm)

#b1 <- metplots(2015, precip_daily, maxx = 80)
#b2 <- metplots(2016, precip_daily, maxx = 80)
#b3 <- metplots(2017, precip_daily, maxx = 80)
#b4 <- metplots(2018, precip_daily, maxx = 80)
#b5 <- metplots(2019, precip_daily, maxx = 80)
#b6 <- metplots(2020, precip_daily, maxx = 80)
#b7 <- metplots(2021, precip_daily, maxx = 80)
#b8 <- metplots(2022, precip_daily, maxx = 80)
#b9 <- metplots(2023, precip_daily, maxx = 80)

#precips<- plot_grid(
#  b1, b2, b3,
#  b4, b5, b6, 
#  b7, b8, b9,
#  ncol = 3
#)
#print(b1)

#print(precips)

#### Air temps #### 

#b1 <- metplots(2015, AirTemp_C_Average, maxx = 50)
#b2 <- metplots(2016, AirTemp_C_Average, maxx = 50)
#b3 <- metplots(2017, AirTemp_C_Average, maxx = 50)
#b4 <- metplots(2018, AirTemp_C_Average, maxx = 50)
#b5 <- metplots(2019, AirTemp_C_Average, maxx = 50)
#b6 <- metplots(2020, AirTemp_C_Average, maxx = 50)
#b7 <- metplots(2021, AirTemp_C_Average, maxx = 50)
#b8 <- metplots(2022, AirTemp_C_Average, maxx = 50)
#b9 <- metplots(2023, AirTemp_C_Average, maxx = 50)


#temps<- plot_grid(
#  b1, b2, b3,
#  b4, b5, b6, 
#  b7, b8, b9,
#  ncol = 3
#)

#print(temps)

#### dailyaverage and dailymax for temps #### 
metdatatemps <- metdataprecip |> 
  group_by(Date, year(DateTime))|> 
  mutate(daily_airtempavg = mean(AirTemp_C_Average, na.rm = TRUE))|>
  mutate(maxdaily_airtemp = max(AirTemp_C_Average, na.rm = TRUE))|>
  mutate(mindaily_airtemp = min(AirTemp_C_Average, na.rm = TRUE))|>
  ungroup()|>
  relocate(daily_airtempavg, .before = AirTemp_C_Average)|>
  relocate(maxdaily_airtemp, .before = AirTemp_C_Average)|>
  relocate(mindaily_airtemp, .before = AirTemp_C_Average)|>
  select(Date,daily_airtempavg, maxdaily_airtemp, mindaily_airtemp, precip_daily)


#### precip and temp to final_datanpratio #### 

metdata_join <- metdatatemps |> 
  group_by(Date) |> 
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) |> 
  distinct()

final_datamet <- final_datanpratio|>
  left_join(metdata_join, by = c("Date"), relationship = "many-to-many")

  
####thermocline ####

# Dataframe with thermocline
just_thermocline <- pwmgslyccn |>
  filter(!is.na(Temp_C)) |>
  group_by(Date) |>
  group_modify(~ {
    # Calculate max depth where Temp_C is not NA
    max_depth <- max(.x$Depth_m[!is.na(.x$Temp_C)], na.rm = TRUE)
    
    # Compute thermocline depth
    thermocline_depth <- thermo.depth(
      .x$Temp_C, 
      .x$Depth_m, 
      Smin = 1, 
      seasonal = TRUE, 
      index = FALSE, 
      mixed.cutoff = 0.25 * max_depth  # Use 25% of max depth for mixed layer
    )
    
    # Apply cutoff: Remove thermocline depths deeper than 40% of max depth
    thermocline_depth <- ifelse(thermocline_depth > 0.6 * max_depth, NA, thermocline_depth)
    thermocline_depth <- ifelse(thermocline_depth < 0.25 * max_depth, NA, thermocline_depth)
    
    # Ensure Date is a single value and return the summarised dataframe
    tibble(Date = .x$Date[1], thermocline_depth = thermocline_depth)
  }) |>
  ungroup()

#fixing incorrect thermoclines
# weird_fixed<- pwmgslyccnt|>
#   filter(thermocline_depth<3)|>
#   group_by(Date)|>
#   filter(Depth_m > thermocline_depth)|>
#   mutate(thermocline_depth = thermo.depth(Temp_C, 
#                                           Depth_m, 
#                                           Smin = 2, 
#                                           seasonal = TRUE, 
#                                           index = FALSE,
#                                           mixed.cutoff = 3
#   ))|>
#   ungroup()

#add to frame
# final_datathermocline <- final_datanpratio|>
#   left_join(both_merged, by = c("CastID", "Depth_m", "Temp_C"), relationship = "many-to-many")

pwmgslyccnt <- pwmgslyccn|>
  left_join(just_thermocline, by = c("Date"))|>
  group_by(Date)|>
  fill(thermocline_depth, .direction = "updown")|>
  ungroup()

#####individual date thermocline check####
plot_data <- pwmgslyccnt |>
  filter(Date %in% c("2015-06-10"))|> #change depths here to see a specific day and see if the thermocline matches up
  select(Date, Depth_m, thermocline_depth, Temp_C)
# Extract the thermocline depth for the specific date for the line
thermocline_depth_value <- unique(plot_data$thermocline_depth)
# Create the plot
ggplot(plot_data, aes(x = Temp_C, y = Depth_m)) +
  geom_point() +  # Add points
  geom_hline(yintercept = thermocline_depth_value, linetype = "dashed", color = "red") +  # Add the thermocline line
  scale_y_reverse() +  # Inverts the y-axis
  labs(x = "Temperature (°C)", y = "Depth (m)", title = "Thermocline Depth on ____") +
  theme_minimal()  # Optional: apply a clean theme

####Buoyancy Frequency ####

pwmgslyccntb <- pwmgslyccnt|>
  group_by(Date)|>
  mutate(buoyancy_freq = c(buoyancy.freq(Temp_C, Depth_m), NA))|>#added for padding for the last value
  relocate(buoyancy_freq, .before = thermocline_depth)
#need to make sure this makes sense

####Peak.width####
#use Totals_mean

for_peaks <- pwmgslyccntb|>
  group_by(Date) %>%
  mutate(
    totals_med = median(TotalConc_ugL, na.rm = TRUE),  # Calculate the median, excluding NA values
    totals_sd = sd(TotalConc_ugL, na.rm = TRUE),       # Calculate the standard deviation
    totals_mean = mean(TotalConc_ugL, na.rm = TRUE),   # Calculate the mean
    totals_mean_plus_sd = totals_mean + totals_sd,          # Calculate mean + sd
    peak.top = as.integer(Depth_m <= Totals_DCM_depth & TotalConc_ugL > totals_mean_plus_sd),  # Create binary indicator
    peak.bottom = as.integer(Depth_m >= Totals_DCM_depth & TotalConc_ugL > totals_mean_plus_sd),
    
    # Apply condition: If Totals_DCM_conc < 40, set peak.top and peak.bottom to 0
    peak.top = if_else(Totals_DCM_conc < 40, 0, peak.top),
    peak.bottom = if_else(Totals_DCM_conc < 40, 0, peak.bottom),
    
    # Replace peak.top and peak.bottom with Depth_m if indicator is 1
    peak.top = if_else(peak.top == 1, Depth_m, 0),
    peak.bottom = if_else(peak.bottom == 1, Depth_m, 0),
    
    # Get the minimum peak.top value, replace Inf with NA if all are NA or 0
    peak.top = if_else(any(peak.top != 0), 
                       min(peak.top[peak.top != 0], na.rm = TRUE), 
                       NA_real_),
    
    # Get the maximum peak.bottom value, replace -Inf with NA if all are NA or 0
    peak.bottom = if_else(any(peak.bottom != 0), 
                          max(peak.bottom[peak.bottom != 0], na.rm = TRUE), 
                          NA_real_),
    
    # Calculate peak width and handle infinite values by replacing them with NA
    peak.width = peak.bottom - peak.top,
    peak.width = if_else(is.na(peak.top) | is.na(peak.bottom), NA_real_, peak.width)
  ) %>%
  ungroup()  # Ungroup after mutations


####Peak.magnitude####

final_data_peaks <- for_peaks|>
  group_by(Date)|>
  mutate(peak.magnitude = max(TotalConc_ugL))|>
  ungroup()|>
  select(Date, Depth_m, totals_mean, totals_sd, totals_mean_plus_sd, peak.top, peak.bottom, peak.width, peak.magnitude) #this is unnecessary. saying how many totals there are at the DCM for total_conc


library(lubridate)
conflicts_prefer(dplyr::filter)
library(dplyr)

pwmgslyccntbp <- pwmgslyccntb |> #with peak calculations
  left_join(final_data_peaks, by = c("Date", "Depth_m")) |>
  mutate(peak.width = if_else(peak.width < .3*WaterLevel_m, peak.width, NA_real_)) |>
  group_by(Date, Depth_m) |>
  mutate(DCM = if_else(Depth_m == Totals_DCM_depth, TRUE, FALSE))|>
  filter(DOY > 133, DOY < 286)

####final dataframe####

#write.csv(final_data0,"./final_data0.csv",row.names = FALSE)

####DCM depth correlations####
#removed buoyancy_freq for now bc had -inf will come back to

# Create a vector of variable names that need to be summarized
depth_variables <- c("Temp_C", "np_ratio", "SFe_mgL", "TFe_mgL", 
                     "SMn_mgL", "SCa_mgL", "TCa_mgL", 
                     "TCu_mgL", "SBa_mgL", "TBa_mgL", 
                     "CO2_umolL", "CH4_umolL", "DO_mgL", 
                     "DOsat_percent", "TN_ugL", "TP_ugL", 
                     "NH4_ugL", "NO3NO2_ugL", "SRP_ugL", 
                     "DOC_mgL", "DIC_mgL", "DC_mgL")

# Initialize an empty list to store results
max_depths <- list()


DCM_final <- pwmgslyccntbp |>
  filter(Totals_DCM_conc > 20)

# Loop through the depth_variables to calculate the depth at which the maximum value occurs for each date
for (var in depth_variables) {
  DCM_final <- DCM_final |>
    left_join(
      pwmgslyccntbp |>
        filter(DOY > 133, DOY < 286)|>
        group_by(Date) |>
        summarise(
          !!paste0("max_depth_", var) := {
            # Check for non-NA values
            if (any(!is.na(.data[[var]]))) {
              Depth_m[which.max(.data[[var]])]  # Get depth of maximum value
            } else {
              NA_real_  # Return NA if all values are NA
            }
          },
          .groups = "drop"
        ),
      by = "Date"
    )
}

# Loop through the depth_variables to calculate the depth at which the minimum value occurs for each date
for (var in depth_variables) {
  DCM_final <- DCM_final |>
    left_join(
      pwmgslyccntbp |>
        filter(month(Date) >= 4, month(Date) < 10) |>
        group_by(Date) |>
        summarise(
          !!paste0("min_depth_", var) := {
            # Check for non-NA values
            if (any(!is.na(.data[[var]]))) {
              Depth_m[which.min(.data[[var]])]  # Get depth of minimum value
            } else {
              NA_real_  # Return NA if all values are NA
            }
          },
          .groups = "drop"
        ),
      by = "Date"
    )
}


# Finalize the DCM_final data frame
DCM_final <- DCM_final |>
  rename_with(~ gsub("max_depth_(.*)", "max_\\1_depth", .), starts_with("max_depth_"))|>
  rename_with(~ gsub("min_depth_(.*)", "min_\\1_depth", .), starts_with("min_depth_"))

#write.csv(DCM_final,"./DCM_final.csv",row.names = FALSE)

####correlation function####

correlations <- function(year1, year2) {
  DCM_final_cor <- DCM_final |>
    filter(year(Date) >= {{year1}}, year(Date) <= {{year2}}) |>
    filter(month(Date) > 4, month(Date) < 10) |>
    filter(Totals_DCM_conc > 20)
  
  drivers_cor <- cor(DCM_final_cor[,c(6:66)],
                     method = "spearman", use = "pairwise.complete.obs")
 
  list(drivers_cor = drivers_cor, DCM_final_cor = DCM_final_cor)

}

#cutoff 0.7
results <- correlations(2014, 2019)
final_data_cor_results <- results$drivers_cor
final_data_cor_results[lower.tri(final_data_cor_results)] = ""
final_data_cor <- results$DCM_final_cor
final_data_cor_results <- results$drivers_cor

final_data_cor_results[lower.tri(final_data_cor_results)] <- NA
diag(final_data_cor_results) <- NA

# Flatten the correlation matrix into a long format
final_data_cor_long <- as.data.frame(as.table(final_data_cor_results)) |>
  filter(!is.na(Freq))  # Remove NAs introduced by setting the lower triangle to NA

final_data_cor_long$Freq <- as.numeric(as.character(final_data_cor_long$Freq))

significant_correlations <- final_data_cor_long |> # Filter correlations based on the cutoff of 0.65
  filter(abs(Freq) >= 0.65) |>  # Apply cutoff for correlation
  arrange(desc(abs(Freq)))# Sort by absolute correlation values

colnames(significant_correlations) <- c("Variable1", "Variable2", "Correlation") # Rename columns for clarity

significant_correlations <- significant_correlations |>
  filter(Variable1 %in% c("Totals_DCM_depth"))|>
  filter(!Variable2 %in% c("peak.top", "peak.bottom"))|>
  mutate(Combined = paste(Variable1, "vs", Variable2))

# Plot to visualize correlations
ggplot(significant_correlations, aes(x = Correlation, y = reorder(Combined, Correlation))) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Use a bar plot
  geom_text(aes(label = round(Correlation, 2)), 
            position = position_stack(vjust = 0.5), 
            color = "white") +  # Add correlation values as text on bars
  labs(title = "Significant Correlations (Cutoff 0.65) 2014-2023",
       x = "Correlation Value",
       y = "Variable Pairs") +
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.y = element_text(size = 10),  # Adjust y-axis text size
        plot.title = element_text(hjust = 0.5))  # Center the title

####correlations across years,max day each year####

DCM_final_maxdays_cor<- DCM_final|>
  filter(Date %in% c("2014-08-13", "2015-08-08", "2016-06-16", "2017-07-20", "2018-08-16", "2019-06-06", "2020-09-16", "2021-08-09", "2022-08-01", "2023-07-31"))

maxdayscor <- cor(DCM_final_maxdays_cor[,c(2:64)], method = "spearman", use = "pairwise.complete.obs")

maxdayscor[lower.tri(maxdayscor)] <- NA
diag(maxdayscor) <- NA

# Flatten the correlation matrix into a long format
maxdayscor_long <- as.data.frame(as.table(maxdayscor)) |>
  filter(!is.na(Freq))  # Remove NAs introduced by setting the lower triangle to NA

maxdayscor_long$Freq <- as.numeric(as.character(maxdayscor_long$Freq))

significant_correlations <- maxdayscor_long |> # Filter correlations based on the cutoff of 0.65
  filter(abs(Freq) >= 0.7) |>  # Apply cutoff for correlation
  arrange(desc(abs(Freq)))# Sort by absolute correlation values

colnames(significant_correlations) <- c("Variable1", "Variable2", "Correlation") # Rename columns for clarity

significant_correlations <- significant_correlations |>
  filter(Variable1 %in% c("Totals_DCM_depth"))|>
  filter(!Variable2 %in% c("peak.top", "peak.bottom"))|>
  mutate(Combined = paste(Variable1, "vs", Variable2))


significant_correlations$Combined <- paste(significant_correlations$Variable1, significant_correlations$Variable2, sep = " - ")

ggplot(significant_correlations, aes(x = Correlation, y = reorder(Combined, Correlation))) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Use a bar plot
  geom_text(aes(label = round(Correlation, 2)), 
            position = position_stack(vjust = 0.5), 
            color = "white") +  # Add correlation values as text on bars
  labs(title = "Significant Correlations (Cutoff 0.65) 2014-2023",
       x = "Correlation Value",
       y = "Variable Pairs") +
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.y = element_text(size = 10),  # Adjust y-axis text size
        plot.title = element_text(hjust = 0.5))  # Center the title

####daily correlation, for choosing specific day####

#these are the days that the max TotalConc_ugL occurs. The biggest bloom. 
blooms <- DCM_final|>
  group_by(year(Date))|>
  mutate(bloommax = if_else(TotalConc_ugL == max(TotalConc_ugL), TRUE, NA_real_))|>
  ungroup()|>
  filter(bloommax == TRUE)|>
  group_by(Date)|>
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) |>
  ungroup()
  

#"2014-08-13" "2015-08-08" "2016-06-16" "2017-07-20" "2018-08-16" "2019-06-06" "2020-09-16" "2021-08-09" "2022-08-01" "2023-07-31"
#change date to see correlations for the singular day that max was the biggest
daily_cor <- pwmgslyccntbp|>
  filter(Date %in% c(""))|>#change this
  select("Depth_m", "TotalConc_ugL", "TotalConc_ugL", "SFe_mgL", "TFe_mgL", "SMn_mgL", "SCa_mgL",
         "TCa_mgL", "TCu_mgL", "SBa_mgL", "TBa_mgL",
         "CO2_umolL", "CH4_umolL", "DO_mgL",
         "DOsat_percent", "np_ratio", "TN_ugL", "TP_ugL", 
         "NH4_ugL", "NO3NO2_ugL", "SRP_ugL", "DOC_mgL", "DIC_mgL", 
         "DC_mgL", "Temp_C", "buoyancy_freq")

daily_cor_result <- cor(daily_cor[,c(6:32)], method = "spearman", use = "pairwise.complete.obs")
  
daily_cor_result[lower.tri(daily_cor_result)] = ""

daily_cor_long <- as.data.frame(as.table(daily_cor_result)) |>
  filter(!is.na(Freq))  # Remove NAs introduced by setting the lower triangle to NA

daily_cor_long$Freq <- as.numeric(as.character(daily_cor_long$Freq))

significant_correlations <- daily_cor_long |> # Filter correlations based on the cutoff of 0.65
  filter(abs(Freq) >= 0.65) |>  # Apply cutoff for correlation
  arrange(desc(abs(Freq)))# Sort by absolute correlation values

colnames(significant_correlations) <- c("Variable1", "Variable2", "Correlation") # Rename columns for clarity
significant_correlations_sorted <- significant_correlations[order(significant_correlations$Variable1), ] #variable 1 sorted alphabetically 

significant_correlations$Combined <- paste(significant_correlations$Variable1, significant_correlations$Variable2, sep = " - ")

significant_correlations <- significant_correlations |>
  filter(Variable1 %in% c("TotalConc_ugL", "TotalConc_ugL"))|>
  filter(!Variable2 %in% c("Depth_m"))|>
  mutate(Combined = paste(Variable1, "vs", Variable2))

ggplot(significant_correlations, aes(x = Correlation, y = reorder(Combined, Correlation))) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Use a bar plot
  geom_text(aes(label = round(Correlation, 2)), 
            position = position_stack(vjust = 0.5), 
            color = "white") +  # Add correlation values as text on bars
  labs(title = "Significant Correlations (Cutoff 0.65) 2019-06-06 (max 2019 conc, DCM depth 8.6)",
       x = "Correlation Value",
       y = "Variable Pairs") +
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.y = element_text(size = 10),  # Adjust y-axis text size
        plot.title = element_text(hjust = 0.5))  # Center the title


looking<- final_data0|>
  filter(Date %in% c("2019-06-06"))



####DCM depth every year####
# Find the maximum TotalConc_ugL value for each day

max_totals_per_day <- plot_dat %>%
  group_by(Date) %>%
  slice(which.max(TotalConc_ugL)) %>%
  filter(DayOfYear > 133, DayOfYear < 285, TotalConc_ugL > 20) |>
  ungroup()

plot <- ggplot(max_totals_per_day, aes(x = DayOfYear, y = Depth_m, group = Year)) +
  geom_line() +
  geom_point(data = max_totals_per_year, aes(x = DayOfYear, y = Depth_m), 
             color = "red", size = 6) +  # Increase point size
  geom_text(data = max_totals_per_year, 
            aes(x = DayOfYear, y = Depth_m, 
                label = paste0("Max: ", round(TotalConc_ugL, 2), " µg/L\nDepth: ", Depth_m, " m")), 
            vjust = -0.5, hjust = 0.5, color = "black", size = 4) +  # Increase text size
  theme_bw() +
  labs(x = "Day of Year", y = "Depth (m)", title = "DCM Depths Across Years (Only Showing Data with totals > 20)") +
  scale_y_reverse(limits = c(10, 0)) +  # Invert y-axis from 0 to 10
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +  # Adjust x-axis breaks
  facet_wrap(~ Year, ncol = 2) +  # Create separate panels for each year
  theme(
    text = element_text(size = 32),  # Double the size of all text
    axis.title = element_text(size = 34),  # Increase axis title size
    axis.text = element_text(size = 28),  # Increase axis label size
    strip.text = element_text(size = 20),  # Increase facet label size
    plot.title = element_text(size = 38, face = "bold"),  # Increase title size
    legend.text = element_text(size = 10),  # Increase legend text size
    legend.title = element_text(size = 28),  # Increase legend title size
    panel.grid.minor = element_blank()  # Optional: remove minor grid lines
  )

ggsave("DCM_depth_across_years.png", plot, width = 20, height = 15, dpi = 300)

####peak width every year####
plot <- ggplot(DCM_final, aes(x = DOY, y = peak.width, group = Year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of Year", y = "Peak Width (m)", title = "DCM Widths Across Years (Only Showing Data with totals > 20)") +
  scale_y_continuous(limits = c(0, 4)) +  
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +  # Adjust x-axis breaks
  facet_wrap(~ Year, ncol = 2) +  # Create separate panels for each year
  theme(
    text = element_text(size = 32),  # Double the size of all text
    axis.title = element_text(size = 34),  # Increase axis title size
    axis.text = element_text(size = 28),  # Increase axis label size
    strip.text = element_text(size = 20),  # Increase facet label size
    plot.title = element_text(size = 38, face = "bold"),  # Increase title size
    legend.text = element_text(size = 10),  # Increase legend text size
    legend.title = element_text(size = 28),  # Increase legend title size
    panel.grid.minor = element_blank()  # Optional: remove minor grid lines
  )

ggsave("Peak_width_across_years.png", plot, width = 20, height = 15, dpi = 300)


####peak magnitude####
max_totals_per_year <- DCM_final %>%
  group_by(Year) %>%
  slice(which.max(TotalConc_ugL)) %>%
  filter(DOY > 133, DOY < 285, TotalConc_ugL > 20) |>
  ungroup()  

plot <- ggplot(DCM_final, aes(x = DOY, y = Totals_DCM_conc, group = Year)) +
  geom_line() +
  theme_bw() +
  geom_point(data = max_totals_per_year, aes(x = DOY, y = TotalConc_ugL), 
             color = "red", size = 6) +  # Increase point size
  geom_text(data = max_totals_per_year, 
            aes(x = DOY, y = Depth_m, 
                label = paste0("Max: ", round(TotalConc_ugL, 2))), 
            vjust = -0.5, hjust = 0.5, color = "black", size = 5) +  # Adjust text position further
  labs(x = "Day of Year", y = "Peak Magnitude (m)", title = "Peak Magnitude Across Years (Only Showing Data with totals > 20)") +
  scale_y_continuous(limits = c(0, 400)) +  
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +  # Adjust x-axis breaks
  facet_wrap(~ Year, ncol = 5) +  # Create separate panels for each year
  theme(
    text = element_text(size = 32),  # Double the size of all text
    axis.title = element_text(size = 34),  # Increase axis title size
    axis.text = element_text(size = 28),  # Increase axis label size
    strip.text = element_text(size = 20),  # Increase facet label size
    plot.title = element_text(size = 38, face = "bold"),  # Increase title size
    legend.text = element_text(size = 10),  # Increase legend text size
    legend.title = element_text(size = 28),  # Increase legend title size
    panel.grid.minor = element_blank()  # Optional: remove minor grid lines
  )

ggsave("Peak Magnitude_across_years.png", plot, width = 30, height = 10, dpi = 300)

looking <- DCM_final|>
  filter(!is.na(peak.magnitude))|>
  select(Date, peak.magnitude)
####boxplots depth of DCM####

#need to use raw data for this to work 

#for june, july, august
boxplot_Data <- DCM_final |>
  filter(Totals_DCM_conc > 20) |>
  filter(month(Date)>5, month(Date)<9) |>
  mutate(Year = year(Date), Month = month(Date))|>
  group_by(Year, Month)|>
  mutate(monthly_avg = mean(Totals_DCM_conc))

# Calculate max_legend_value for the color scale limits
max_legend_value <- max(boxplot_Data$Totals_DCM_conc, na.rm = TRUE)

# Create the multi-panel boxplot with an overlay of colored points for Totals_DCM_conc
ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
                         y = Totals_DCM_depth, 
                         fill = monthly_avg)) +
  geom_boxplot() +  # Boxplot with filled colors based on Totals_DCM_conc
  facet_wrap(~ Year) +  # Create a panel for each year
  scale_fill_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to boxes
  scale_y_reverse(name = "DCM Depth (inverted)") +  # Reverse the y-axis
  ylim(10, 0) +  # Set the y-axis limits, reversing the range
  labs(x = "Month", y = "DCM Depth", fill = "Total's µg/L") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
#visualizing just one box per year

boxplot_Data <- DCM_final |>
  filter(Totals_DCM_conc > 20) |>
  mutate(DayOfYear = yday(Date))|>
  filter(DayOfYear>133, DayOfYear<286) |>
  mutate(Year = year(Date), Month = month(Date))

label_data <- boxplot_Data %>%
  group_by(Year) %>%
  summarise(n = n())  # Calculate the number of data points per year

# Plot with labels for the number of data points
ggplot(boxplot_Data, aes(x = factor(Year), y = Totals_DCM_depth)) +
  geom_boxplot() +
  geom_point(aes(color = Totals_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  scale_y_reverse(name = "DCM Depth (inverted)") +  # Reverse the y-axis
  ggtitle(label = "DCM Depths only displaying totals > 20") +
  ylim(10, 0) +  # Set the y-axis limits, reversing the range
  labs(x = "Year", y = "DCM Depth", color = "totals ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = label_data, aes(x = factor(Year), y = 0.5, label = paste0("n = ", n)), 
            vjust = -0.5)  # Add labels at the top of each column

####boxplot width of DCM####

boxplot_Data <- DCM_final |>
  filter(Totals_DCM_conc > 20) |>
  filter(month(Date)>5, month(Date)<9) |>
  mutate(Year = year(Date), Month = month(Date))|>
  filter(peak.width<2.5)

# Create the multi-panel boxplot with an overlay of colored points for Totals_DCM_conc
ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
                         y = peak.width)) +
  geom_boxplot() +
  geom_point(aes(color = Totals_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  facet_wrap(~ Year) +  # Create a panel for each year
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  labs(x = "Month", y = "Peak Width", color = "totals ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#one box per year
boxplot_Data <- DCM_final |>
  filter(Totals_DCM_conc > 20) |>
  mutate(DayOfYear = yday(Date))|>
  filter(DayOfYear>133, DayOfYear<286) |>
  mutate(Year = year(Date), Month = month(Date))|>
  filter(peak.width<2.5)

label_data <- boxplot_Data %>%
  group_by(Year) %>%
  summarise(n = n())  # Calculate the number of data points per year

ggplot(boxplot_Data, aes(x = factor(Year), y = peak.width)) +
  geom_boxplot() +
  geom_point(aes(color = Totals_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  ggtitle(label = "Peak Width only displaying totals > 20") +
  ylim(0, 5) +  # Set the y-axis limits
  labs(x = "Year", y = "Peak Width", color = "totals ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = label_data, aes(x = factor(Year), y = 4.5, label = paste0("n = ", n)), 
            vjust = -0.5)  # Adjust y position for labels at the top


####boxplots magnitude of DCM####

#for June-August

boxplot_Data <- DCM_final |>
    filter(Totals_DCM_conc > 20) |>
  filter(month(Date)>5, month(Date)<9) |>
  mutate(Year = year(Date), Month = month(Date))

ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
                         y = peak.magnitude)) +
  geom_boxplot() +
  geom_point(aes(color = Totals_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  facet_wrap(~ Year) +  # Create a panel for each year
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  labs(x = "Month", y = "Peak Magnitude", color = "totals ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#visualizing just one box per year

boxplot_Data <- DCM_final |>
  filter(Totals_DCM_conc > 20) |>
  mutate(DayOfYear = yday(Date))|>
  filter(DayOfYear>133, DayOfYear<286) |>
  mutate(Year = year(Date), Month = month(Date))

ggplot(boxplot_Data, aes(x = factor(Year), y = peak.magnitude)) +
  geom_boxplot() +
  geom_point(aes(color = Totals_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  ggtitle(label = "Peak Magnitudes only displaying totals > 20")+
  ylim(0, 150) +  # Set the y-axis limits, reversing the range
  labs(x = "Year", y = "Peak Magnitude", color = "totals ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### plot TC, DCM, and PZ####

ggplot(DCM_final, aes(x = Date)) +
  scale_y_reverse() +
  geom_line(aes(y = thermocline_depth, color = "Thermocline Depth")) +  # Line for thermocline_depth
  geom_line(aes(y = WaterLevel_m, color = "Water Level")) +  # Line for water_level
  geom_line(aes(y = PZ, color = "Photic Zone (PZ)")) +  # Line for PZ
  geom_line(aes(y = Totals_DCM_depth, color = "DCM Depth"), size = 1) +  # Line for Totals_DCM_depth
  facet_wrap(~ Year, scales = "free_x") +  # Facet by Year to avoid connecting across years
  scale_color_manual(values = c("Thermocline Depth" = "red", 
                                "Water Level" = "blue", 
                                "Photic Zone (PZ)" = "orange", 
                                "DCM Depth" = "green")) +  # Set custom colors
  labs(color = "Variable")  # Add a label to the legend


####RandomForest Anually####

#"We constructed a RF of 1500 trees for each of the two response
#variables using 1% PAR depth (m), DOC concentration (mg L21),
#thermocline depth (m), metalimnion thickness (m),
#buoyancy frequency at the thermocline (s21),
#lake surface area (log10(km2)), and maximum depth (log10(m))
#as predictors included in each analysis."
#Leach Patterns and Drivers

#prepare data for random forest

DCM_RF <- DCM_final |> 
  ungroup()|>
  select(Date, WaterLevel_m, thermocline_depth, Totals_DCM_conc, Totals_DCM_depth, PZ, 
         max_Temp_C_depth, max_np_ratio_depth, max_SFe_mgL_depth, max_TFe_mgL_depth, 
         max_SMn_mgL_depth, max_SCa_mgL_depth, max_TCa_mgL_depth, max_TCu_mgL_depth, 
         max_SBa_mgL_depth, max_TBa_mgL_depth, max_CO2_umolL_depth, max_CH4_umolL_depth, 
         max_DO_mgL_depth, max_DOsat_percent_depth, max_TN_ugL_depth, max_TP_ugL_depth, 
         max_NH4_ugL_depth, max_NO3NO2_ugL_depth, max_SRP_ugL_depth, max_DOC_mgL_depth, 
         max_DIC_mgL_depth, max_DC_mgL_depth, min_Temp_C_depth, min_np_ratio_depth, 
         min_SFe_mgL_depth, min_TFe_mgL_depth, min_SMn_mgL_depth, min_SCa_mgL_depth, 
         min_TCa_mgL_depth, min_TCu_mgL_depth, min_SBa_mgL_depth, min_TBa_mgL_depth, 
         min_CO2_umolL_depth, min_CH4_umolL_depth, min_DO_mgL_depth, min_DOsat_percent_depth, 
         min_TN_ugL_depth, min_TP_ugL_depth, min_NH4_ugL_depth, min_NO3NO2_ugL_depth, 
         min_SRP_ugL_depth, min_DOC_mgL_depth) |>
  group_by(Date) |>
  summarise(across(everything(), mean, na.rm = TRUE)) |>
  ungroup()

library(randomForest)
library(missForest)

set.seed(123)  # Setting seed for reproducibility

# Splitting data into training (70%) and testing (30%)
index <- sample(1:nrow(DCM_RF), size = 0.7 * nrow(DCM_RF))  # 70% training data
train_data <- DCM_RF[index, ]
test_data <- DCM_RF[-index, ]


#should run model on test dataset , currently mine is running on training. 
#training and test RMSE, MAE and Rsquare value important to show. 
#GLM-AED

# Remove non-numeric columns (excluding Date, Depth_m, Year, etc.)
non_numeric_columns <- sapply(train_data, function(x) !is.numeric(x) & !is.factor(x))
train_data_no_non_numeric <- train_data %>%
  select(-which(non_numeric_columns)) 

# Replace Inf and NaN with NA in all numeric columns
train_data_no_non_numeric <- train_data_no_non_numeric %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.) | is.nan(.), NA, .)))

# Remove columns with more than 75% NA values
train_data_no_na <- train_data_no_non_numeric %>%
  select(where(~ mean(is.na(.)) <= 0.25))  # Keep columns with ≤ 25% NA

# Remove remaining rows with any NA values
train_data_imputed_z <- train_data_no_na %>%
  na.omit()
# Removes any rows with NAs
  #for 2022 leave out waterlevel bc too many NAs
  #train_data_imputed_z <- train_data_imputed_z|>
  #  select(-WaterLevel_m)
  
  # Add the excluded non-numeric columns (e.g., Date) back to the imputed dataset
  model_rf <- randomForest(Totals_DCM_depth ~ ., data = train_data_imputed_z, ntree = 500, importance = TRUE)
  ######grid search CV function
  #look at hyperparameters for RandomForest tune as many as 
  #n_estimators, max_depth, minimum samples split, max_leaf nodes, and min_samples leaf
  #change ntrees to 100 when playing with parameters
  
  
  importance(model_rf)
  
  #visualize
  importance_df <- as.data.frame(importance(model_rf))
  importance_df <- rownames_to_column(importance_df, var = "Variable") # Convert row names to a column
  
  filtered_importance_df <- importance_df %>%
    filter(!is.na(`%IncMSE`), `%IncMSE` > 0)# Filter for positive %IncMSE values
  
#for all years
  filtered_importance_df <- filtered_importance_df |>
    filter(!Variable %in% c("peak.top", "pH", "min_pH_depth",  "min_Cond_uScm_depth", "max_Cond_uScm_depth", "peak.magnitude", "peak.bottom", "secchi_PZ", "PAR_PZ", "max_Cond_uScm_depth", "DayOfYear", "DOY", "min_Cond_uScm_dept", "min_CH4_umolL_depth", "max_CH4_umolL_depth"))
  
  
  
  # Create the plot
  ggplot(filtered_importance_df, aes(x = `%IncMSE`, y = reorder(Variable, `%IncMSE`))) +
    geom_point(color = "blue", size = 3) +
    labs(
      title = "Variable Importance based on % IncMSE 2014-2023",
      x = "% IncMSE",
      y = "Variables"
    ) +
    theme_minimal()
  
  
  
  
 
  
####making predictions####
  # Predict Totals_DCM_depth for test data
  test_predictions <- predict(model_rf, newdata = train_data_imputed_z)
  
  rmse <- sqrt(mean((test_predictions - train_data_imputed_z$Totals_DCM_depth)^2, na.rm = TRUE))
  print(paste("RMSE:", rmse))
  
  # Plot predicted vs. actual values
  ggplot(data = NULL, aes(x = train_data_imputed_z$Totals_DCM_depth, y = test_predictions)) +
    geom_point(color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = "Predicted vs Actual Totals_DCM_depth",
      x = "Actual Totals_DCM_depth",
      y = "Predicted Totals_DCM_depth"
    ) +
    theme_minimal()
  
  
 


####RF totals magnitude####
#####make dataframe#####
  
  # Create a vector of variable names that need to be summarized
  depth_variables <- c("Temp_C", "np_ratio", "SFe_mgL", "TFe_mgL", 
                       "SMn_mgL", "SCa_mgL", "TCa_mgL", 
                       "TCu_mgL", "SBa_mgL", "TBa_mgL", 
                       "CO2_umolL", "CH4_umolL", "DO_mgL", 
                       "DOsat_percent", "Cond_uScm", "ORP_mV", 
                       "pH", "TN_ugL", "TP_ugL", 
                       "NH4_ugL", "NO3NO2_ugL", "SRP_ugL", 
                       "DOC_mgL", "DIC_mgL", "DC_mgL")
  
  DCM_depths <- list()
  
  # Process to obtain the concentrations at Totals_DCM_depth for each variable
  DCM_final_conc <- final_data0 |>
    mutate(Date = as.Date(Date)) |>
    filter(month(Date) >= 4, month(Date) < 10) |>
    group_by(Date) |>
    mutate(DCM_buoyancy_freq = if_else(DCM == TRUE, buoyancy_freq, NA_real_)) |>
    fill(DCM_buoyancy_freq, .direction = "updown") |>
    summarise(
      DCM_buoyancy_freq = mean(DCM_buoyancy_freq, na.rm = TRUE),
      Totals_DCM_depth = mean(Totals_DCM_depth, na.rm = TRUE),
      Totals_DCM_conc = mean(Totals_DCM_conc, na.rm = TRUE),
      peak.top = mean(peak.top, na.rm = TRUE),
      peak.bottom = mean(peak.bottom, na.rm = TRUE),
      peak.width = mean(peak.width, na.rm = TRUE),
      peak.magnitude = mean(peak.magnitude, na.rm = TRUE),
      secchi_PZ = mean(secchi_PZ, na.rm = TRUE),
      PAR_PZ = mean(PAR_PZ, na.rm = TRUE),
      PZ = mean(PZ, na.rm = TRUE),
      Zeu = mean(Zeu, na.rm = TRUE),
      thermocline_depth = mean(thermocline_depth, na.rm = TRUE),
      WaterLevel_m = mean(WaterLevel_m, na.rm = TRUE),
      .groups = "drop"  # Ungroup to prevent grouping issues in following steps
    ) |>
    filter(Totals_DCM_conc > 20)
  
  # Loop through depth_variables to join the concentration at Totals_DCM_depth
  for (var in depth_variables) {
    DCM_final_conc <- DCM_final_conc |>
      left_join(
        final_data0 |>
          filter(month(Date) >= 4, month(Date) < 10) |>
          select(Date, Depth_m, !!sym(var)) |>
          rename(Concentration = !!sym(var)) |>
          group_by(Date) |>
          filter(Depth_m == DCM_final_conc$Totals_DCM_depth[match(Date, DCM_final_conc$Date)]) |>
          summarise(!!paste0(var, "_at_DCM") := mean(Concentration, na.rm = TRUE)),
        by = "Date"
      )
  }
  
  
  DCM_final_conc|>
    select(-peak.top, -peak.bottom, -peak.width, )
  
  
  
  #####run RF#####
  library(randomForest)
  library(missForest)
  
  #trying within a year
  yearDCM_final <- DCM_final_conc |>
    # filter(year(Date) == 2023) |>
    mutate(DOY = yday(Date)) |>
    select(where(~ mean(is.na(.)) < 0.5))
  
  
  # select(-DCM_buoyancy_freq)#just for 2021
  
  # List of columns to apply the interpolation to (excluding Date and DOY)
  cols_to_interpolate <- c()
  
  # Loop through each column name in yearDCM_final
  for (col in colnames(yearDCM_final)) {
    # Check if the column has at least 3 non-NA observations
    if (sum(!is.na(yearDCM_final[[col]])) >= 3) {
      cols_to_interpolate <- c(cols_to_interpolate, col)  # Add column to list if condition is met
    }
  }
  
  # If you want to exclude specific columns (e.g., Date and DOY):
  cols_to_exclude <- c("Date", "DOY")
  
  # Loop through and filter out the excluded columns
  cols_to_interpolate <- cols_to_interpolate[!cols_to_interpolate %in% cols_to_exclude]
  
  library(pracma)
  
  #this is not appropriate
  # Loop through each column and apply pchip interpolation
  yearDCM_final <- yearDCM_final[order(yearDCM_final$DOY), ]
  
  for (col in cols_to_interpolate) {
    # Identify rows with non-NA values for the current column
    non_na_rows <- !is.na(yearDCM_final[[col]])
    
    # Perform PCHIP interpolation only on non-NA values
    yearDCM_final[[col]][!non_na_rows] <- pracma::pchip(
      yearDCM_final$DOY[non_na_rows],           # DOY values where the column is not NA
      yearDCM_final[[col]][non_na_rows],        # Column values where not NA
      yearDCM_final$DOY[!non_na_rows]           # DOY values where the column is NA
    )
  }
  # 
  
  set.seed(123) # Setting seed for reproducibility
  index <- sample(1:nrow(yearDCM_final), size = 0.7 * nrow(yearDCM_final))  # 70% training data
  train_data <- yearDCM_final[index, ]
  test_data <- yearDCM_final[-index, ]
  non_numeric_columns <- sapply(train_data, function(x) !is.numeric(x) & !is.factor(x))
  train_data_no_non_numeric <- train_data %>% select(-which(non_numeric_columns))
  
  # Apply na.roughfix() to impute missing values in numeric and factor columns
  train_data_imputed <- na.roughfix(train_data_no_non_numeric)
  
  #z-transform
  train_data_imputed_z <- train_data_imputed %>%
    mutate(across(everything(), ~ scale(.)))
  
  #when running RF for all years remove buoyancy freq
  train_data_imputed_z<- train_data_imputed_z|>
    select(-DCM_buoyancy_freq)
  
  train_data_imputed_z <- train_data_imputed_z %>%
    na.omit()  # Removes any rows with NAs
  
  #for 2022 leave out waterlevel bc too many NAs
  #train_data_imputed_z <- train_data_imputed_z|>
  #  select(-WaterLevel_m)
  
  # Add the excluded non-numeric columns (e.g., Date) back to the imputed dataset
  model_rf <- randomForest(Totals_DCM_conc ~ ., data = train_data_imputed_z, ntree = 500, importance = TRUE)
  
  importance(model_rf)
  
  #visualize
  importance_df <- as.data.frame(importance(model_rf))
  importance_df <- rownames_to_column(importance_df, var = "Variable") # Convert row names to a column
  
  filtered_importance_df <- importance_df %>%
    filter(!is.na(`%IncMSE`), `%IncMSE` > 0)# Filter for positive %IncMSE values
  
  #for all years
  filtered_importance_df <- filtered_importance_df |>
    filter(!Variable %in% c("peak.top", "pH_at_DCM", "min_pH_depth",  "min_Cond_uScm_depth", "max_Cond_uScm_depth", "peak.magnitude", "peak.bottom", "secchi_PZ", "PAR_PZ", "max_Cond_uScm_depth", "DayOfYear", "DOY", "min_Cond_uScm_dept", "min_CH4_umolL_depth", "max_CH4_umolL_depth"))
  
  
  
  # Create the plot
  ggplot(filtered_importance_df, aes(x = `%IncMSE`, y = reorder(Variable, `%IncMSE`))) +
    geom_point(color = "blue", size = 3) +
    labs(
      title = "Variable Importance based on % IncMSE 2014-2023",
      x = "% IncMSE",
      y = "Variables"
    ) +
    theme_minimal()
  
  



