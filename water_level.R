#Maria Popescu
#Water level for Beaverdam 

#waterlevel data
wtrlvl <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/725/4/43476abff348c81ef37f5803986ee6e1") 

#waterlevel data using the pressure sensor (platform data) https://portal.edirepository.org/nis/metadataviewer?packageid=edi.725.4
#for past 2020
BVRplatform <- read.csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.725.4&entityid=9adadd2a7c2319e54227ab31a161ea12")


#list of DOY for interpolation purpose
DOY_list <- 32:334  # DOYs from February 1 to November 30
years <- unique(year(wtrlvl$Date))
DOY_year_ref <- expand.grid(Year = years, DOY = DOY_list)|>
  arrange(Year, DOY)


#add water level to data frame to use as the max depth for creating sequence of depths to interpolate each cast to
####Waterlevel####
wtrlvl <- wtrlvl |> 
  mutate(Date = as_date(DateTime))

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
water_level <- expanded_dates|>
  left_join(water_levelscoalesced, by = c("Year", "DOY"))|>
  filter(!is.na(WaterLevel_m))

ggplot(water_level, aes(x = Date, y = WaterLevel_m)) +
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

write.csv(water_level, "water_level.csv", row.names = FALSE)
