#### section to start editing: 2/7 Adding PAR, DO, DOsat_percent, cond, ORP, pH, temp ####

# Maria DCM BVR
#data includes chlorophyll maxima, PAR, secchi, light attenuation, metals, ghgs, nutrients


pacman::p_load(tidyverse, lubridate, akima, reshape2, pracma,
               gridExtra, grid, colorRamps, RColorBrewer, rLakeAnalyzer,
               reader, cowplot, dplyr, tidyr, ggplot2, zoo, purrr, beepr, forecast, ggthemes, splines)

source("interpolate_variable.R")
source("data_availability_function.R")
source("weekly_sum_variables.R")

#need to update links for all the data now

#### Loading Data  ####

#ctd data https://portal.edirepository.org/nis/codeGeneration?packageId=edi.200.15&statisticalFileType=r
#updated 2025
CTD <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/200/15/9d741c9cced69cfd609c473ada2812b1")

#flora data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.272.9
#nonupdated
current_df <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/272/8/0359840d24028e6522f8998bd41b544e")

#published 2025
current_df <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/272/9/f246b36c591a888cc70ebc87a5abbcb7")

# metals data https://portal.edirepository.org/nis/codeGeneration?packageId=edi.455.9&statisticalFileType=r
#updated 2025
metalsdf <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/455/9/9a072c4e4af39f96f60954fc4f7d8be5")
#removed flags for 68 as per Cece's advice

#secchi data https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=198&revision=13
#updated 2025
secchiframe <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/198/13/3ee0ddb9f2183ad4d8c955d50d1b8fba")

#ysi https://portal.edirepository.org/nis/mapbrowse?scope=edi&identifier=198&revision=13
#updated 2025
ysi_profiles <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/198/13/e50a50d062ee73f4d85e4f20b360ce4f")

#data from here https://portal.edirepository.org/nis/mapbrowse?packageid=edi.199.12
chemistry <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/199/12/a33a5283120c56e90ea414e76d5b7ddb")

#meteorological data from FCR https://portal.edirepository.org/nis/mapbrowse?packageid=edi.389.8
#options(timeout = 300)
#metdata <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/389/8/d4c74bbb3b86ea293e5c52136347fbb0")

#bathymetry data for BVR https://portal.edirepository.org/nis/metadataviewer?packageid=edi.1254.1
#bath <- read.csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.1254.1&entityid=f7fa2a06e1229ee75ea39eb586577184")


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

#write.csv(phytos, "phytos.csv", row.names = FALSE)

####flora instrument data availability####

#days on the x axis, years on the y axis
plot_dat <- phytos %>%
  filter(!is.na(TotalConc_ugL)) %>%
  mutate(Year = year(Date), 
         DayOfYear = yday(Date))|> # Extract year and day of the year
  select(Date, Year, Week, DayOfYear, TotalConc_ugL, Depth_m)

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

#see that the data is dispersed at random intervals
####choosing casts and calculating peaks####
#1. Look at every cast for every year and remove casts that do not make sense
#2. Calculate peak metrics for each cast (peak depth, width, and magnitude)
#3. Visually check all casts and each metric to make sure it makes sense
#4. Average the peak metrics together (if appropriate)
#5. New data frame with one set of peak metrics for each week that we have data


#plot each one for cast selection 

# Make sure the Figs directory exists
if (!dir.exists("Figs")) {
  dir.create("Figs")
}

# Prepare your data with FacetID
DCM_metrics <- phytos |>
  select(Reservoir, Site, Date, Week, CastID, Depth_m, TotalConc_ugL) |>
  group_by(Reservoir, Date, Site) |>
  mutate(FacetID = paste(CastID, Reservoir, Site, Date, "Week", Week, sep = " ")) |>
  ungroup()

# Get unique years in the dataset
years <- unique(year(DCM_metrics$Date))

# Loop over each year
for (yr in years) {
  
  # Filter data for the year
  test <- DCM_metrics |>
    filter(year(Date) == yr)
  
  # Skip if there's no data
  if (nrow(test) == 0) next
  
  # Create plot
  plot_casts <- ggplot(test, aes(x = TotalConc_ugL, y = Depth_m)) +
    geom_path() +
    scale_y_reverse() +
    theme_bw() +
    facet_wrap(vars(FacetID)) +
    xlab("micrograms per liter") +
    ylab("Depth (m)") +
    ggtitle(paste(yr, "raw casts"))
  
  # Save plot
  ggsave(filename = paste0("Figs/", yr, "_raw_casts.png"),
         plot = plot_casts,
         width = 12,
         height = 10,
         dpi = 300)
}

#notes on casts
#casts to remove: 467, 814, 856, 920, 1149 

DCM_metrics_filtered <- DCM_metrics |>
  filter(!CastID %in% c(467, 814, 856, 920, 1149)) |>
  mutate(CastID = case_when(
    CastID == 485 ~ 484,  # Change 485 to 484
    CastID == 492 ~ 493,  # Combine 492 and 493
    CastID == 499 ~ 500,  # Combine 499 and 500
    CastID == 603 ~ 604,  # Combine 603 and 604
    TRUE ~ CastID          # Keep other values unchanged
  ))|>
  mutate(DOY = yday(Date), Year = year(Date))|>
  filter(DOY > 133, DOY < 285)

#question about how to address casts that don't qualify as a bloom. They won't calculate metrics
#will RandomForest assume a bloom just didn't happen?
#does there need to be a component that first predicts whether or not a bloom will happen?
#and then if yes 

#join waterlevel because this will be important for peak metrics
water_level <- read.csv("water_level.csv") |>
  mutate(Date = as_date(Date)) |>
  select(Week, Year, WaterLevel_m) |>
  group_by(Week, Year) |>
  summarise(WaterLevel_m = mean(WaterLevel_m, na.rm = TRUE), .groups = "drop")

DCM_metrics_filtered <- DCM_metrics_filtered|>
  left_join(water_level, by = c("Week", "Year"))

####Peak.depth and max_conc####
DCM_metrics_depth <- DCM_metrics_filtered|>
  group_by(CastID) %>%
  mutate(max_conc = max(TotalConc_ugL, na.rm = TRUE))|> #concentration of totals at totals DCM
  mutate(DCM_depth = ifelse(TotalConc_ugL == max_conc, Depth_m, NA_real_))|>
  fill(max_conc, .direction = "downup")|>
  fill(DCM_depth, .direction = "downup")|>
  ungroup()

####Peak.width####
#use Totals_mean
#calculate the metrics on the actual observed profiles 
#don't interpolate 

#Using mean as per Mary's recommendation

for_peaks <- DCM_metrics_depth|>
  group_by(CastID) %>%
  mutate(
    totals_med = median(TotalConc_ugL, na.rm = TRUE),  # Calculate the median, excluding NA values
    totals_mean = mean(TotalConc_ugL, na.rm = TRUE),   # Calculate the mean
    peak.top = as.integer(Depth_m <= DCM_depth & TotalConc_ugL >= totals_mean),  # Create binary indicator
    peak.bottom = as.integer(Depth_m >= DCM_depth & TotalConc_ugL >= totals_mean),
    # Apply condition: If Totals_DCM_conc < 30, set peak.top and peak.bottom to 0
    peak.top = if_else(max_conc < 30, 0, peak.top),
    peak.bottom = if_else(max_conc < 30, 0, peak.bottom),
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
  ungroup()|>  # Ungroup after mutations
  select(CastID, Depth_m, peak.top, peak.bottom, peak.width, totals_mean, totals_med)

final_DCM_metrics <- DCM_metrics_depth |> #with peak calculations
  left_join(for_peaks, by = c("CastID", "Depth_m")) |>
  #mutate(peak.width = if_else(peak.width < .3*WaterLevel_m, peak.width, NA_real_)) |>
  filter(DOY > 133, DOY < 286)

####visualize DCM metrics####
# Loop over each year
for (yr in years) {
  
  # Filter data for the year
  test <- final_DCM_metrics |>
    filter(year(Date) == yr)
  
  # Skip if there's no data
  if (nrow(test) == 0) next
  
  # Create plot
  plot_casts <- ggplot(test, aes(x = TotalConc_ugL, y = Depth_m)) +
    geom_path() +
    # Light blue grid lines for every whole meter
    geom_hline(yintercept = seq(0, max(test$Depth_m, na.rm = TRUE), by = 1), 
               color = "lightblue", linetype = "dotted", linewidth = 0.3) +
    # Horizontal lines for depths
    geom_hline(aes(yintercept = peak.top), linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = DCM_depth), linetype = "solid", color = "red") +
    geom_hline(aes(yintercept = peak.bottom), linetype = "dashed", color = "red") +
    # Vertical lines for concentrations
    geom_vline(aes(xintercept = totals_mean), color = "red") +
    scale_y_reverse(breaks = seq(0, max(test$Depth_m, na.rm = TRUE), by = 1)) +
    theme_bw() +
    facet_wrap(vars(FacetID)) +
    xlab("micrograms per liter") +
    ylab("Depth (m)") +
    ggtitle(paste(yr, "raw casts"))+
    geom_text(aes(label = round(DCM_depth, 1), x = Inf, y = DCM_depth), 
              color = "black", hjust = 1.1, size = 3)
  
  # Save plot
  ggsave(filename = paste0("Figs/", yr, "_raw_casts.png"),
         plot = plot_casts,
         width = 12,
         height = 10,
         dpi = 300)
}

#for now not going to keep peak calculations in analysis 
#will come back to this at a later date

final_DCM_metrics<- final_DCM_metrics|>
  select(-peak.top, -peak.bottom, -peak.width)|>#can remove this once have peak calculations figured out
  group_by(CastID)|>
  mutate(Q3 = quantile(TotalConc_ugL, 0.75)) #25% of data falls above this value 

####boxplots depth of DCM####

#need to use raw data for this to work 

#for june, july, august
boxplot_Data <- final_DCM_metrics |>
  filter(max_conc > 20) |>
  filter(month(Date)>5, month(Date)<9) |>
  mutate(Year = year(Date), Month = month(Date))|>
  group_by(Year, Month)|>
  mutate(monthly_avg = mean(max_conc))

# Calculate max_legend_value for the color scale limits
max_legend_value <- max(boxplot_Data$max_conc, na.rm = TRUE)

# Create the multi-panel boxplot with an overlay of colored points for Totals_DCM_conc
ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
                         y = DCM_depth, 
                         fill = monthly_avg)) +
  geom_boxplot() +  # Boxplot with filled colors based on Totals_DCM_conc
  facet_wrap(~ Year) +  # Create a panel for each year
  scale_fill_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to boxes
  scale_y_reverse(name = "DCM Depth (inverted)") +  # Reverse the y-axis
  ylim(10, 0) +  # Set the y-axis limits, reversing the range
  labs(x = "Month", y = "DCM Depth", fill = "Total's µg/L") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
#visualizing just one box per year

boxplot_Data <- final_DCM_metrics|>
  filter(max_conc > 20) |>
  mutate(DayOfYear = yday(Date))|>
  filter(DayOfYear>133, DayOfYear<286) |>
  mutate(Year = year(Date), Month = month(Date))|>
  select(CastID, Date,Year,Month, DCM_depth, max_conc)|>
  group_by(CastID, Date,Year,Month,)|>
  summarise(
    DCM_depth = mean(DCM_depth, na.rm = TRUE),
    max_conc = mean(max_conc, na.rm = TRUE),
    .groups = "drop"
  )
  
label_data <- boxplot_Data %>%
  group_by(Year) %>%
  summarise(n = n())  # Calculate the number of data points per year

# Plot with labels for the number of data points
ggplot(boxplot_Data, aes(x = factor(Year), y = DCM_depth)) +
  geom_boxplot() +
  geom_point(aes(color = max_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  scale_color_gradientn(
    colours = c("blue","cyan", "green","yellow", "red", "red3"),
    values = scales::rescale(c(0,40, 75, 100, 200, max_legend_value)),  # Make red hit at 150 µg/L
    na.value = "gray",
    limits = c(20, max_legend_value), 
    breaks = c(20, 100, 200, 300, 380)
  ) +  scale_y_reverse(name = "DCM Depth (inverted)") +  # Reverse the y-axis
  ggtitle(label = "DCM Depths only displaying totals > 20") +
  ylim(10, 0) +  # Set the y-axis limits, reversing the range
  labs(x = "Year", y = "DCM Depth", color = "totals ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(data = label_data, aes(x = factor(Year), y = 0.5, label = paste0("n = ", n)), 
            vjust = -0.5)  # Add labels at the top of each column

####boxplot width of DCM####

# boxplot_Data <- DCM_final |>
#   filter(Totals_DCM_conc > 20) |>
#   filter(month(Date)>5, month(Date)<9) |>
#   mutate(Year = year(Date), Month = month(Date))|>
#   filter(peak.width<2.5)

# # Create the multi-panel boxplot with an overlay of colored points for width of DCM
# ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
#                          y = peak.width)) +
#   geom_boxplot() +
#   geom_point(aes(color = Totals_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
#   facet_wrap(~ Year) +  # Create a panel for each year
#   scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
#   labs(x = "Month", y = "Peak Width", color = "totals ugL") +  # Label the legend
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# #one box per year
# boxplot_Data <- DCM_final |>
#   filter(Totals_DCM_conc > 20) |>
#   mutate(DayOfYear = yday(Date))|>
#   filter(DayOfYear>133, DayOfYear<286) |>
#   mutate(Year = year(Date), Month = month(Date))|>
#   filter(peak.width<2.5)
# 
# label_data <- boxplot_Data %>%
#   group_by(Year) %>%
#   summarise(n = n())  # Calculate the number of data points per year
# 
# ggplot(boxplot_Data, aes(x = factor(Year), y = peak.width)) +
#   geom_boxplot() +
#   geom_point(aes(color = Totals_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
#   scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
#   ggtitle(label = "Peak Width only displaying totals > 20") +
#   ylim(0, 5) +  # Set the y-axis limits
#   labs(x = "Year", y = "Peak Width", color = "totals ugL") +  # Label the legend
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   geom_text(data = label_data, aes(x = factor(Year), y = 4.5, label = paste0("n = ", n)), 
#             vjust = -0.5)  # Adjust y position for labels at the top
# 

####boxplots magnitude of DCM####

#for June-August

boxplot_Data <- final_DCM_metrics |>
  filter(max_conc > 20) |>
  filter(month(Date)>5, month(Date)<9) |>
  mutate(Year = year(Date), Month = month(Date))|>
  select(CastID, Date,Year,Month, DCM_depth, max_conc)|>
  group_by(CastID, Date,Year,Month)|>
  summarise(
    DCM_depth = mean(DCM_depth, na.rm = TRUE),
    max_conc = mean(max_conc, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
                         y = max_conc)) +
  geom_boxplot() +
  geom_point(aes(color = max_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  facet_wrap(~ Year) +  # Create a panel for each year
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  labs(x = "Month", y = "Peak Magnitude", color = "totals ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#visualizing just one box per year

boxplot_Data <- final_DCM_metrics |>
  filter(max_conc > 20) |>
  mutate(DayOfYear = yday(Date))|>
  filter(DayOfYear>133, DayOfYear<286) |>
  mutate(Year = year(Date), Month = month(Date))

ggplot(boxplot_Data, aes(x = factor(Year), y = max_conc)) +
  geom_boxplot() +
  geom_point(aes(color = max_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  ggtitle(label = "Peak Magnitudes only displaying totals > 20")+
  ylim(0, 385) +  # Set the y-axis limits, reversing the range
  labs(x = "Year", y = "Peak Magnitude", color = "totals ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#frame that will be added to for RF analysis at the end
#one measurement for each week that we have data available for 
RF_frame <- final_DCM_metrics|>
  group_by(Year,Week, WaterLevel_m)|>
  summarise(
    DCM_depth = mean(DCM_depth, na.rm = TRUE),
    max_conc = mean(max_conc, na.rm = TRUE),
    totals_mean = mean(totals_mean),
    totals_med = mean(totals_med),
    .groups = "drop"
  )

#metrics for each variable that needs to be calculated 
#depth_var_max
#depth_var_min
#var_max_val
#var_min_val
#var_mean
#var_range (max-main)


#### metals  ####
metalsdf <- metalsdf|>
  mutate(Date = as.Date(DateTime))


variables <- c("SFe_mgL", "TFe_mgL", "SMn_mgL")

metals_summarised <- weekly_sum_variables(metalsdf, variables)

RF_frame_m <- RF_frame|> #random forest frame with metals
  left_join(metals_summarised, by = c("Week", "Year"))
  
#how much raw metal data do we have?This is just to visualize
# #variables <- c(
#   "depth_SFe_mgL_max", "depth_SFe_mgL_min", "SFe_mgL_max_val", "SFe_mgL_min_val", "SFe_mgL_range",
#   "depth_TFe_mgL_max", "depth_TFe_mgL_min", "TFe_mgL_max_val", "TFe_mgL_min_val", "TFe_mgL_range",
#   "depth_SMn_mgL_max", "depth_SMn_mgL_min", "SMn_mgL_max_val", "SMn_mgL_min_val", "SMn_mgL_range")
# data_availability(metals_summarised, variables)
#metals_interpolated <- interpolate_variable(metalsdf_filtered, variables, expanded_dates) no need for this anymore
#data_availability(metals_interpolated, variables)

#### secchi PZ  ####
{
  
  secchi_df <- secchiframe |>
    mutate(Date = as_date(DateTime)) |>
    group_by(Date, Reservoir, Site) |>
    summarise(Secchi_m = mean(Secchi_m, na.rm = TRUE), .groups = "drop") |>
    filter(Reservoir == "BVR" & Site == 50) |>
    mutate(Year = year(Date), DOY = yday(Date))
  
  variables <- c("Secchi_m")
  
  data_availability(secchi_df, variables) #see how much raw secchi data is available
  
  # Ensure DOY_year_ref contains all DOY values for each Year
  DOY_year_ref <- expand.grid(Year = unique(secchi_df$Year), DOY = 1:366) |> # Handle leap years
    filter(!(Year %% 4 != 0 & DOY == 366)) |>  # Remove DOY 366 for non-leap years
    mutate(Date = as.Date(DOY - 1, origin = paste0(Year, "-01-01")))  # Ensure Date matches DOY
  
  # Perform interpolation so that we have secchi for each day
  secchi_interpolated <- DOY_year_ref %>%
    mutate(Week = week(Date))|>
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
    select(Year, Week, Secchi_m)|>
    ungroup()|>
    group_by(Year ,Week)|>
    summarise(Secchi_m = mean(Secchi_m, na.rm = TRUE))
           

# Adding Secchi
RF_frame_ms <- RF_frame_m|> #RF frame w metals and secchi
  left_join(secchi_interpolated, by = c("Year", "Week"))

data_availability(RF_frame_ms, variables)

# Calculating K_d and light availability from secchi

RF_frame_msl <- RF_frame_ms |> #add light
  mutate(sec_K_d = 1.7/Secchi_m) |>
  mutate(PZ = 4.065 /sec_K_d)|>
  group_by(Year, Week)|>
  mutate(PZ = if_else(PZ > WaterLevel_m, WaterLevel_m, PZ))
}

ggplot(RF_frame_msl, aes(x = Date, y = PZ)) +
  geom_line(aes(group = factor(year(Date)))) +
  scale_y_reverse()

####Adding PAR, DO, DOsat_percent, cond, ORP, pH, temp ####
#I don't end up using any of these variables, because we don't have consistent data through all years, 
#but also they are very correlated to other variables that are included
#####YSI#####
ysi_profiles <- ysi_profiles|>
  mutate(Date = as_date(DateTime))

#come back to remove flags
variables <- c("DO_mgL", "PAR_umolm2s", "DOsat_percent", "Cond_uScm", "ORP_mV", "pH", "Temp_C")

data_availability(ysi_profiles, variables)
# Generate the plot
#plot <- data_availability(ysi_profiles, variables)  
# Save the plot with specific dimensions
#ggsave("raw_ysi_availability.png", plot = plot, width = 20, height = 15, dpi = 300)

#removing PAR, ORP, cond, and pH due to limited data availability
#keeping temp because YSI has the most temp

variables <- c("DO_mgL","DOsat_percent", "Temp_C")
ysi <- ysi_profiles|>
  select(-PAR_umolm2s, -ORP_mV, -Cond_uScm, -pH)|>
  filter(Reservoir == "BVR", Site == 50)|>
  filter((hour(DateTime) >= 8), (hour(DateTime) <= 18))
  #remove flags
data_availability(ysi, variables)

#####CTD#####
CTDfiltered <- CTD |>
  filter(Reservoir == "BVR", Site == 50) |>
  filter(!if_any(starts_with("Flag"), ~ . == 68)) |>
  mutate(
    DateTime = ymd_hms(DateTime, tz = "UTC"),
    Date = as_date(DateTime)
  ) |>
  filter(hour(DateTime) >= 8, hour(DateTime) <= 18)
  
variables <- c("DO_mgL", "PAR_umolm2s", "DOsat_percent", "Cond_uScm", "ORP_mV", 
               "pH", "Temp_C")
data_availability(CTDfiltered, variables)

#can't use many of the variables because not enough data for every year
#can use Temp from CTD for 2014, 2015, 2016, 2019, 2022, 2023, and 2024
#can use Temp from YSI for 2017, 2018, 2020

CTDtemp<- CTDfiltered|>
  mutate(Year = year(Date), Week = week(Date))|>
  filter(Year %in% c(2014, 2015, 2016, 2019, 2022, 2023, 2024))|> #remove flags
  select(Date, Year, Week, Temp_C, Depth_m)

ysitemp<- ysi%>%
  mutate(Year = year(Date), Week = week(Date))|>
  filter(Year %in% c(2017, 2018, 2020))|>
  select(Date, Year, Week, Temp_C, Depth_m)

#coalesce 

temp_depths_coalesced <- full_join(ysitemp, CTDtemp, by = c("Date", "Year", "Depth_m", "Week"))|>
  group_by(Date)|>
  mutate(Temp_C = coalesce(Temp_C.x, Temp_C.y))|>
  ungroup()|>
  filter(Depth_m > 0.09)|>
  select(-Temp_C.y, -Temp_C.x)

#variables <- ("Temp_C")
#temp_depths_coalesced_summarised <- weekly_sum_variables(temp_depths_coalesced, variables)

#probably this isn't the most useful information for Temp^^^
#instead I will calculate: VWT, SurfTemp (average of first meter),
#thermocline depth, epilimnion, metalimnion, hypolimnion mean temperatures
#temp range (mean(top 1m)- mean(bottom 1m))
#standard deviation of temp 
#max and min temp
#schmidt stability
#buoyancy freq
####Temp calculations####
temp_depths_cleaned <- temp_depths_coalesced |> #adding buoyancy freq here
  filter(!is.na(Temp_C)) |>
  group_by(Date, Depth_m) |>
  summarise(Temp_C = mean(Temp_C, na.rm = TRUE), .groups = "drop")|>
  mutate(buoyancy_freq = c(buoyancy.freq(Temp_C, Depth_m), NA))#added for padding for the last value
  

temp_weekly_sum <- weekly_sum_variables(temp_depths_cleaned, "Temp_C")

####Thermocline####

# Dataframe with thermocline
just_thermocline <- temp_depths_cleaned |>
  filter(!is.na(Temp_C)) |>
  group_by(Date, Depth_m)|>
  mutate(Temp_C = mean(Temp_C), na.rm = TRUE)|>
  ungroup()|>
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
#    thermocline_depth <- ifelse(thermocline_depth > 0.6 * max_depth, NA, thermocline_depth)
#    thermocline_depth <- ifelse(thermocline_depth < 0.25 * max_depth, NA, thermocline_depth)
    
    # Ensure Date is a single value and return the summarised dataframe
    tibble(Date = .x$Date[1], thermocline_depth = thermocline_depth)
  }) |>
  
  ungroup()|>
  mutate(Week = week(Date), 
         Year = year(Date))

#####individual date thermocline check####
#join thermocline to temp profiles so that I can plot them 
thermocline_and_depth_profiles <- temp_depths_cleaned|>
  left_join(just_thermocline, by = c("Date"))|>
  group_by(Date)|>
  fill(thermocline_depth, .direction = "updown")|>
  ungroup()|>
  mutate(year = year(Date))
  
  
for (yr in years) {
  
  # Filter data for the year
  test <- thermocline_and_depth_profiles |>
    filter(year(Date) == yr)
  
  # Skip if there's no data
  if (nrow(test) == 0) next
  
  # Create plot
  plot_casts <- ggplot(test, aes(x = Temp_C, y = Depth_m)) +
    geom_path() +
    geom_point(size = 0.8, alpha = 0.8) +  # Small points at each observation
    # Light blue grid lines for every whole meter
    geom_hline(yintercept = seq(0, max(test$Depth_m, na.rm = TRUE), by = 1), 
               color = "lightblue", linetype = "dotted", linewidth = 0.3) +
    # Horizontal lines for depths
    geom_hline(aes(yintercept = thermocline_depth), linetype = "dashed", color = "red") +
    # Vertical lines for concentrations
    scale_y_reverse(breaks = seq(0, max(test$Depth_m, na.rm = TRUE), by = 1)) +
    theme_bw() +
    facet_wrap(vars(Date)) +
    xlab("Temp") +
    ylab("Depth (m)") +
    ggtitle(paste(yr, "Temp Profiles"))+
    geom_text(aes(label = round(thermocline_depth, 1), x = Inf, y = thermocline_depth), 
              color = "black", hjust = 1.1, size = 3)
  
  # Save plot
  ggsave(filename = paste0("Figs/Thermocline/", yr, "_raw_casts.png"),
         plot = plot_casts,
         width = 12,
         height = 10,
         dpi = 300)
}

#need to recalculate some of them will come back later for now just add to main dataframe

just_thermocline <- just_thermocline|>
  group_by(Week, Year)|>
  summarise(thermocline_depth = mean(thermocline_depth, na.rm = TRUE))
  
RF_frame_mslt <- RF_frame_msl|>
  left_join(just_thermocline, by = c("Week", "Year"))

####Buoyancy Frequency ####
buoyancy_frame <- temp_depths_cleaned|>
  select(Date, buoyancy_freq, Depth_m)|>
  mutate(Week = week(Date), 
         Year = year(Date))

joined_df <- buoyancy_frame |>
  left_join(RF_frame_mslt |> select(Date, DCM_depth), by = c("Week", "Year"))|>
  filter(!is.na(DCM_depth))


buoyancy_with_dcm <- joined_df |>
  group_by(Week, Year) |>
  mutate(
    depth_diff = abs(Depth_m - DCM_depth),
    N_at_DCM = buoyancy_freq[which.min(depth_diff)]
  ) |>
  ungroup() |>
  select(Week, Year, Depth_m, buoyancy_freq, DCM_depth, N_at_DCM)|>
  select(Week, Year, N_at_DCM)|>
  group_by(Week, Year)|>
  summarise(N_at_DCM = mean(N_at_DCM, na.rm = TRUE))

RF_frame_msltb <- RF_frame_mslt|>
  left_join(buoyancy_with_dcm, by = c("Week", "Year"))

#### Nutrients  ####

chemistry_filtered <- chemistry |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime), 
         DateTime = as.POSIXct(DateTime))|>
  filter((hour(DateTime) >= 8), (hour(DateTime) <= 18))|>
  mutate(
    TN_ugL = if_else(Flag_TN_ugL == 9, NA_real_, TN_ugL),
    TP_ugL = if_else(Flag_TP_ugL == 9, NA_real_, TP_ugL),
    NH4_ugL = if_else(Flag_NH4_ugL == 9, NA_real_, NH4_ugL),
    NO3NO2_ugL = if_else(Flag_NO3NO2_ugL == 9, NA_real_, NO3NO2_ugL),
    DIC_mgL = if_else(Flag_DIC_mgL == 9, NA_real_, DIC_mgL)
  ) |>
  select(Date, Depth_m, TN_ugL, TP_ugL, NH4_ugL, NO3NO2_ugL, DIC_mgL)
 
   
variables <- c("TN_ugL", "TP_ugL", "NH4_ugL", "NO3NO2_ugL", 
               "DIC_mgL")

#raw data availability 
plot <- data_availability(chemistry_filtered, variables)

ggsave("raw_chem_availability.png", plot = plot, width = 20, height = 15, dpi = 300)

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
chemistry_filtered_np <- chemistry_filtered %>%
  mutate(np_ratio = calculate_np_ratio(TN_ugL,TP_ugL))|>
  relocate(np_ratio, .before = TN_ugL)

variables <- c("TN_ugL", "TP_ugL", "NH4_ugL", "NO3NO2_ugL", 
               "DIC_mgL", "np_ratio")

chem_weekly_sum <- weekly_sum_variables(chemistry_filtered_np, variables)

#join to RF frame
RF_frame_msltbc <- RF_frame_msltb|>
  left_join(chem_weekly_sum, by = c("Week", "Year"))

# Visualizing metdata####
# 
# metdata0 <- metdata|>
#   mutate(Date = as_date(DateTime))|>
#   mutate(DOY = yday(Date))|>
#   relocate(DOY, .before = DateTime)|>
#   relocate(Date, .before = DateTime)
# 
# ##### Function for plotting meteorological variables #### 
# metplots <- function(yearz, variable, maxx = NULL){
#   
#   metviz <- metdata0|>
#     filter(year(DateTime) == yearz) #filtering for the year specified
#   
#   
#   ggplot(metviz, aes(x = DateTime, y = {{variable}}))+
#     geom_path()+
#     ggtitle(paste(deparse(substitute(variable)), yearz))+
#     theme_minimal()+
#     scale_y_continuous(limits = c(0, maxx))  # setting consistent y-axis limits
#   
# }
# 
# ##### Precipitation ####
# #not including all of this for now
# metdataprecip <- metdata0 |> 
#   group_by(Date, year(DateTime))|> 
#   mutate(precip_daily = sum(Rain_Total_mm, na.rm = TRUE))|>
#   ungroup()|>
#   relocate(Date, .before = DateTime)|>
#   relocate(precip_daily, .before = Rain_Total_mm)
# 
# #b1 <- metplots(2015, precip_daily, maxx = 80)
# #b2 <- metplots(2016, precip_daily, maxx = 80)
# #b3 <- metplots(2017, precip_daily, maxx = 80)
# #b4 <- metplots(2018, precip_daily, maxx = 80)
# #b5 <- metplots(2019, precip_daily, maxx = 80)
# #b6 <- metplots(2020, precip_daily, maxx = 80)
# #b7 <- metplots(2021, precip_daily, maxx = 80)
# #b8 <- metplots(2022, precip_daily, maxx = 80)
# #b9 <- metplots(2023, precip_daily, maxx = 80)
# 
# #precips<- plot_grid(
# #  b1, b2, b3,
# #  b4, b5, b6, 
# #  b7, b8, b9,
# #  ncol = 3
# #)
# #print(b1)
# 
# #print(precips)
# 
# ##### Air temps #### 
# 
# #b1 <- metplots(2015, AirTemp_C_Average, maxx = 50)
# #b2 <- metplots(2016, AirTemp_C_Average, maxx = 50)
# #b3 <- metplots(2017, AirTemp_C_Average, maxx = 50)
# #b4 <- metplots(2018, AirTemp_C_Average, maxx = 50)
# #b5 <- metplots(2019, AirTemp_C_Average, maxx = 50)
# #b6 <- metplots(2020, AirTemp_C_Average, maxx = 50)
# #b7 <- metplots(2021, AirTemp_C_Average, maxx = 50)
# #b8 <- metplots(2022, AirTemp_C_Average, maxx = 50)
# #b9 <- metplots(2023, AirTemp_C_Average, maxx = 50)
# 
# 
# #temps<- plot_grid(
# #  b1, b2, b3,
# #  b4, b5, b6, 
# #  b7, b8, b9,
# #  ncol = 3
# #)
# 
# #print(temps)
# 
# ##### dailyaverage and dailymax for temps #### 
# metdatatemps <- metdataprecip |> 
#   group_by(Date, year(DateTime))|> 
#   mutate(daily_airtempavg = mean(AirTemp_C_Average, na.rm = TRUE))|>
#   mutate(maxdaily_airtemp = max(AirTemp_C_Average, na.rm = TRUE))|>
#   mutate(mindaily_airtemp = min(AirTemp_C_Average, na.rm = TRUE))|>
#   ungroup()|>
#   relocate(daily_airtempavg, .before = AirTemp_C_Average)|>
#   relocate(maxdaily_airtemp, .before = AirTemp_C_Average)|>
#   relocate(mindaily_airtemp, .before = AirTemp_C_Average)|>
#   select(Date,daily_airtempavg, maxdaily_airtemp, mindaily_airtemp, precip_daily)
# 
# 
# ##### precip and temp to final_datanpratio #### 
# 
# metdata_join <- metdatatemps |> 
#   group_by(Date) |> 
#   summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) |> 
#   distinct()
# 
# final_datamet <- final_datanpratio|>
#   left_join(metdata_join, by = c("Date"), relationship = "many-to-many")
# 
#  
#final dataframe----
Final_RF_frame <- RF_frame_msltbc |>
  mutate(Date = coalesce(Date.x, Date.y))|>
  select(-Date.x, -Date.y, -Secchi_m, -sec_K_d)|>
  relocate(Date, .before = "Year")|>
  rename(buoyancy_freq = N_at_DCM)

#####correlation function----

correlations <- function(year1, year2) {
  DCM_final_cor <- Final_RF_frame |>
    filter(year(Date) >= {{year1}}, year(Date) <= {{year2}}) |>
    filter(month(Date) > 4, month(Date) < 10) |>
    filter(max_conc > 20)
  
  drivers_cor <- cor(DCM_final_cor[,c(4:65)],
                     method = "spearman", use = "pairwise.complete.obs")
 
  list(drivers_cor = drivers_cor, DCM_final_cor = DCM_final_cor)

}

#cutoff 0.7
results <- correlations(2014, 2024)
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
  filter(Variable1 %in% c("DCM_depth"))|>
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

#####correlations across years,max day each year####

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
  filter(Variable1 %in% c("DCM_depth"))|>
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

#####daily correlation, for choosing specific day####

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

####RandomForest Anually####

#"We constructed a RF of 1500 trees for each of the two response
#variables using 1% PAR depth (m), DOC concentration (mg L21),
#thermocline depth (m), metalimnion thickness (m),
#buoyancy frequency at the thermocline (s21),
#lake surface area (log10(km2)), and maximum depth (log10(m))
#as predictors included in each analysis."
#Leach Patterns and Drivers

#prepare data for random forest

library(randomForest)
library(missForest)

set.seed(123)  # Setting seed for reproducibility

# Splitting data into training (70%) and testing (30%)
index <- sample(1:nrow(Final_RF_frame), size = 0.7 * nrow(Final_RF_frame))  # 70% training data
train_data <- Final_RF_frame[index, ]
test_data <- Final_RF_frame[-index, ]

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
  model_rf <- randomForest(DCM_depth ~ ., data = train_data_imputed_z, ntree = 500, importance = TRUE)
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
#  filtered_importance_df <- filtered_importance_df |>
#    filter(!Variable %in% c("peak.top", "pH", "min_pH_depth",  "min_Cond_uScm_depth", "max_Cond_uScm_depth", "peak.magnitude", "peak.bottom", "secchi_PZ", "PAR_PZ", "max_Cond_uScm_depth", "DayOfYear", "DOY", "min_Cond_uScm_dept", "min_CH4_umolL_depth", "max_CH4_umolL_depth"))
  
  
  
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
  
  rmse <- sqrt(mean((test_predictions - train_data_imputed_z$DCM_depth)^2, na.rm = TRUE))
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
  
  



