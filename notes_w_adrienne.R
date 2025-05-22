
metalsdf_filtered <- metalsdf |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  filter(!if_any(starts_with("Flag"), ~. == 68))

variables <- c("SFe_mgL", "TFe_mgL", "TMn_mgL", "SMn_mgL")

#test
test <- interpolate_variable(metalsdf_filtered, variables, expanded_dates)  

looking <- final_result|>
  filter(!is.na(SFe_mgL))

#first clean up dataframes

#pmap purr package
#the argument would be the data frame and in the data
#frame you have a column for each argument in your function that you made

#all_instruments_interpolated <- df_witharguments|>
#purrr::pmap(function)
#then afterwards join all(all_instruments_interpolated)
#look at fcr catwalk in edi already published



#####OTHER NOTES#####
#to do:
#change the way I am interpolating
#met data
#add formula to thermocline function to make sure it doesn't go super deep 

#other variables to add 
#Radiation
#Albedo_Average_W_m2, Calculated from ShortwaveRadiationDown_Average_W_m2 divided by ShortwaveRadiationUp_Average_W_m2	
#higher albedo indicateds greater conductivi;.ty

#ShortwaveRadiationUp_Average_W_m2 (incoming solar radiation), ShortwaveRadiationDown_Average_W_m2 (reflected radiation)
#important though, because this includes PAR, UV poriton (which can have beneficial and harmful effects).
#can drive processes like photosynthesis and photodegradation
#InfraredRadiationUp_Average_W_m2,   InfraredRadiationDown_Average_W_m2


#othermetdata variables to potentially look at
#WindSpeed_Average_m_s, Wind speed averaged over measurement interval	
#WindDir_degrees, Direction of wind at time of measurement




####metalimnion####
#SKIP THIS FOR NOW GO TO BUOYANCY FREQ
metalimnion_prep <- pwmgslyccnt %>%
  group_by(Date, Depth_m) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  group_by(Date) %>%
  distinct(Depth_m, .keep_all = TRUE)|>
  ungroup()

#mutate(Temp_C = if_else(duplicated(Temp_C), Temp_C + runif(n(), min = 0.0001, max = 0.001), Temp_C))
#adjust the temp by .001 if temp duplicated otherwise the metalimnion function will not run

final_datameta <- metalimnion_prep |>
  group_by(Date, ) |>
  arrange(Date, Depth_m) |>
  summarise(
    # Calculate the metalimnion depths based on the temperature
    metalimnion_depths = list(meta.depths(Temp_C, Depth_m, slope = 0.1, seasonal = TRUE, mixed.cutoff = 1)), 
    .groups = 'drop'
  ) |>
  mutate(
    # Extract the upper and lower metalimnion temperatures
    metalimnion_upper= map_dbl(metalimnion_depths, 1),  
    metalimnion_lower = map_dbl(metalimnion_depths, 2)  
  )

#now will join to the final_datathermo
final_datathermometa <- final_datathermo|>
  left_join(final_datameta, by = c("Date", "CastID"))

#visualizing temps at a specific date with the thermocline 
#does the metalimnion upper and lower make sense?
plot_dat<- final_datathermometa|>
  select(Date, CastID, Depth_m, Temp_C, metalimnion_upper, metalimnion_lower, thermocline_depth)|>
  filter(Date %in% c("2014-07-02"))

ggplot(plot_dat, aes(x = Temp_C, y = Depth_m))+
  geom_point()+
  scale_y_reverse()+
  labs(x = "Temp_C", y = "Depth_m")+
  geom_hline(yintercept = unique(plot_dat$thermocline_depth),  # Add horizontal line at thermocline depth
             color = "red",  # Color of the line
             size = 1,       # Line thickness
             linetype = "dashed") +  # Line type (dashed, solid, etc.)
  theme_minimal()

#calculate thickness of metalimnion
alldata_andmetacalc <- final_datathermometa|>
  mutate(metalimnion_upper = if_else(is.na(metalimnion_upper), NA_real_, metalimnion_upper))|>
  mutate(metalimnion_lower = if_else(is.na(metalimnion_lower), NA_real_, metalimnion_lower))|>
  mutate(metalimnion_lower = if_else(thermocline_depth > metalimnion_lower | thermocline_depth < metalimnion_upper, NA_real_, metalimnion_lower), #if thermocline falls outside the indicated metalimnion, the metalimnion is calculated incorrectly
         metalimnion_upper = if_else(thermocline_depth > metalimnion_lower | thermocline_depth < metalimnion_upper, NA_real_, metalimnion_upper))|>
  mutate(meta_width = (metalimnion_lower-metalimnion_upper))|>
  select(-Date.y, -na.rm, -thermocline_depth.y, -thermocline_depth.x, -Date.x.x, -Date.x)|>
  mutate(meta_width = if_else(is.na(thermocline_depth), NA_real_, meta_width))
#SKIP THIS FOR NOW GO TO BUOYANCY FREQ
#need to add the bathymetry here for surface area at different depths

BVRbath <- bath|>
  filter(Reservoir == "BVR")

library(signal)

new_depths <- seq(0, 14, by = 0.01)
interpolated_SA <- signal::pchip(BVRbath$Depth_m, BVRbath$SA_m2, new_depths)
interpolated_Volume_layer <- signal::pchip(BVRbath$Depth_m, BVRbath$Volume_layer_L, new_depths)
interpolated_Volume_below <- signal::pchip(BVRbath$Depth_m, BVRbath$Volume_below_L, new_depths)

#new bathymetry dataframe with finer sequence
BVRbath_interpolated <- data.frame(
  Depth_m = new_depths,
  SA_m2 = interpolated_SA,
  Volume_layer_L = interpolated_Volume_layer,
  Volume_below_L = interpolated_Volume_below
)

BVRbath_interpolated<- BVRbath_interpolated|>
  filter(SA_m2 != 0, Depth_m != 0)|>
  filter(Depth_m >= 13.4)

bathytest <- final_data_water|>
  group_by(Date)|>
  mutate(Dadjust = 13.4-WaterLevel_m)|> 
  mutate(tempbathdepths = Depth_m + Dadjust)|> #I will use this depth to extract the surface area from BVRbath_interpolated
  ungroup()

final_bathy <- bathytest |>
  mutate(
    SA_m2 = approx(BVRbath_interpolated$Depth_m, BVRbath_interpolated$SA_m2, tempbathdepths, rule = 2)$y,
    Volume_layer_L = approx(BVRbath_interpolated$Depth_m, BVRbath_interpolated$Volume_layer_L, tempbathdepths, rule = 2)$y,
    Volume_below_L = approx(BVRbath_interpolated$Depth_m, BVRbath_interpolated$Volume_below_L, tempbathdepths, rule = 2)$y
  )

question<- final_bathy|>
  select(Date, Depth_m, TotalConc_ugL, tempbathdepths, WaterLevel_m, Dadjust, SA_m2, Volume_layer_L, Volume_below_L)|>
  group_by(Date)|>
  mutate(difference = max(Depth_m)- WaterLevel_m)

ggplot(question, aes(x = Date, y = difference))+
  geom_point()

#look at BVRbath for bathymetry comparisons

#SKIP THIS FOR NOW GO TO BUOYANCY FREQ
####whole lake temp####
#Calculates volumetrically weighted average whole lake temperature using the supplied water temperature timeseries.

#use tempbathdepths when using packages that require bathymetric data. adjusted to match up the 0-14 bathymetric data
#final_bathy <- final_bathy |>
#  select(-CastID, -DateTime)|>
#  group_by(Date, Depth_m) |>
#  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop')

lake_temp <- final_bathy %>%
  filter(!is.na(Temp_C))|>
  group_by(Date) %>%
  mutate(whole_lake_temp = whole.lake.temperature(Temp_C, tempbathdepths, BVRbath_interpolated$Depth_m, BVRbath_interpolated$SA_m2))|>
  ungroup()

looking<- final_bathy|>
  select(Date, CastID, Depth_m, Temp_C, tempbathdepths, WaterLevel_m)

#I think i did it need to check on this 

looking<- lake_temp|>
  select(Date, Depth_m, whole_lake_temp)
