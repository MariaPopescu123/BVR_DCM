# Maria DCM BVR
#data includes chlorophyll maxima, PAR, secchi, light attenuation, metals, ghgs, nutrients


beep <- function(){
  system("rundll32 user32.dll,MessageBeep") #just putting this here so I can be notified when something that takes a while to run is done
}  

pacman::p_load(tidyverse, lubridate, akima, reshape2, 
               gridExtra, grid, colorRamps, RColorBrewer, rLakeAnalyzer,
               reader, cowplot, dplyr, tidyr, ggplot2, zoo, purrr, beepr, forecast, ggthemes)

#### Loading Data  ####

#CTD
CTD <- read.csv("./CTD.csv")

#flora
current_df <- read.csv("./current_df.csv")

# metals 
metalsdf <- read.csv("./metalsdf.csv")
#removed flags for 68 as per Cece's advice

#ghgss
ghgs <- read.csv("./ghgs.csv")
#not sure whether or not to remove those with flags 3 and 4
#3 = The difference between the reps are above the limit of quantification and >30% and <50% different from each other. Both replicates were retained but flagged
#4 = The difference between the reps are above the limit of quantification and >50% different from each other. Both replicates were retained but flagged

#secchi
secchiframe <- read.csv("./secchiframe.csv")

#ysi
ysi_profiles <- read.csv("ysi_profiles.csv")

#chemistry
chemistry <- read.csv("./chemistry.csv")

#meteorological data from FCR https://portal.edirepository.org/nis/mapbrowse?packageid=edi.389.8
options(timeout = 300)
metdata <- read.csv("./metdata.csv")

#bathymetry
bath <- read.csv("./bath.csv")

#waterlevel data
wtrlvl <- read.csv("./wtrlvl.csv") 

#waterlevel platform data (for past 2020)
BVRplatform <- read.csv("./BVRplatform.csv")

#adding columns with total_conc max and the depth at which it occurs
phytos <- current_df %>% 
  filter(Reservoir == "BVR", Site == 50)%>%
  mutate(Date  = as_date(DateTime)) |> 
  filter((hour(DateTime) >= 8), (hour(DateTime) <= 15))|>
  filter(!(CastID == 592))|> #filter out weird drop in 2017
  filter(!(CastID == 395))|> #weird drop in 2016
#  filter(Flag_TotalConc_ugL != 2,Flag_TotalConc_ugL != 3)|> #2 is instrument malfunction and #3 is low transmission value
  group_by(CastID)|>
  mutate(Totals_DCM_conc = max(TotalConc_ugL, na.rm = TRUE))|> #concentration of totals at totals DCM
  mutate(Totals_DCM_depth = ifelse(TotalConc_ugL == Totals_DCM_conc, Depth_m, NA_real_))|>
  ungroup()

DCM_BVRdata <- phytos %>%
  group_by(CastID) %>%
  fill(Totals_DCM_conc, .direction = "downup")|>
  fill(Totals_DCM_depth, .direction = "downup")|>
  ungroup()%>%
  mutate(DOY = yday(DateTime))%>%
  mutate(Date  = as_date(DateTime)) |> 
  select(Date, DateTime, CastID, Depth_m, Totals_DCM_conc, Totals_DCM_depth, TotalConc_ugL,  Temp_C)

#### metals  ####
{
metalsdf_filtered <- metalsdf |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  filter(!if_any(starts_with("Flag"), ~. == 68))|>
  group_by(Depth_m, Date)|>
  summarise(SFe_mgL = mean(SFe_mgL, na.rm = TRUE),
            TFe_mgL = mean(TFe_mgL, na.rm = TRUE),
            SMn_mgL = mean(SMn_mgL, na.rm = TRUE), 
            SCa_mgL = mean(SCa_mgL, na.rm = TRUE),
            TCa_mgL = mean(TCa_mgL, na.rm = TRUE),
            TCu_mgL = mean(TCu_mgL, na.rm = TRUE),
            SCu_mgL = mean(SCu_mgL, na.rm = TRUE),
            SBa_mgL = mean(SBa_mgL, na.rm = TRUE), 
            TBa_mgL = mean(TBa_mgL, na.rm = TRUE))


interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(metalsdf_filtered) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(interp_SFe_mgL = zoo::na.approx(SFe_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_TFe_mgL = zoo::na.approx(TFe_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_SMn_mgL = zoo::na.approx(SMn_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_SCa_mgL = zoo::na.approx(SCa_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_TCa_mgL = zoo::na.approx(TCa_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_TCu_mgL = zoo::na.approx(TCu_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_SCu_mgL = zoo::na.approx(SCu_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_SBa_mgL = zoo::na.approx(SBa_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_TBa_mgL = zoo::na.approx(TBa_mgL, x = Depth_m, na.rm = FALSE)) |>
  ungroup() # Ensure the grouping is removed for the final merge

#Merge with DCM BVR data and keep only relevant rows
DCM_BVRwmetals <- DCM_BVRdata %>%
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  select(-SFe_mgL, -TFe_mgL, -SMn_mgL, -SCa_mgL, -TCa_mgL, -TCu_mgL, -SCu_mgL, -SBa_mgL, -TBa_mgL) %>% # Remove unnecessary columns
  filter(Depth_m %in% DCM_BVRdata$Depth_m) # Keep only rows with depths present in DCMdata
}

#### ghgs  ####
{
ghgs_filtered <- ghgs |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(CO2_umolL = mean(CO2_umolL, na.rm = TRUE)
            , CH4_umolL = mean(CH4_umolL, na.rm = TRUE))

interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(ghgs_filtered) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(interp_CO2_umolL = zoo::na.approx(CO2_umolL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_CH4_umolL = zoo::na.approx(CH4_umolL, x = Depth_m, na.rm = FALSE)) |>
  ungroup() # Ensure the grouping is removed for the final merge

DCM_BVRwmetalsghgs <- DCM_BVRwmetals |>
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many")|>
  select(-CO2_umolL, -CH4_umolL)|>
  filter(Depth_m %in% DCM_BVRwmetals$Depth_m)|> #filtering to make sure only the values with the depths in fluoroprobe data are brought in
  mutate(DCM = (Totals_DCM_depth==Depth_m))|>
  relocate(DCM, .before = 8)
}

#### secchi and attenuation coefficient  ####
{
secchi <- secchiframe |>
  mutate(Date = as_date(DateTime)) |>
  group_by(Date, Reservoir, Site) |>
  summarise(Secchi_m = mean(Secchi_m, na.rm = TRUE))

BVRsecchi <- secchi |>
  filter(Reservoir == "BVR" & Site == 50)

# Adding Secchi
DCM_BVRwmetalsghgssecchi <- DCM_BVRwmetalsghgs|>
  left_join(BVRsecchi, by = c("Date"))|>
  group_by(Date)|>
  fill(Secchi_m, .direction = "updown")|>
  ungroup()

# Calculating K_d and light availability from secchi

DCM_BVRwmetalsghgssecchilight <- DCM_BVRwmetalsghgssecchi |>
  mutate(sec_K_d = 1.7/Secchi_m) |>
  mutate(light_availability_fraction = exp(-sec_K_d * Depth_m)) |>
  mutate(sec_LAP = light_availability_fraction * 100) #light availability percentage calculated from secchi
}

####Adding PAR, DO, DOsat_percent, cond, ORP, pH, temp ####

ysi_profiles_filtered <- ysi_profiles |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(DO_mgL = mean(DO_mgL, na.rm = TRUE),
            PAR_umolm2s = mean(PAR_umolm2s, na.rm = TRUE), 
            DOsat_percent = mean(DOsat_percent, na.rm = TRUE),
            Cond_uScm = mean(Cond_uScm, na.rm = TRUE),
            ORP_mV = mean(ORP_mV, na.rm = TRUE),
            pH = mean(pH, na.rm = TRUE),
            Temp_C = mean(Temp_C, na.rm = TRUE))

CTDfiltered <- CTD|> #flag 2, instrument malfunction. haven't removed flags yet
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(DO_mgL = mean(DO_mgL, na.rm = TRUE),
            PAR_umolm2s = mean(PAR_umolm2s, na.rm = TRUE), 
            DOsat_percent = mean(DOsat_percent, na.rm = TRUE),
            Cond_uScm = mean(Cond_uScm, na.rm = TRUE),
            ORP_mV = mean(ORP_mV, na.rm = TRUE),
            pH = mean(pH, na.rm = TRUE),
            Temp_C = mean(Temp_C))

#now join both together
combined_df <- full_join(ysi_profiles_filtered, CTDfiltered, by = c("Date", "Depth_m"))

combined_df2 <- combined_df %>%
  filter(year(Date) != 2013) %>%
  mutate(PAR_umolm2s = coalesce(PAR_umolm2s.y, PAR_umolm2s.x)) %>%
  mutate(DO_mgL = coalesce(DO_mgL.y, DO_mgL.x)) %>%
  mutate(DOsat_percent = coalesce(DOsat_percent.y, DOsat_percent.x)) %>%
  mutate(Cond_uScm = coalesce(Cond_uScm.y, Cond_uScm.x)) %>%
  mutate(ORP_mV = coalesce(ORP_mV.y, ORP_mV.x)) %>%
  mutate(pH = coalesce(pH.y, pH.x)) %>%
  mutate(Temp_C = coalesce(Temp_C.y, Temp_C.x)) %>%
  select(-PAR_umolm2s.x, -PAR_umolm2s.y,
         -DO_mgL.x, -DO_mgL.y,
         -DOsat_percent.x, -DOsat_percent.y,
         -Cond_uScm.x, -Cond_uScm.y,
         -ORP_mV.x, -ORP_mV.y,
         -pH.x, -pH.y,
         -Temp_C.x, -Temp_C.y)

#this was used from Lewis et al. 2023
{
atten_k <- combined_df2%>%
  filter(!is.na(PAR_umolm2s),
         !is.na(Depth_m))%>%
  group_by(Date)%>%
  filter(sum(Depth_m<0)>0)%>%
  mutate(I0 = mean(PAR_umolm2s[Depth_m<0], na.rm = T),
         PAR_umolm2s = ifelse(PAR_umolm2s==0,0.001,PAR_umolm2s))%>%
  filter(Depth_m>0,
         !I0==0)%>%
  summarize(I0 = unique(I0),
            k = coef(lm(I(log(PAR_umolm2s)-log(I0))~ 0 + Depth_m)),
            r2 = summary(lm(I(log(PAR_umolm2s)-log(I0)) ~ 0 + Depth_m))$r.squared,
            Zeu = min(Depth_m[PAR_umolm2s<I0/100]),
            Zeu_0.1 = min(Depth_m[PAR_umolm2s<I0/1000]))%>%
  filter(r2>0.9)|>
  select(Date, Zeu, Zeu_0.1)
}#ends here


interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(combined_df2) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(interp_PAR_umolm2s = zoo::na.approx(PAR_umolm2s, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_DO_mgL = zoo::na.approx(DO_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_DOsat_percent = zoo::na.approx(DOsat_percent, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_Cond_uScm = zoo::na.approx(Cond_uScm, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_ORP_mV = zoo::na.approx(ORP_mV, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_pH = zoo::na.approx(pH, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_Temp_C = zoo::na.approx(Temp_C, x = Depth_m, na.rm = FALSE)) |>
  ungroup()|>
  select(-DO_mgL, -PAR_umolm2s, -DOsat_percent, -Cond_uScm, -ORP_mV, -pH, -Temp_C)

final_PAR <- DCM_BVRwmetalsghgssecchilight %>%
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  filter(Depth_m %in% DCM_BVRdata$Depth_m) # Keep only rows with depths present in flora data

final_PAR <- final_PAR|> #joining lewis calculations
  left_join(atten_k, by = "Date", relationship = "many-to-many")

library(conflicted)
conflicts_prefer(dplyr::lag)

PAR_Kd_PZ <- final_PAR |>
  group_by(CastID) |>
  mutate(interp_PAR_umolm2s = if_else(interp_PAR_umolm2s == 0, 0.001, interp_PAR_umolm2s))|>
  mutate(PAR_kd = (log(lag(interp_PAR_umolm2s)) - log(interp_PAR_umolm2s)) / (Depth_m - lag(Depth_m)),
         PAR_kd = ifelse(is.na(PAR_kd) | Depth_m == lag(Depth_m), NA, PAR_kd)) |>
  mutate(mean_Kd = mean(PAR_kd, na.rm = TRUE)) |>
  mutate(PAR_PZ = ifelse(is.na(mean_Kd) | mean_Kd == 0, NA, log(100) / mean_Kd))

# Now calculating light availability percentage (LAP) using PAR
final_dataLAP <- PAR_Kd_PZ |>
  group_by(CastID) |>
  mutate(PAR_LAP = 100 * exp(-PAR_kd * Depth_m))


####Secchi PZ####
#using secchi_PZ because the data for PAR between CTD and YSI is too different (for example look at ysi_profiles filtered and CTD filtered for 2018-5-4)
#if I want to calculate PZ for specific years it would be ok but across all years no
final_datasecPZ <- final_dataLAP |>
  mutate(secchi_PZ = 2.7*Secchi_m)|>
  relocate(secchi_PZ, .before = sec_LAP)|>
  relocate(PAR_LAP, .after = sec_LAP)|>
  group_by(Date, Depth_m, CastID)|>
  mutate(Temp_C = rowMeans(cbind(Temp_C, interp_Temp_C), na.rm = TRUE))|>
  select(-interp_Temp_C)


#decided to use secchi based on exploring data availability and values across years)


#### Nutrients  ####
{
chemistry_filtered <- chemistry |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Depth_m, Date)|>
  summarise(TN_ugL = mean(TN_ugL, na.rm = TRUE),
            TP_ugL = mean(TP_ugL, na.rm = TRUE),
            NH4_ugL = mean(NH4_ugL, na.rm = TRUE), 
            NO3NO2_ugL = mean(NO3NO2_ugL, na.rm = TRUE),
            SRP_ugL = mean(SRP_ugL, na.rm = TRUE),
            DOC_mgL = mean(DOC_mgL, na.rm = TRUE),
            DIC_mgL = mean(DIC_mgL, na.rm = TRUE),
            DC_mgL = mean(DC_mgL, na.rm = TRUE))

interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(chemistry_filtered) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(interp_TN_ugL = zoo::na.approx(TN_ugL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_TP_ugL = zoo::na.approx(TP_ugL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_NH4_ugL = zoo::na.approx(NH4_ugL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_NO3NO2_ugL = zoo::na.approx(NO3NO2_ugL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_SRP_ugL = zoo::na.approx(SRP_ugL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_DOC_mgL = zoo::na.approx(DOC_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_DIC_mgL = zoo::na.approx(DIC_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_DC_mgL = zoo::na.approx(DC_mgL, x = Depth_m, na.rm = FALSE)) |>
  ungroup() # Ensure the grouping is removed for the final merge

#Merge with DCM BVR data and keep only relevant rows
final_datanutrients <- final_datasecPZ %>%
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  select(-TN_ugL,-TP_ugL, -NH4_ugL, -NO3NO2_ugL,-SRP_ugL, -DOC_mgL, -DIC_mgL, -DC_mgL) %>% # Remove unnecessary columns
  filter(Depth_m %in% DCM_BVRdata$Depth_m)|> # Keep only rows with depths present in DCMdata
  mutate(Reservoir = "BVR")|>
  mutate(Site = 50)|>
  relocate(Reservoir, .before = 1)|>
  relocate(Site, .before = 1)


conflicts_prefer(dplyr::filter)

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

conflicts_prefer(dplyr::filter)
# added np ratio to dataframe
final_datanpratio <- final_datanutrients %>%
  mutate(np_ratio = calculate_np_ratio(interp_TN_ugL,interp_TP_ugL))|>
  relocate(np_ratio, .before = interp_TN_ugL)
}

####metdata####
metdata0 <- metdata|>
  mutate(Date = as_date(DateTime))|>
  mutate(DOY = yday(Date))|>
  relocate(DOY, .before = DateTime)|>
  relocate(Date, .before = DateTime)

##### Precipitation #####
metdataprecip <- metdata0 |> 
  group_by(Date, year(DateTime))|> 
  mutate(precip_daily = sum(Rain_Total_mm, na.rm = TRUE))|>
  ungroup()|>
  relocate(Date, .before = DateTime)|>
  relocate(precip_daily, .before = Rain_Total_mm)

##### dailyaverage and dailymax for temps #####
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

####viz met variables ####
#air_tempavg
#maxdaily_airtemp
#mindaily_airtemp
#precip_daily

metplots <- function(yearz, variable, maxx = NULL){
  
  metviz <- metdatatemps|>
    filter(year(Date) == yearz) #filtering for the year specified
  
  
  ggplot(metviz, aes(x = Date, y = {{variable}}))+
    geom_path()+
    ggtitle(paste(deparse(substitute(variable)), yearz))+
    theme_minimal()+
    scale_y_continuous(limits = c(0, maxx))  # setting consistent y-axis limits
  
}



#####viz Air temps ##### 

b1 <- metplots(2015, daily_airtempavg, maxx = 50)
b2 <- metplots(2016, daily_airtempavg, maxx = 50)
b3 <- metplots(2017, daily_airtempavg, maxx = 50)
b4 <- metplots(2018, daily_airtempavg, maxx = 50)
b5 <- metplots(2019, daily_airtempavg, maxx = 50)
b6 <- metplots(2020, daily_airtempavg, maxx = 50)
b7 <- metplots(2021, daily_airtempavg, maxx = 50)
b8 <- metplots(2022, daily_airtempavg, maxx = 50)
b9 <- metplots(2023, daily_airtempavg, maxx = 50)


temps<- plot_grid(
  b1, b2, b3,
  b4, b5, b6, 
  b7, b8, b9,
  ncol = 3
)

print(temps)


#max daily temp
b1 <- metplots(2015, maxdaily_airtemp, maxx = 50)
b2 <- metplots(2016, maxdaily_airtemp, maxx = 50)
b3 <- metplots(2017, maxdaily_airtemp, maxx = 50)
b4 <- metplots(2018, maxdaily_airtemp, maxx = 50)
b5 <- metplots(2019, maxdaily_airtemp, maxx = 50)
b6 <- metplots(2020, maxdaily_airtemp, maxx = 50)
b7 <- metplots(2021, maxdaily_airtemp, maxx = 50)
b8 <- metplots(2022, maxdaily_airtemp, maxx = 50)
b9 <- metplots(2023, maxdaily_airtemp, maxx = 50)


temps<- plot_grid(
  b1, b2, b3,
  b4, b5, b6, 
  b7, b8, b9,
  ncol = 3
)

print(temps)



#####viz precip#####
b1 <- metplots(2015, precip_daily, maxx = 80)
b2 <- metplots(2016, precip_daily, maxx = 80)
b3 <- metplots(2017, precip_daily, maxx = 80)
b4 <- metplots(2018, precip_daily, maxx = 80)
b5 <- metplots(2019, precip_daily, maxx = 80)
b6 <- metplots(2020, precip_daily, maxx = 80)
b7 <- metplots(2021, precip_daily, maxx = 80)
b8 <- metplots(2022, precip_daily, maxx = 80)
b9 <- metplots(2023, precip_daily, maxx = 80)


precips<- plot_grid(
  b1, b2, b3,
  b4, b5, b6, 
  b7, b8, b9,
  ncol = 3
)

print(precips)



#### merge metdata to final_datanpratio #### 

metdata_join <- metdatatemps |> 
  group_by(Date) |> 
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) |> 
  distinct()

final_datamet <- final_datanpratio|>
  left_join(metdata_join, by = c("Date"), relationship = "many-to-many")

  
####thermocline ####

# Dataframe with thermocline
thermocline_df <- final_datamet |>
  group_by(CastID, Date, Depth_m) |>
  summarize(Temp_C = mean(Temp_C, na.rm = TRUE)) |>
  ungroup() |>
  group_by(CastID) |>
  mutate(thermocline_depth = thermo.depth(
    Temp_C, 
    Depth_m, 
    Smin = 2, 
    seasonal = TRUE, 
    index = FALSE,
    mixed.cutoff = 3
  )) |>
  ungroup()

#fixing incorrect thermoclines
thermocline_df_unique<- thermocline_df|>
  filter(thermocline_depth<3)|>
  group_by(Date)|>
  filter(Depth_m > thermocline_depth)|>
  mutate(thermocline_depth = thermo.depth(Temp_C, 
                                          Depth_m, 
                                          Smin = 2, 
                                          seasonal = TRUE, 
                                          index = FALSE,
                                          mixed.cutoff = 3
  ))|>
  ungroup()

both_merged<- thermocline_df|>
  left_join(thermocline_df_unique, by = c("CastID", "Depth_m", "Temp_C"), relationship = "many-to-many")|>
  mutate(thermocline_depth = coalesce(thermocline_depth.y, thermocline_depth.x))|>
  mutate(Date = Date.x)|>
  select(-Date.y, -thermocline_depth.y, -thermocline_depth.x)

#add to frame
final_datathermocline <- final_datamet|>
  left_join(both_merged, by = c("CastID", "Depth_m", "Temp_C"), relationship = "many-to-many")

initial_thermo <- final_datathermocline|>
  group_by(CastID)|>
  fill(thermocline_depth, .direction = "updown")|>
  ungroup()|>
  relocate(thermocline_depth, .before = Temp_C)

final_datathermo <- initial_thermo|>
  mutate(Date = Date.x.x)|>
  #thermocline added manually via plot inspection
  mutate(thermocline_depth = if_else(Date %in% c("2022-06-27", "2022-08-08"), 3, thermocline_depth))|>
  mutate(thermocline_depth = if_else(!Date %in% c("2016-05-19", "2016-08-23","2021-08-09", "2019-09-04"),thermocline_depth, 5.2))|>
  mutate(thermocline_depth = if_else(!Date %in% c("2018-06-21", "2016-10-25"), thermocline_depth, NA_real_))|> #watercolumn mixed
  mutate(thermocline_depth = if_else(Date %in% c("2020-08-06"), 3.75, thermocline_depth))|>
  group_by(Date)|>
  mutate(thermocline_depth = mean(thermocline_depth), na.rm = TRUE)|>
  select(-Date.x, -Date.y, -Date.x.x)|>
  relocate(Date, .before = DateTime)

#checking to make sure the thermocline is where I expect
conflicts_prefer(dplyr::filter)

plot_data <- final_datathermo |>
  filter(Date %in% c("2019-09-04"))|> #change depths here to see a specific day and see if the thermocline matches up
  select(Date, Depth_m, thermocline_depth, Temp_C)
# Extract the thermocline depth for the specific date for the line
thermocline_depth_value <- unique(plot_data$thermocline_depth)

# Create the plot
ggplot(plot_data, aes(x = Temp_C, y = Depth_m)) +
  geom_point() +  # Add points
  geom_hline(yintercept = thermocline_depth_value, linetype = "dashed", color = "red") +  # Add the thermocline line
  scale_y_reverse() +  # Inverts the y-axis
  labs(x = "Temperature (°C)", y = "Depth (m)", title = "Thermocline Depth on 2016-10-25") +
  theme_minimal()  # Optional: apply a clean theme


#visualize thermocline depth across all years
ggplot(final_datathermo, aes(x = Date, y = thermocline_depth))+
  scale_y_reverse()+
  geom_point()
#2016-10-11 looks very deep but it is real



#very deep thermoclines not fixed yet
#"2021-09-06" recalculate this thermocline. currently at 6.75 and it looks closer to 5

#2019-09-04


####metalimnion####

final_datametaprep <- final_datathermo %>%
  group_by(Date, CastID, Depth_m) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = 'drop') %>%
  group_by(Date, CastID) %>%
  distinct(Depth_m, .keep_all = TRUE)|>
  ungroup()

final_datatempadjust <- final_datametaprep %>%
  group_by(Date, CastID) %>%
  mutate(Temp_C = if_else(duplicated(Temp_C), Temp_C + runif(n(), min = 0.0001, max = 0.001), Temp_C))
#adjust the temp by .001 if temp duplicated otherwise the metalimnion function will not run

final_datameta <- final_datatempadjust |>
  group_by(Date, CastID) |>
  arrange(Date, CastID, Depth_m) |>
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
  select(-na.rm)|>
  mutate(meta_width = if_else(is.na(thermocline_depth), NA_real_, meta_width))


####Buoyancy Frequency ####

final_databuoy <- alldata_andmetacalc|>
  group_by(CastID)|>
  mutate(buoyancy_freq = c(buoyancy.freq(Temp_C, Depth_m), NA))|>#added for padding for the last value
  relocate(buoyancy_freq, .before = thermocline_depth)
#need to make sure this makes sense

####Waterlevels####

wtrlvl2<- wtrlvl|>
  mutate(Date = as.Date(DateTime))|>
  filter(Reservoir == "BVR", Site == 50)

#list of DOY for interpolation purpose
DOY_list <- 32:334  # DOYs from February 1 to November 30
years <- unique(year(wtrlvl2$Date))

DOY_year_ref <- expand.grid(Year = years, DOY = DOY_list)|>
  arrange(Year, DOY)

#Add DOY and Year columns to wtrlvl2, then join with DOY_year_ref
wtrlvl2 <- wtrlvl2 |>
  mutate(Year = year(Date), DOY = yday(Date))

#join and interpolate WaterLevel_m for each DOY in each year
wtrlvl2_interpolated <- DOY_year_ref |>
  left_join(wtrlvl2, by = c("Year" = "Year", "DOY" = "DOY")) |>
  group_by(Year) |>
  mutate(
    WaterLevel_m = na.approx(WaterLevel_m, x = DOY, na.rm = FALSE)
  )|>
  filter(Year > 2013)|>
  arrange(Year, DOY)

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
    LvlDepth_m_13 = na.approx(LvlDepth_m_13, x = DOY, na.rm = FALSE)
  )|>
  filter(Year > 2019, Site == 50)|>
  arrange(Year, DOY)|>
  select(Year, DOY, DateTime, LvlDepth_m_13)

water_levelsjoined <- DOY_year_ref|>
  left_join(BVRplatform2_interpolated, by = c("Year", "DOY"), relationship = "many-to-many")|>
  left_join(wtrlvl2_interpolated, by = c("Year", "DOY"), relationship = "many-to-many")|>
  filter(Year>2013)

water_levelscoalesced<- water_levelsjoined|>
  mutate(WaterLevel_m = coalesce(LvlDepth_m_13,WaterLevel_m))|>
  select(Year, DOY, WaterLevel_m)

final_data_water<- final_databuoy|>
  mutate(DOY = yday(Date))|>
  mutate(Year = year(Date))|>
  left_join(water_levelscoalesced, by = c("DOY", "Year"), relationship = "many-to-many" )

#ggplot(final_data_water, aes(x = Date, y = WaterLevel_m))+
#  geom_point()+
#  labs(title = "BVR Water Level")+
#  scale_x_date(date_labels = "%Y")


#SKIP THIS FOR NOW. GO TO PEAK WIDTH
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

#SKIP THIS FOR NOW GO TO PEAK
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

####Peak.width####
#focusing on totals

#separate data frame for peak widths, depths, and magnitude calculations
for_peaks <- final_data_water|> #type in here the last frame that it matches up with (coming back to this because can't get whole lake to work)
  select(-Site, -Reservoir, -DateTime, -CastID, -DCM)|>
  group_by(Date, Depth_m) |>
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

peaks_calculated <- for_peaks %>%
  group_by(Date) %>%
  mutate(
    total_med = median(TotalConc_ugL, na.rm = TRUE),  # Calculate the median, excluding NA values
    total_sd = sd(TotalConc_ugL, na.rm = TRUE),       # Calculate the standard deviation
    total_mean = mean(TotalConc_ugL, na.rm = TRUE),   # Calculate the mean
    total_mean_plus_sd = total_mean + total_sd,          # Calculate mean + sd
    peak.top = as.integer(Depth_m <= Totals_DCM_depth & TotalConc_ugL > total_mean_plus_sd),  # Create binary indicator
    peak.bottom = as.integer(Depth_m >= Totals_DCM_depth & TotalConc_ugL > total_mean_plus_sd),
    
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

final_data_peaks <- peaks_calculated|>
  group_by(Date)|>
  mutate(peak.magnitude = max(TotalConc_ugL)-mean(TotalConc_ugL))|>
  ungroup()|>
  select(Date, Depth_m, total_mean, total_sd, total_mean_plus_sd, peak.top, peak.bottom, peak.width, peak.magnitude) #this is unnecessary. saying how many there are at the DCM for total_conc

library(lubridate)
conflicts_prefer(dplyr::filter)
library(dplyr)

final_data0 <- final_data_water |>
  left_join(final_data_peaks, by = c("Date", "Depth_m")) |>
  mutate(peak.width = if_else(peak.width < 3, peak.width, NA_real_)) |>
  group_by(Date, CastID, Depth_m) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), 
            .groups = 'drop')|>
  mutate(DayOfYear = yday(Date))|>
  mutate(DCM = if_else(Depth_m == Totals_DCM_depth, TRUE, FALSE))|>
  mutate(PAR_PZ = if_else(PAR_PZ<0, NA_real_, PAR_PZ))|>
  mutate(PZ = if_else(!is.na(secchi_PZ), 
                      secchi_PZ, 
                      rowMeans(select(cur_data(), PAR_PZ, Zeu), na.rm = TRUE)))|>
  mutate(PZ = if_else(PZ>10, 9.5, PZ))|>
  mutate(secchi_PZ = if_else(secchi_PZ>10, 9.5, secchi_PZ))|>
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) |>
  mutate(DayOfYear = yday(Date)) |>
  filter(DayOfYear > 133, DayOfYear < 286)|>  # Timeframe filtering     May 13 - October 12
  rename_with(~ gsub("^interp_", "", .), starts_with("interp_"))  # Remove "interp_" prefix

#same thing but including all data to see flora data availability
final_data_alldates <- final_data_water|>
  left_join(final_data_peaks, by = c("Date", "Depth_m")) |>
  mutate(peak.width = if_else(peak.width < 3, peak.width, NA_real_)) |>
  group_by(Date, CastID, Depth_m) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), 
            .groups = 'drop')|>
  mutate(DayOfYear = yday(Date))|>
  mutate(DCM = if_else(Depth_m == Totals_DCM_depth, TRUE, FALSE))|>
  mutate(PAR_PZ = if_else(PAR_PZ<0, NA_real_, PAR_PZ))|>
  mutate(PZ = if_else(!is.na(secchi_PZ), 
                      secchi_PZ, 
                      rowMeans(select(cur_data(), PAR_PZ, Zeu), na.rm = TRUE)))|>
  mutate(PZ = if_else(PZ>10, 9.5, PZ))|>
  mutate(secchi_PZ = if_else(secchi_PZ>10, 9.5, secchi_PZ))|>
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) |>
  mutate(DayOfYear = yday(Date)) |>
  rename_with(~ gsub("^interp_", "", .), starts_with("interp_"))  # Remove "interp_" prefix


####Schmidt_stability####

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



####final dataframe####

write.csv(final_data0,"./final_data0.csv",row.names = FALSE)


write.csv(final_data_alldates, "./final_data_alldates.csv", row.names = FALSE)


looking<- final_data0|>
  filter(TotalConc_ugL>380)

#"2022-08-18"

####DCM depth correlations####
#removed buoyancy_freq for now bc had -inf will come back to

# Create a vector of variable names that need to be summarized
depth_variables <- c("Temp_C", "np_ratio", "SFe_mgL", "TFe_mgL", 
                     "SMn_mgL", "SCa_mgL", "TCa_mgL", 
                     "TCu_mgL", "SBa_mgL", "TBa_mgL", 
                     "CO2_umolL", "CH4_umolL", "DO_mgL", 
                     "DOsat_percent", "Cond_uScm", "ORP_mV", 
                     "pH", "TN_ugL", "TP_ugL", 
                     "NH4_ugL", "NO3NO2_ugL", "SRP_ugL", 
                     "DOC_mgL", "DIC_mgL", "DC_mgL")

# Initialize an empty list to store results
max_depths <- list()


DCM_final <- final_data0 |>
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
    meta_width = mean(meta_width, .na.rm = TRUE),
    WaterLevel_m = mean(WaterLevel_m, na.rm = TRUE),
    maxdaily_airtemp = mean(maxdaily_airtemp, na.rm = TRUE),
    mindaily_airtemp = mean(mindaily_airtemp, na.rm = TRUE),
    precip_daily = mean(precip_daily, na.rm = TRUE),
    .groups = "drop"  # Ungroup to prevent grouping issues in the following steps
  )|>
  filter(Totals_DCM_conc > 20)

# Loop through the depth_variables to calculate the depth at which the maximum value occurs for each date
for (var in depth_variables) {
  DCM_final <- DCM_final |>
    left_join(
      final_data0 |>
        filter(month(Date) >= 4, month(Date) < 10) |>
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
      final_data0 |>
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
  mutate(DayOfYear = yday(Date)) |>
  select(Date, Totals_DCM_conc, Totals_DCM_depth, meta_width, peak.magnitude,
         secchi_PZ, PAR_PZ, PZ, Zeu, DCM_buoyancy_freq, thermocline_depth, WaterLevel_m, 
         everything())|>  # Include all max and min depth columns added
  rename_with(~ gsub("max_depth_(.*)", "max_\\1_depth", .), starts_with("max_depth_"))|>
  rename_with(~ gsub("min_depth_(.*)", "min_\\1_depth", .), starts_with("min_depth_"))

#write.csv(DCM_final,"./DCM_final.csv",row.names = FALSE)




####correlation function####

correlations <- function(year1, year2) {
  DCM_final_cor <- DCM_final |>
    filter(year(Date) >= {{year1}}, year(Date) <= {{year2}}) |>
    filter(month(Date) > 4, month(Date) < 10) |>
    filter(Totals_DCM_conc > 20)
  
  drivers_cor <- cor(DCM_final_cor[,c(2:69)],
                     method = "spearman", use = "pairwise.complete.obs")
 
  list(drivers_cor = drivers_cor, DCM_final_cor = DCM_final_cor)

}

#cutoff 0.7
results <- correlations(2014, 2023)
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
#these are the days that the max TotalConc_ugL occurs. The biggest bloom. 
blooms <- final_data0|>
  group_by(year(Date))|>
  mutate(bloommax = if_else(TotalConc_ugL == max(TotalConc_ugL), TRUE, NA_real_))|>
  ungroup()|>
  filter(bloommax == TRUE)|>
  group_by(Date)|>
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) |>
  ungroup()


DCM_final_maxdays_cor<- DCM_final|>
  filter(Date %in% c("2014-09-25", "2015-08-08", "2016-06-16", "2017-07-20", "2018-05-24", "2019-05-16", "2020-08-20", "2021-08-09","2022-09-19", "2023-06-19"))

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
blooms <- final_data0|>
  group_by(year(Date))|>
  mutate(bloommax = if_else(TotalConc_ugL == max(TotalConc_ugL), TRUE, NA_real_))|>
  ungroup()|>
  filter(bloommax == TRUE)|>
  group_by(Date)|>
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) |>
  ungroup()
  

#"2014-08-13" "2015-08-08" "2016-06-16" "2017-07-20" "2018-08-16" "2019-06-06" "2020-09-16" "2021-08-09" "2022-08-01" "2023-07-31"
#change date to see correlations for the singular day that max was the biggest
daily_cor <- final_data0|>
  filter(Date %in% c("2019-06-06"))|>
  select("Depth_m", "TotalConc_ugL", "TotalConc_ugL", "SFe_mgL", "TFe_mgL", "SMn_mgL", "SCa_mgL",
         "TCa_mgL", "TCu_mgL", "SBa_mgL", "TBa_mgL",
         "CO2_umolL", "CH4_umolL", "DO_mgL",
         "DOsat_percent", "Cond_uScm", "ORP_mV", "pH", "np_ratio", "TN_ugL", "TP_ugL", 
         "NH4_ugL", "NO3NO2_ugL", "SRP_ugL", "DOC_mgL", "DIC_mgL", 
         "DC_mgL", "PAR_LAP", "PAR_umolm2s", "sec_LAP", "Temp_C", "buoyancy_freq")

daily_cor_result <- cor(daily_cor[,c(1:31)], method = "spearman", use = "pairwise.complete.obs")
  
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




####flora data availability####

#days on the x axis, years on the y axis
plot_dat <- final_data_alldates %>%
  filter(!is.na(TotalConc_ugL)) %>%
  mutate(Year = year(Date), 
         DayOfYear = yday(Date))|> # Extract year and day of the year
  select(Date, Year, DayOfYear, TotalConc_ugL, Depth_m, peak.width, peak.magnitude)

# Find the maximum TotalConc_ugL value for each year
max_total_per_year <- plot_dat %>%
  group_by(year(Date)) %>%
  slice(which.max(TotalConc_ugL)) %>%
  ungroup()

# Plot: x-axis is DayOfYear, y-axis is Year, with a line and highlighted points
ggplot(plot_dat, aes(x = DayOfYear, y = as.factor(Year), group = Year)) +
  geom_line() +  # Line for each year
  geom_point() +  # Data points
  geom_point(data = max_total_per_year, aes(x = DayOfYear, y = as.factor(Year)), 
             color = "red", size = 3) +  # Highlight max points in red
  geom_text(data = max_total_per_year, 
            aes(x = DayOfYear, y = as.factor(Year), 
                label = paste0("Max: ", round(TotalConc_ugL, 2), " µg/L\nDepth: ", Depth_m, " m")), 
            vjust = 1.5, hjust = 0.5, color = "black", size = 3) +  # Smaller text and place below the point
  theme_bw() +
  labs(x = "Day of Year", y = "Year", title = "Fluoroprobe Data Availability") +
  scale_x_continuous(breaks = seq(1, 365, by = 30), limits = c(1, 365)) +  # Set x-axis limits and breaks
  theme(panel.grid.minor = element_blank())  # Optional: remove minor grid lines

####DCM depth every year####
# Find the maximum TotalConc_ugL value for each day
max_total_per_day <- plot_dat %>%
  group_by(Date) %>%
  slice(which.max(TotalConc_ugL)) %>%
  filter(DayOfYear > 133, DayOfYear < 285, TotalConc_ugL >20)|>
  ungroup()

ggplot(max_total_per_day, aes(x = DayOfYear, y = Depth_m, group = Year)) +
  geom_line() +
  geom_point(data = max_total_per_year, aes(x = DayOfYear, y = Depth_m), 
             color = "red", size = 3) +  # Highlight max points in red
  geom_text(data = max_total_per_year, 
            aes(x = DayOfYear, y = Depth_m, 
                label = paste0("Max: ", round(TotalConc_ugL, 2), " µg/L\nDepth: ", Depth_m, " m")), 
            vjust = -.5, hjust = 0.5, color = "black", size = 2.5) +  # Adjust text position
  theme_bw() +
  labs(x = "Day of Year", y = "Depth (m)", title = "DCM Depths Across Years (Only Showing Data with Totals > 20)") +
  scale_y_reverse(limits = c(10, 0)) +  # Invert y-axis from 0 to 10
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +  # Adjust x-axis breaks
  facet_wrap(~ Year, ncol = 2) +  # Create separate panels for each year
  theme(panel.grid.minor = element_blank())  # Optional: remove minor grid lines

####peak width every year####

ggplot(max_total_per_day, aes(x = DayOfYear, y = peak.width, group = Year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of Year", y = "Peak Width (m)", title = "DCM Widths Across Years (Only Showing Data with Totals > 20)") +
  scale_y_continuous(limits = c(0, 4)) +  
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +  # Adjust x-axis breaks
  facet_wrap(~ Year, ncol = 2) +  # Create separate panels for each year
  theme(panel.grid.minor = element_blank())  # Optional: remove minor grid lines

####peak magnitude####
ggplot(max_total_per_day, aes(x = DayOfYear, y = peak.magnitude, group = Year)) +
  geom_line() +
  theme_bw() +
  labs(x = "Day of Year", y = "Peak Magnitude (m)", title = "Peak Magnitude Across Years (Only Showing Data with Totals > 20)") +
  scale_y_continuous(limits = c(0, 150)) +  
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +  # Adjust x-axis breaks
  facet_wrap(~ Year, ncol = 2) +  # Create separate panels for each year
  theme(panel.grid.minor = element_blank())  # Optional: remove minor grid lines



####boxplots depth of DCM####

#for june, july, august
boxplot_Data <- DCM_final |>
  filter(Totals_DCM_conc > 20) |>
  filter(month(Date)>5, month(Date)<9) |>
  mutate(Year = year(Date), Month = month(Date))

# Calculate max_legend_value for the color scale limits
max_legend_value <- max(boxplot_Data$Totals_DCM_conc, na.rm = TRUE)

# Create the multi-panel boxplot with an overlay of colored points for Totals_DCM_conc
ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
                         y = Totals_DCM_depth)) +
  geom_boxplot() +
  geom_point(aes(color = Totals_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  facet_wrap(~ Year) +  # Create a panel for each year
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  scale_y_reverse(name = "DCM Depth (inverted)") +  # Reverse the y-axis
  ylim(10, 0) +  # Set the y-axis limits, reversing the range
  labs(x = "Month", y = "DCM Depth", color = "TotalConc ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

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
  ggtitle(label = "DCM Depths only displaying Totals > 20") +
  ylim(10, 0) +  # Set the y-axis limits, reversing the range
  labs(x = "Year", y = "DCM Depth", color = "TotalConc ugL") +  # Label the legend
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
  labs(x = "Month", y = "Peak Width", color = "Total ugL") +  # Label the legend
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
  ggtitle(label = "Peak Width only displaying TotalConc > 20") +
  ylim(0, 5) +  # Set the y-axis limits
  labs(x = "Year", y = "Peak Width", color = "TotalConc ugL") +  # Label the legend
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
  labs(x = "Month", y = "Peak Magnitude", color = "TotalConc ugL") +  # Label the legend
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
  ggtitle(label = "Peak Magnitudes only displaying TotalConc > 20")+
  labs(x = "Year", y = "Peak Magnitude", color = "TotalConc ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### plot TC, DCM, and PAR####
ggplot(DCM_final, aes(x = thermocline_depth, y = Totals_DCM_depth, color = PZ))+
  geom_point()

looking<- DCM_final%>%
  filter(thermocline_depth <3)|>
  select(Totals_DCM_depth, thermocline_depth, Date, WaterLevel_m)

ggplot(DCM_final, aes(x = Totals_DCM_depth, y = PZ, color = thermocline_depth)) +
  geom_point() +
#  scale_x_log10() +
#  scale_y_log10() +
  labs(title = "Thermocline vs PZ",
       x = "Totals_DCM_depth",
       y = "PZ ")


#"2016-05-19" "2016-08-04" dates with Totals_DCM_depth >9.5

####RandomForest Anually####

#"We constructed a RF of 1500 trees for each of the two response
#variables using 1% PAR depth (m), DOC concentration (mg L21),
#thermocline depth (m), metalimnion thickness (m),
#buoyancy frequency at the thermocline (s21),
#lake surface area (log10(km2)), and maximum depth (log10(m))
#as predictors included in each analysis."
#Leach Patterns and Drivers
library(randomForest)
library(missForest)

#trying within a year
yearDCM_final <- DCM_final |>
 # filter(year(Date) == 2023) |>
  mutate(DOY = yday(Date)) |>
  select(-peak.top, -max_pH_depth, -min_pH_depth, -min_Cond_uScm_depth, -max_Cond_uScm_depth, -peak.magnitude, -peak.bottom, -secchi_PZ, -PAR_PZ, -max_Cond_uScm_depth, -DayOfYear, -min_Cond_uScm_depth, -min_CH4_umolL_depth, -max_CH4_umolL_depth)|>
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
  model_rf <- randomForest(Totals_DCM_depth ~ ., data = train_data_imputed_z, ntree = 500, importance = TRUE)
  
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
  
  
  
  
  


####RF TotalConc magnitude####
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
  
  
  

####ARIMA####
  # adapted from Carly's code
  
  
  #pchip for interpolation - signal package - will need to specify to use this package 
  # pracma::pchip is another 
  # zoo::na.approx
  rm(list=ls(all=TRUE))
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(forecast)
  library(ggplot2)
  library(fable)
  library(tsibble)
  
  metals <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/455/8/9c8c61b003923f4f03ebfe55cea8bbfd")
  
  #metals <- read_csv("https://raw.githubusercontent.com/carlybauer/Reservoirs/refs/heads/master/Data/DataAlreadyUploadedToEDI/EDIProductionFiles/MakeEMLmetals/2024/Metals_2014_2023.csv")
  # changes to the proper format
  metals$DateTime <- as.POSIXct(metals$DateTime)
  
  metals <- metals %>% 
    mutate(Year = year(DateTime)) %>% # mutate creates new Year variable from DateTime variable to be able to look at specific years
    mutate(Date_real = date(DateTime)) %>% 
    filter(Reservoir == "FCR", Site == 50, Year >= 2020, Depth_m == 1.6 )
  
  m_data <- metals %>% 
    select(-Reservoir, -Site, -starts_with("Flag"), -DateTime)
  
  m_data<- m_data %>%
    select(Date_real, Year, Depth_m, SBa_mgL, TBa_mgL, 
           SCu_mgL, TCu_mgL, SSr_mgL, TSr_mgL, SAl_mgL, TAl_mgL) %>% 
    drop_na()
  
  # make sure it's recognizing data as a time series
  df <- as_tsibble(m_data)
  
  # put into week intervals 
  df <- df %>%
    mutate(Week = week(Date_real))
  
  start_date <- as.Date("2020-01-01")
  end_date <- as.Date("2023-12-31")
  
  weekly_dates <- data.frame(
    Date_fake = seq.Date(from = start_date, to = end_date, by = "week")
  ) %>%
    mutate(Year = year(Date_fake),
           Week = week(Date_fake))
  
  df_weekly <- weekly_dates %>%
    left_join(df, by = c("Year", "Week"))
  
  # need to gap fill - possible solution below 
  # linear interpolation, carrying forward the last observation
  # df_weekly <- df_weekly %>%
  #   arrange(Date) %>%
  #   tidyr::fill(SBa_mgL:TAl_mgL, .direction = "downup")
  
  
  # attr(df, 'interval') <- tsibble::new_interval(.regular = FALSE)
  # df <- df %>% 
  #   tsibble::update_tsibble(.data, regular = TRUE)
  
  # Check for duplicates in Date_fake
  duplicates(df_weekly, index = Date_fake)
  # Summarize by mean for each Date_fake if there are duplicates
  df_weekly <- df_weekly %>%
    group_by(Date_fake) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
    ungroup()
  
  # Convert to a tsibble, specifying Date_fake as the index
  df_weekly <- df_weekly %>%
    as_tsibble(index = Date_fake)
  
  # Fill gaps in the data
  df_weekly <- df_weekly %>%
    fill_gaps()
  
  # df_weekly <- as_tsibble(df_weekly)
  # df_weekly <- df_weekly %>% 
  #    tsibble::fill_gaps()
  
  my.arima <- df_interp %>%
    model(Ba_ARIMA = fable::ARIMA(TBa_mgL)) 
  my.arima #TBa_mgL (3,1,0)
  
  # Interpolate only the TBa_mgL column and add it back to df_weekly
  df_interp <- df_weekly %>%
    mutate(across(SBa_mgL:TAl_mgL, ~ forecast::na.interp(.)))
  
  # Extract forecast data with mean and confidence intervals
  forecast_data <- forecast(my.arima, h = 20) %>%
    as_tibble() %>%  # Convert forecast to a regular data frame
    transmute(
      Date_fake,
      mean = as.numeric(.mean),
      lower_80 = as.numeric(hilo(TBa_mgL, level = 80)[[1]]$lower),
      upper_80 = as.numeric(hilo(TBa_mgL, level = 80)[[1]]$upper),
      lower_95 = as.numeric(hilo(TBa_mgL, level = 95)[[1]]$lower),
      upper_95 = as.numeric(hilo(TBa_mgL, level = 95)[[1]]$upper)
    )
  
  # Prepare the actual (known) data
  known_data <- df_weekly %>%
    filter(!is.na(TBa_mgL)) %>%
    select(Date_fake, TBa_mgL) %>%
    rename(mean = TBa_mgL)  # Rename for consistency in plotting
  
  # Combine known and forecasted data
  combined_data <- bind_rows(
    known_data %>% mutate(type = "Known"),
    forecast_data %>% mutate(type = "Forecast")
  )
  
  # Plot known data and forecast with confidence intervals
  ggplot(combined_data, aes(x = Date_fake, y = mean, color = type)) +
    geom_line() +
    geom_ribbon(data = filter(combined_data, type == "Forecast"),
                aes(ymin = lower_80, ymax = upper_80), fill = "blue", alpha = 0.2) +
    geom_ribbon(data = filter(combined_data, type == "Forecast"),
                aes(ymin = lower_95, ymax = upper_95), fill = "blue", alpha = 0.1) +
    scale_color_manual(values = c("Known" = "black", "Forecast" = "blue")) +
    labs(title = "TBa_mgL Forecast with Historical Data", y = "TBa_mgL", x = "Date") +
    theme_minimal()
  
  # Generate the forecast and extract the mean and confidence intervals
  forecast_data <- forecast(my.arima, h = 20) %>%
    mutate(
      mean = as.numeric(.mean),
      lower_80 = hilo(TBa_mgL, level = 80)$lower,
      upper_80 = hilo(TBa_mgL, level = 80)$upper,
      lower_95 = hilo(TBa_mgL, level = 95)$lower,
      upper_95 = hilo(TBa_mgL, level = 95)$upper
    )
  
  # Plot the forecast
  ggplot(forecast_data, aes(x = Date_fake)) +
    geom_line(aes(y = mean), color = "blue") +
    geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "blue", alpha = 0.2) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "blue", alpha = 0.1) +
    labs(title = "Forecast of TBa_mgL", y = "TBa_mgL", x = "Date") +
    theme_minimal()
  
  # using forecast package
   my.arima <-auto.arima(df$TBa_mgL)
   plot(forecast(my.arima, h=20))
  
  #finding time differences to maybe interpolate
  # diff(df$DateTime)
  
  # then check residuals using acf 
  # and portmanteau test 
  # Acf(my.arima$residuals)
  
  # mydata <- m_data
  # colnames(mydata)
  # 
  # driver_cor <- cor(mydata, method = "spearman", use = "pairwise.complete.obs")
  # driver_cor[lower.tri(driver_cor)]=""
  # driver_cor
  # write.csv(driver_cor, file = "Metal_driver_cor_ref.csv", row.names = FALSE)
  
  
  
  # ACF -> checking for stationarity
  ## drop quickly to 0 = stationary 
  # ## decreases slowly = non stationary 
  # Acf(my.fp.data$SBa_mgL)
  
   
####LSTM####
  library(keras)
  
  # Define number of timesteps and features
  timesteps <- # number of past observations per sequence
    n_features <- length(cols_to_interpolate)  # number of features
  
  # Build the model
  model <- keras_model_sequential()
  model %>%
    layer_lstm(units = 50, return_sequences = TRUE, input_shape = c(timesteps, n_features)) %>%
    layer_lstm(units = 25) %>%
    layer_dense(units = 1)  # adjust units based on your prediction target
  
  # Compile the model
  model %>% compile(loss = "mse", optimizer = "adam")  # adjust loss and optimizer
  
  # Train the model
  model %>% fit(train_data_sequences, train_targets, epochs = 100, batch_size = 32)
  
#layer_lstm: Defines LSTM layers with the number of units (neurons) and options like return_sequences for feeding information through layers.
#layer_dense: Final layer for predicting the target variable (e.g., Totals_DCM_conc).
#compile: Sets the loss function (e.g., mean squared error for regression) and optimizer (e.g., Adam).
#fit: Trains the model with your prepared training sequences and target values (train_data_sequences, train_targets).
  
  
  #prediction
  predicted_values <- predict(model, test_data_sequences)
  
  
  




