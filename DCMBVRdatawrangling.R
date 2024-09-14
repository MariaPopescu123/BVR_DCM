# Maria DCM BVR data set
#data includes chlorophyll maxima, PAR, secchi, light attenuation, metals, ghgs, nutrients

beep <- function(){
  system("rundll32 user32.dll,MessageBeep") #just putting this here so I can be notified when something that takes a while to run is done
}  

pacman::p_load(tidyverse, lubridate, akima, reshape2, 
               gridExtra, grid, colorRamps, RColorBrewer, rLakeAnalyzer,
               reader, cowplot, dplyr, tidyr, ggplot2, zoo, purrr, beepr, forecast)

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

#PAR
PAR_profiles <- read.csv("PAR_profiles.csv")

#chemistry
chemistry <- read.csv("./chemistry.csv")

#meteorological data from FCR https://portal.edirepository.org/nis/mapbrowse?packageid=edi.389.8
options(timeout = 300)
metdata <- read.csv("./metdata.csv")

#bathymetry
bath <- read.csv("./bath.csv")

#waterlevel data
wtrlvl <- read.csv("./wtrlvl.csv") 

#adding columns with total_conc max and the depth at which it occurs
DCM_BVRdata <- current_df %>% 
  filter(Reservoir == "BVR")%>%
  mutate(Date  = as_date(DateTime)) |> 
  group_by(CastID) %>%
  summarise(DCM_totalconc = max(TotalConc_ugL, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(current_df, by = "CastID") %>% # Join back the original data to keep all columns
  mutate(DCM_depth = ifelse(TotalConc_ugL == DCM_totalconc, Depth_m, NA_real_))%>% # Add DCM_depth column
  mutate(cyanosatDCM = ifelse(TotalConc_ugL == DCM_totalconc, Bluegreens_ugL, NA_real_))|>  #the concentration of cyanos at DCM
  group_by(CastID)|>
  mutate(Bluegreens_DCM_conc = max(Bluegreens_ugL, na.rm = TRUE))|> #concentration of bluegreens at bluegreens DCM
  mutate(Bluegreens_DCM_depth = ifelse(Bluegreens_ugL == Bluegreens_DCM_conc, Depth_m, NA_real_))
  

DCM_BVRdata <- DCM_BVRdata %>%
  group_by(CastID) %>%
  mutate(DCM_depth = ifelse(!is.na(DCM_depth), DCM_depth, NA_real_)) %>%  # Ensure DCM_depth is NA where condition didn't match
  fill(DCM_depth, .direction = "downup") %>%  # Fill NA values in DCM_depth column within each CastID
  mutate(cyanosatDCM = ifelse(!is.na(cyanosatDCM), cyanosatDCM, NA_real_)) %>%  # Ensure cyanosatDCM is NA where condition didn't match
  fill(cyanosatDCM, .direction = "downup") %>%  # Fill NA values in cyanosatDCM column within each CastID
  fill(Bluegreens_DCM_conc, .direction = "downup")|>
  fill(Bluegreens_DCM_depth, .direction = "downup")|>
  ungroup()%>%
  relocate(DCM_depth, .before = 3)%>%
  relocate(cyanosatDCM, .before = 3)%>%
  mutate(DOY = yday(DateTime))%>%
  mutate(Date  = as_date(DateTime)) |> 
  select(Date, DateTime, CastID, DCM_totalconc, DCM_depth, Depth_m, cyanosatDCM, Bluegreens_DCM_conc, Bluegreens_DCM_depth, Bluegreens_ugL, TotalConc_ugL,  Temp_C)

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
  mutate(DCM = (DCM_depth==Depth_m))|>
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

####Adding PAR, DO, DOsat_percent, cond, ORP, pH ####

PAR_profiles_filtered <- PAR_profiles |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(DO_mgL = mean(DO_mgL, na.rm = TRUE),
            PAR_umolm2s = mean(PAR_umolm2s, na.rm = TRUE), 
            DOsat_percent = mean(DOsat_percent, na.rm = TRUE),
            Cond_uScm = mean(Cond_uScm, na.rm = TRUE),
            ORP_mV = mean(ORP_mV, na.rm = TRUE),
            pH = mean(pH, na.rm = TRUE))

DCM_BVRdata <- DCM_BVRdata %>%
    mutate(Date = as.Date(Date, format = "%Y-%m-%d")) # Adjust format as needed
  
  
interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(PAR_profiles_filtered) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(interp_PAR_umolm2s = zoo::na.approx(PAR_umolm2s, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_DO_mgL = zoo::na.approx(DO_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_DOsat_percent = zoo::na.approx(DOsat_percent, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_Cond_uScm = zoo::na.approx(Cond_uScm, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_ORP_mV = zoo::na.approx(ORP_mV, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_pH = zoo::na.approx(pH, x = Depth_m, na.rm = FALSE)) |>
  ungroup()|> # Ensure the grouping is removed for the final merge'
  filter(year(Date) < 2018)

final_data <- DCM_BVRwmetalsghgssecchilight %>%
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  select(-DO_mgL, -PAR_umolm2s, -DOsat_percent, -Cond_uScm, -ORP_mV, -pH) %>% # Remove unnecessary columns
  filter(Depth_m %in% DCM_BVRdata$Depth_m) # Keep only rows with depths present in flora data

CTDfiltered <- CTD|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(DO_mgL = mean(DO_mgL, na.rm = TRUE),
            PAR_umolm2s = mean(PAR_umolm2s, na.rm = TRUE), 
            DOsat_percent = mean(DOsat_percent, na.rm = TRUE),
            Cond_uScm = mean(Cond_uScm, na.rm = TRUE),
            ORP_mV = mean(ORP_mV, na.rm = TRUE),
            pH = mean(pH, na.rm = TRUE))|>
  filter(year(Date) > 2017)

interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(CTDfiltered) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(interp_PAR_umolm2s = zoo::na.approx(PAR_umolm2s, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_DO_mgL = zoo::na.approx(DO_mgL, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_DOsat_percent = zoo::na.approx(DOsat_percent, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_Cond_uScm = zoo::na.approx(Cond_uScm, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_ORP_mV = zoo::na.approx(ORP_mV, x = Depth_m, na.rm = FALSE)) |>
  mutate(interp_pH = zoo::na.approx(pH, x = Depth_m, na.rm = FALSE)) |>
  ungroup()

final_data2 <- final_data %>% #remember to remove the 2 after this works so i can call it in the next thing
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  select(-DO_mgL, -PAR_umolm2s, -DOsat_percent, -Cond_uScm, -ORP_mV, -pH) %>% # Remove unnecessary columns
  filter(Depth_m %in% DCM_BVRdata$Depth_m) # Keep only rows with depths present in flora data

final_data <- final_data2|>
  mutate(interp_PAR_umolm2s = coalesce(interp_PAR_umolm2s.x, interp_PAR_umolm2s.y)) %>%
  mutate(interp_DO_mgL = coalesce(interp_DO_mgL.x, interp_DO_mgL.y)) %>%
  mutate(interp_DOsat_percent = coalesce(interp_DOsat_percent.x, interp_DOsat_percent.y)) %>%
  mutate(interp_Cond_uScm = coalesce(interp_Cond_uScm.x, interp_Cond_uScm.y)) %>%
  mutate(interp_ORP_mV = coalesce(interp_ORP_mV.x, interp_ORP_mV.y)) %>%
  mutate(interp_pH = coalesce(interp_pH.x, interp_pH.y))|>
  select(-interp_PAR_umolm2s.x, -interp_PAR_umolm2s.y,
         -interp_DO_mgL.x, -interp_DO_mgL.y,
         -interp_DOsat_percent.x, -interp_DOsat_percent.y,
         -interp_Cond_uScm.x, -interp_Cond_uScm.y, 
         -interp_ORP_mV.x, -interp_ORP_mV.y, 
         -interp_pH.x, -interp_pH.y)|>
  relocate(Reservoir, .before = 1)|>
  relocate(Site, .before = 1)

####calculating PAR_LAP (Light availability percentage using interpolated PAR)####
final_data$log_PAR <- log(final_data$interp_PAR_umolm2s)

final_data <- final_data|>
  relocate(log_PAR, .after = interp_PAR_umolm2s)

#first calculating new Kd value from provided PAR (that I interpolated from CTD data and PAR_profiles)
final_data <- final_data |> 
  group_by(CastID) |> 
  filter(!all(is.na(log_PAR)) & !all(is.na(Depth_m))) |>  # Remove groups with all NAs
  filter(!is.na(log_PAR) & !is.na(Depth_m) & is.finite(log_PAR) & is.finite(Depth_m)) |>  # Remove rows with NA, NaN, or Inf in relevant columns
  mutate(PAR_K_d = abs(map_dbl(list(lm(log_PAR ~ Depth_m, data = cur_data())), ~coef(.x)["Depth_m"])))

#now calculating light availability percentage from PAR_K_d
final_data <- final_data |>
  group_by(CastID)|>
  mutate(PAR_LAP = 100* exp(-PAR_K_d * Depth_m))

####Secchi PZ####
#using secchi_PZ because the data for PAR between CTD and YSI is too different (for example look at PAR_profiles filtered and CTD filtered for 2018-5-4)
#if I want to calculate PZ for specific years it would be ok but across all years no
final_data <- final_data |>
  mutate(secchi_PZ = 2.8*Secchi_m)|>
  relocate(secchi_PZ, .before = sec_LAP)|>
  relocate(PAR_LAP, .after = sec_LAP)

####Temps from 2017-2019  ####
CTDfiltered <- CTD |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(Temp_C = mean(Temp_C, na.rm = TRUE))|>
  filter(year(Date)>2016, year(Date)<2020)

interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(CTDfiltered) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(Temp_C = zoo::na.approx(Temp_C, x = Depth_m, na.rm = FALSE)) |>
  ungroup()

final_data2 <- final_data %>% 
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  filter(Depth_m %in% DCM_BVRdata$Depth_m) # Keep only rows with depths present in flora data

final_data <- final_data2|>
  mutate(Temp_C = coalesce(Temp_C.x, Temp_C.y))|>
  select(-Temp_C.x, -Temp_C.y)

looking <- final_data|>
  filter(year(DateTime) == 2019)|>
  select(Date, Depth_m, Temp_C)

#### pH ####
# Adding pH for 2017
CTDfiltered <- CTD |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(pH = mean(pH, na.rm = TRUE))|>
  filter(year(Date) == 2017|year(Date) == 2020)|>
  filter(!is.na(pH))

interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(CTDfiltered) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(pH = zoo::na.approx(pH, x = Depth_m, na.rm = FALSE)) |>
  ungroup()

final_data2 <- final_data %>% 
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  filter(Depth_m %in% DCM_BVRdata$Depth_m) # Keep only rows with depths present in flora data

final_data <- final_data2|>
  mutate(interp_pH = coalesce(interp_pH, pH))|>
  select(-pH)

# pH for 2021 from the PAR profiles  #
PAR_profiles_filtered <- PAR_profiles |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(pH = mean(pH, na.rm = TRUE))|>
  filter(year(Date) == 2021)|>
  filter(!is.na(pH))

interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(PAR_profiles_filtered) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(pH = zoo::na.approx(pH, x = Depth_m, na.rm = FALSE)) |>
  ungroup()

final_data2 <- final_data %>% 
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  filter(Depth_m %in% DCM_BVRdata$Depth_m) # Keep only rows with depths present in flora data

final_data <- final_data2|>
  mutate(interp_pH = coalesce(interp_pH, pH))|>
  select(-pH)

#### ORP  ####
#there just is no data for 2021
CTDfiltered <- CTD |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(ORP_mV = mean(ORP_mV, na.rm = TRUE))|>
  filter(year(Date) == 2017|year(Date) == 2021)|>
  filter(!is.na(ORP_mV))

interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(CTDfiltered) |> 
  arrange(Date, Depth_m) |> # Sort data by Date and Depth_m
  group_by(Date)  |> # Group data by Date for interpolation
  mutate(ORP_mV = zoo::na.approx(ORP_mV, x = Depth_m, na.rm = FALSE)) |>
  ungroup()

final_data2 <- final_data %>% 
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  filter(Depth_m %in% DCM_BVRdata$Depth_m) # Keep only rows with depths present in flora data

final_data <- final_data2|>
  mutate(interp_ORP_mV = coalesce(interp_ORP_mV, ORP_mV))|>
  select(-ORP_mV)

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
final_data0 <- final_data %>%
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  select(-TN_ugL,-TP_ugL, -NH4_ugL, -NO3NO2_ugL,-SRP_ugL, -DOC_mgL, -DIC_mgL, -DC_mgL) %>% # Remove unnecessary columns
  filter(Depth_m %in% DCM_BVRdata$Depth_m)|> # Keep only rows with depths present in DCMdata
  mutate(Reservoir = "BVR")|>
  mutate(Site = 50)|>
  relocate(Reservoir, .before = 1)|>
  relocate(Site, .before = 1)

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
final_data0 <- final_data0 %>%
  mutate(np_ratio = calculate_np_ratio(interp_TN_ugL,interp_TP_ugL))|>
  relocate(np_ratio, .before = interp_TN_ugL)
}



#### Visualizing metdata  ####

metdata0 <- metdata|>
  mutate(Date = as_date(DateTime))|>
  mutate(DOY = yday(Date))|>
  relocate(DOY, .before = DateTime)|>
  relocate(Date, .before = Date)|>
  mutate(DateTime = ymd_hms(DateTime))

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
metdata0 <- metdata0 |> 
  group_by(DOY, year(DateTime)) |> 
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
metdata0 <- metdata0 |> 
  group_by(DOY, year(DateTime)) |> 
  mutate(daily_airtempavg = mean(AirTemp_C_Average, na.rm = TRUE))|>
  mutate(maxdaily_airtemp = max(AirTemp_C_Average, na.rm = TRUE))|>
  mutate(mindaily_airtemp = min(AirTemp_C_Average, na.rm = TRUE))|>
  ungroup()|>
  relocate(daily_airtempavg, .before = AirTemp_C_Average)

b1 <- metplots(2015, daily_airtempavg, maxx = 40)
b2 <- metplots(2016, daily_airtempavg, maxx = 40)
b3 <- metplots(2017, daily_airtempavg, maxx = 40)
b4 <- metplots(2018, daily_airtempavg, maxx = 40)
b5 <- metplots(2019, daily_airtempavg, maxx = 40)
b6 <- metplots(2020, daily_airtempavg, maxx = 40)
b7 <- metplots(2021, daily_airtempavg, maxx = 40)
b8 <- metplots(2022, daily_airtempavg, maxx = 40)
b9 <- metplots(2023, daily_airtempavg, maxx = 40)


dailyaveragetemps<- plot_grid(
  b1, b2, b3,
  b4, b5, b6, 
  b7, b8, b9,
  ncol = 3
)

print(dailyaveragetemps)

#### precip and temp to final_data0 #### 

metdata_join <- metdata0|>
  select(Date, precip_daily, daily_airtempavg, maxdaily_airtemp, mindaily_airtemp)|>
  distinct()

final_data0 <- final_data0|>
  left_join(metdata_join, by = c("Date"), relationship = "many-to-many")|>
  filter(Date %in% final_data0$Date)


####thermocline ####

# Dataframe with thermocline
thermocline_df <- final_data0 |>
  group_by(CastID, Depth_m) |>
  summarize(Temp_C = mean(Temp_C, na.rm = TRUE)) |>
  ungroup() |>
  group_by(CastID) |>
  mutate(thermocline_depth = thermo.depth(
    Temp_C, 
    Depth_m, 
    Smin = 0.1, 
    seasonal = TRUE, 
    index = FALSE,
    mixed.cutoff = 1
  )) 

#add to final_data0
final_datatest <- final_data0|>
  left_join(thermocline_df, by = c("CastID", "Depth_m", "Temp_C"), relationship = "many-to-many")

final_data0 <- final_datatest|>
  group_by(CastID)|>
  fill(thermocline_depth, .direction = "updown")|>
  ungroup()|>
  relocate(thermocline_depth, .before = Temp_C)

####Buoyancy Frequency ####

BVRbath <- bath|>
  filter(Reservoir == "BVR")

final_data0 <- final_data0|>
  group_by(CastID)|>
  mutate(buoyancy_freq = c(buoyancy.freq(Temp_C, Depth_m), NA))|>#added for padding for the last value
  relocate(buoyancy_freq, .before = thermocline_depth)
#need to make sure this makes sense

####Waterlevels####

wtrlvl2 <- wtrlvl|>
  mutate(Date = as.Date(DateTime))|>
  select(Date, WaterLevel_m)

final_data_water <- final_data0|>
  left_join(wtrlvl2, by = c("Date"), relationship = "many-to-many")|>
  mutate(
    # For rows where 'value' from wtrlvl2 is NA after the join,
    # find the closest date in wtrlvl2 and get the corresponding value
    WaterLevel_m = ifelse(
      is.na(WaterLevel_m),
      sapply(Date, function(d) {
        closest_date <- wtrlvl2$Date[which.min(abs(difftime(wtrlvl2$Date, d, units = "days")))]
        wtrlvl2$WaterLevel_m[wtrlvl2$Date == closest_date]
      }),
      WaterLevel_m
    )
  )

####Peak.width####
#use blue_mean not blue_median
#focusing on bluegreens

final_data0test <- final_data_water %>%
  group_by(CastID) %>%
  mutate(
    blue_med = median(Bluegreens_ugL, na.rm = TRUE),  # Calculate the mean, excluding NA values if necessary
    blue_mean = mean(Bluegreens_ugL, na.rm = TRUE),
    blue_sd = sd(Bluegreens_ugL, na.rm = TRUE), #calculate the standard deviation
    peak.top = as.integer(Depth_m <= Bluegreens_DCM_depth & Bluegreens_ugL > blue_mean),  # Create a binary indicator
    peak.bottom = as.integer(Depth_m >= Bluegreens_DCM_depth & Bluegreens_ugL > blue_mean),
    # Apply the condition: If Bluegreens_DCM_conc < 20, set peak.top and peak.bottom to 0
    peak.top = if_else(Bluegreens_DCM_conc < 20, 0, peak.top),
    peak.bottom = if_else(Bluegreens_DCM_conc < 20, 0, peak.bottom),
    peak.top = if_else(peak.top == 1, Depth_m, 0),    # Replace peak.top with Depth_m where peak.top == 1
    peak.bottom = if_else(peak.bottom == 1, Depth_m, 0),
    #Get the minimum peak.top value excluding 0, and replace Inf with NA if all values are NA or 0
    peak.top = min(peak.top[peak.top != 0], na.rm = TRUE),
    # Get the maximum peak.bottom value excluding 0, and replace -Inf with NA if all values are NA or 0
    peak.bottom = max(peak.bottom[peak.bottom != 0], na.rm = TRUE),
    peak.width = peak.bottom - peak.top,
    peak.width = if_else(is.infinite(peak.width), NA_real_, peak.width)
  ) %>%
  ungroup()  # Ungroup after mutations


####Peak.magnitude####

final_data_peaks <- final_data0test|>
  group_by(CastID)|>
  mutate(peak.magnitude = max(Bluegreens_ugL)-mean(Bluegreens_ugL))|>
  ungroup()|>
  select(-cyanosatDCM) #this is unnecessary. saying how many bluegreens there are at the DCM for total_conc

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



#make sure to change final dataframe to this
final_data0 <- final_data_peaks


write.csv(final_data0,"./final_data0.csv",row.names = FALSE)





####correlations####
#removed buoyancy_freq for now bc had -inf will come back to

####daily DCM dataframe with daily averages####
DCM_final <- final_data0 |>
  filter(month(DateTime) >= 4, month(DateTime) < 10) |>
  group_by(Date) |>
  mutate(across(c(interp_SFe_mgL, interp_TFe_mgL, interp_SMn_mgL, interp_SCa_mgL,
                  interp_TCa_mgL, interp_TCu_mgL, interp_SBa_mgL, interp_TBa_mgL,
                  interp_CO2_umolL, interp_CH4_umolL, interp_DO_mgL,
                  interp_DOsat_percent, interp_Cond_uScm, interp_ORP_mV, interp_pH, interp_TN_ugL, interp_TP_ugL, 
                  interp_NH4_ugL, interp_NO3NO2_ugL, interp_SRP_ugL, interp_DOC_mgL, interp_DIC_mgL, 
                  interp_DC_mgL), 
                ~ if_else(DCM == TRUE, .x, NA_real_), .names = "DCM_{.col}")) |>
  fill(starts_with("DCM_interp_"), .direction = "updown") |>
  mutate(DCM_np_ratio = if_else(DCM == TRUE, np_ratio, NA_real_)) |>
  fill(DCM_np_ratio, .direction = "updown") |>
  mutate(DCM_Temp_C = if_else(DCM == TRUE, Temp_C, NA_real_)) |>
  fill(DCM_Temp_C, .direction = "updown") |>
  mutate(DCM_buoyancy_freq = if_else(DCM == TRUE, buoyancy_freq, NA_real_)) |>
  fill(DCM_buoyancy_freq, .direction = "updown") |>
  select(Date, Bluegreens_DCM_conc, Bluegreens_DCM_depth, peak.top, peak.bottom, peak.width, peak.magnitude,
         sec_K_d,PAR_K_d, DCM_buoyancy_freq, thermocline_depth, DCM_Temp_C, DCM_np_ratio, WaterLevel_m,DCM_interp_SFe_mgL,
         DCM_interp_TFe_mgL, DCM_interp_SMn_mgL, DCM_interp_SCa_mgL,
         DCM_interp_TCa_mgL, DCM_interp_TCu_mgL, DCM_interp_SBa_mgL, DCM_interp_TBa_mgL,
         DCM_interp_CO2_umolL, DCM_interp_CH4_umolL,secchi_PZ, DCM_interp_DO_mgL,
         DCM_interp_DOsat_percent, DCM_interp_Cond_uScm, DCM_interp_ORP_mV, DCM_interp_pH, DCM_interp_TN_ugL, DCM_interp_TP_ugL, 
         DCM_interp_NH4_ugL, DCM_interp_NO3NO2_ugL, DCM_interp_SRP_ugL, DCM_interp_DOC_mgL, DCM_interp_DIC_mgL, 
         DCM_interp_DC_mgL)|>
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) |>
  ungroup()|>
  select(Date, Bluegreens_DCM_conc, Bluegreens_DCM_depth, peak.top, peak.bottom, peak.width, peak.magnitude,
         sec_K_d, PAR_K_d, DCM_buoyancy_freq, thermocline_depth, DCM_Temp_C, DCM_np_ratio, WaterLevel_m,DCM_interp_SFe_mgL,
         DCM_interp_TFe_mgL, DCM_interp_SMn_mgL, DCM_interp_SCa_mgL,
         DCM_interp_TCa_mgL, DCM_interp_TCu_mgL, DCM_interp_SBa_mgL, DCM_interp_TBa_mgL,
         DCM_interp_CO2_umolL, DCM_interp_CH4_umolL,secchi_PZ, DCM_interp_DO_mgL,
         DCM_interp_DOsat_percent, DCM_interp_Cond_uScm, DCM_interp_ORP_mV, DCM_interp_pH, DCM_interp_TN_ugL, DCM_interp_TP_ugL, 
         DCM_interp_NH4_ugL, DCM_interp_NO3NO2_ugL, DCM_interp_SRP_ugL, DCM_interp_DOC_mgL, DCM_interp_DIC_mgL, 
         DCM_interp_DC_mgL)


####correlation function####

correlations <- function(year1, year2) {
  DCM_final_cor <- DCM_final |>
    filter(year(Date) >= {{year1}}, year(Date) <= {{year2}}) |>
    filter(month(Date) > 4, month(Date) < 10) |>
    filter(Bluegreens_DCM_conc > 20)
  
  drivers_cor <- cor(DCM_final_cor[,c(2:38)],
                     method = "spearman", use = "pairwise.complete.obs")
 
  list(drivers_cor = drivers_cor, DCM_final_cor = DCM_final_cor)

}

#cutoff 0.7
results <- correlations(2022, 2022)
final_data_cor_results <- results$drivers_cor
final_data_cor_results[lower.tri(final_data_cor_results)] = ""
final_data_cor <- results$DCM_final_cor


####daily correlation, for choosing specific day####

#these are the days that the max Bluegreens_ugL occurs. The biggest bloom. 
blooms <- final_data0|>
  group_by(year(DateTime))|>
  mutate(bloommax = if_else(Bluegreens_ugL == max(Bluegreens_ugL), TRUE, NA_real_))|>
  ungroup()|>
  filter(bloommax == TRUE)|>
  group_by(Date)|>
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) |>
  ungroup()
  
#"2014-07-02" "2015-06-18" "2016-06-30" "2017-07-20" "2018-08-16" "2019-06-27" "2020-09-16" "2021-07-26" "2022-08-01" "2023-07-24"

daily_cor <- final_data0|>
  filter(Date %in% c("2022-08-01"))|>
  select(Depth_m, Bluegreens_ugL, TotalConc_ugL, interp_SFe_mgL, interp_TFe_mgL, interp_SMn_mgL, interp_SCa_mgL,
         interp_TCa_mgL, interp_TCu_mgL, interp_SBa_mgL, interp_TBa_mgL,
         interp_CO2_umolL, interp_CH4_umolL, interp_DO_mgL,
         interp_DOsat_percent, interp_Cond_uScm, interp_ORP_mV, interp_pH, np_ratio ,interp_TN_ugL, interp_TP_ugL, 
         interp_NH4_ugL, interp_NO3NO2_ugL, interp_SRP_ugL, interp_DOC_mgL, interp_DIC_mgL, 
         interp_DC_mgL, PAR_LAP, interp_PAR_umolm2s, sec_LAP, Temp_C, buoyancy_freq)

daily_cor_result <- cor(daily_cor[,c(1:32)], method = "spearman", use = "pairwise.complete.obs")
  
daily_cor_result[lower.tri(daily_cor_result)] = ""

####DCM depth boxplot####

library(ggplot2)
library(lubridate)
library(ggthemes)

# Filter and prepare the data for the years 2015-2023
boxplot_Data <- DCM_final |>
  filter(year(Date) >= 2014, year(Date) <= 2023, Bluegreens_DCM_conc > 30) |>
  filter(month(Date) > 5, month(Date) < 9) |>
  mutate(Year = year(Date), Month = month(Date))

# Calculate max_legend_value for the color scale limits
max_legend_value <- max(boxplot_Data$Bluegreens_DCM_conc, na.rm = TRUE)

# Create the multi-panel boxplot with an overlay of colored points for Bluegreens_DCM_conc
ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
                         y = Bluegreens_DCM_depth)) +
  geom_boxplot() +
  geom_point(aes(color = Bluegreens_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  facet_wrap(~ Year) +  # Create a panel for each year
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  scale_y_reverse(name = "DCM Depth (inverted)") +  # Reverse the y-axis
  ylim(10, 0) +  # Set the y-axis limits, reversing the range
  labs(x = "Month", y = "DCM Depth", color = "Bluegreens ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####DCM peak width box plots####

# Create the multi-panel boxplot with an overlay of colored points for Bluegreens_DCM_conc
ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
                         y = peak.width)) +
  geom_boxplot() +
  geom_point(aes(color = Bluegreens_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  facet_wrap(~ Year) +  # Create a panel for each year
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  labs(x = "Month", y = "Peak Width", color = "Bluegreens ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####DCM peak magnitude box plots####
ggplot(boxplot_Data, aes(x = factor(Month, labels = c("June", "July", "August")), 
                         y = peak.magnitude)) +
  geom_boxplot() +
  geom_point(aes(color = Bluegreens_DCM_conc), position = position_jitter(width = 0.2), size = 2) +  # Add points with color representing concentration
  facet_wrap(~ Year) +  # Create a panel for each year
  scale_color_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +  # Apply color gradient to points
  labs(x = "Month", y = "Peak Magnitude", color = "Bluegreens ugL") +  # Label the legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))









#converting dataframe to timeseries object
#but it says it has to be evenly spaced in time

# DO NOT CHANGE THIS ALSO USED FOR LINE IN HEATMAPS LATER
chlorophyll_data <- final_data0 |>
  filter(DCM_totalconc > 20)|> #choosing greater than 20 as the bloom
  filter(DCM == TRUE)|>
  filter(!(month(DateTime) == 8 & year(DateTime) == 2017 & Bluegreens_ugL < 35))|> #filter out weird drop in 2017
  filter(!(CastID == 395)) #filter out weird drop in 2016

#start with DCM depth 
#this WORKS but the time series is not stationary
#check to make sure it's weekly, do May-October

DCM_deptharima <- chlorophyll_data|>
  filter(year(DateTime) == 2015)|>
  select(DCM_depth)

data <- ts(DCM_deptharima)
arima_model <- auto.arima(data)

summary(arima_model)



#attempt at time series analysis

library(forecast)
library(MuMIn)

#square root transformation?

final_data0$bluegreens_sqrt <- sqrt(final_data0$Bluegreens_ugL)

# Check for autocorrelation
acf(final_data0$bluegreens_sqrt, lag.max = 10)

# Fit an AR(1) model with different environmental predictors
model <- auto.arima(final_data0$bluegreens_sqrt, xreg = final_data0[, c("np_ratio", "Temp_C")])

# Model selection using AICc
dredge(model, rank = "AICc")

# Best model summary
summary(best_model)



