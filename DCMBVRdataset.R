#Maria DCM BVR data set
#data includes chlorophyll maxima, PAR, secchi, light attenuation, metals, ghgs, nutrients

beep <- function(){
  system("rundll32 user32.dll,MessageBeep") #just putting this here so I can be notified when something that takes a while to run is done
}  

pacman::p_load(tidyverse, lubridate, akima, reshape2, 
               gridExtra, grid, colorRamps, RColorBrewer, rLakeAnalyzer,
               reader, cowplot, dplyr, tidyr, ggplot2, zoo, purrr, beepr)

#ctd data
CTD <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/200/14/0432a298a90b2b662f26c46071f66b8a")

#flora data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.272.8
current_df <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/272/8/0359840d24028e6522f8998bd41b544e")

# metals data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.455.8
metalsdf <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/455/8/9c8c61b003923f4f03ebfe55cea8bbfd")
#removed flags for 68 as per Cece's advice

#ghgs data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.551.8
ghgs <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/551/8/454c11035c491710243cae0423efbe7b")
#not sure whether or not to remove those with flags 3 and 4
#3 = The difference between the reps are above the limit of quantification and >30% and <50% different from each other. Both replicates were retained but flagged
#4 = The difference between the reps are above the limit of quantification and >50% different from each other. Both replicates were retained but flagged

#secchi data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.198.11
secchiframe <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052")

#PAR https://portal.edirepository.org/nis/mapbrowse?packageid=edi.198.11
PAR_profiles <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/198/11/6e5a0344231de7fcebbe6dc2bed0a1c3")

#data from here https://portal.edirepository.org/nis/mapbrowse?packageid=edi.199.12
chemistry <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/199/12/a33a5283120c56e90ea414e76d5b7ddb")

#meteorological data from FCR https://portal.edirepository.org/nis/mapbrowse?packageid=edi.389.8
options(timeout = 300)
metdata <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/389/8/d4c74bbb3b86ea293e5c52136347fbb0")

#bathymetry data for BVR https://portal.edirepository.org/nis/metadataviewer?packageid=edi.1254.1
bath <- read.csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.1254.1&entityid=f7fa2a06e1229ee75ea39eb586577184")


#adding columns with total_conc max and the depth at which it occurs
DCM_BVRdata <- current_df %>% 
  filter(Reservoir == "BVR")%>%
  mutate(Date  = as_date(DateTime)) |> 
  group_by(CastID) %>%
  summarise(DCM_totalconc = max(TotalConc_ugL, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(current_df, by = "CastID") %>% # Join back the original data to keep all columns
  mutate(DCM_depth = ifelse(TotalConc_ugL == DCM_totalconc, Depth_m, NA_real_))%>% # Add DCM_depth column
  mutate(cyanosatDCM = ifelse(TotalConc_ugL == DCM_totalconc, Bluegreens_ugL, NA_real_))  #the concentration of cyanos at DCM

DCM_BVRdata <- DCM_BVRdata %>% #change this back so it doesn't say test
  group_by(CastID) %>%
  mutate(DCM_depth = ifelse(!is.na(DCM_depth), DCM_depth, NA_real_)) %>%  # Ensure DCM_depth is NA where condition didn't match
  fill(DCM_depth, .direction = "downup") %>%  # Fill NA values in DCM_depth column within each CastID
  mutate(cyanosatDCM = ifelse(!is.na(cyanosatDCM), cyanosatDCM, NA_real_)) %>%  # Ensure cyanosatDCM is NA where condition didn't match
  fill(cyanosatDCM, .direction = "downup") %>%  # Fill NA values in cyanosatDCM column within each CastID
  ungroup()%>%
  relocate(DCM_depth, .before = 3)%>%
  relocate(cyanosatDCM, .before = 3)%>%
  mutate(DOY = yday(DateTime))%>%
  mutate(Date  = as_date(DateTime)) |> 
  select(Date, DateTime, CastID, DCM_totalconc, DCM_depth, Depth_m, cyanosatDCM, Bluegreens_ugL, TotalConc_ugL,  Temp_C)

ggplot(DCM_BVRdata, aes(x = DateTime, y = DCM_depth, color = DCM_totalconc)) +
  geom_point() +
  scale_y_reverse() +  # Invert y-axis
  labs(x = "DateTime", y = "DCM_depth", color = "Concentration") +  # Labels for axes
  ggtitle("Line Graph of DCM_depth over DateTime") +  # Title of the plot
  theme_minimal()

#Adding metals
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

# Merge with DCM BVR data and keep only relevant rows
DCM_BVRwmetals <- DCM_BVRdata %>%
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  select(-SFe_mgL, -TFe_mgL, -SMn_mgL, -SCa_mgL, -TCa_mgL, -TCu_mgL, -SCu_mgL, -SBa_mgL, -TBa_mgL) %>% # Remove unnecessary columns
  filter(Depth_m %in% DCM_BVRdata$Depth_m) # Keep only rows with depths present in DCMdata
}

#Adding ghgs
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

#Adding secchi and attenuation coefficient
{
secchi <- secchiframe |>
  mutate(Date = as_date(DateTime)) |>
  group_by(Date, Reservoir, Site) |>
  summarise(Secchi_m = mean(Secchi_m, na.rm = TRUE))

BVRsecchi <- secchi |>
  filter(Reservoir == "BVR" & Site == 50)

#added Secchi
DCM_BVRwmetalsghgssecchi <- DCM_BVRwmetalsghgs|>
  left_join(BVRsecchi, by = c("Date"))|>
  group_by(Date)|>
  fill(Secchi_m, .direction = "updown")|>
  ungroup()

#calculating attenuation coefficient (alpha) and light availability

DCM_BVRwmetalsghgssecchilight <- DCM_BVRwmetalsghgssecchi |>
  mutate(alpha = 1.7/Secchi_m) |>
  mutate(light_availability_fraction = exp(-alpha * Depth_m)) |>
  mutate(light_availability_percentage = light_availability_fraction * 100)
}

# Adding PAR, DO, DOsat_percent, cond, ORP, pH (and Temp for 2017-2020)

{
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

#adding temps from 2017-2019
PAR_profiles_filtered <- PAR_profiles |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(Temp_C = mean(Temp_C, na.rm = TRUE))|>
  filter(year(Date)>2016, year(Date)<2020)

interpolated_data <- DCM_BVRdata |> 
  select(Date, Depth_m) |> 
  distinct(Date, Depth_m) |> # Get unique combinations of Date and Depth_m
  bind_rows(PAR_profiles_filtered) |> 
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

#adding pH for 2017
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

#adding pH for 2021 from the PAR profiles
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

#fixing ORP for 2017 and 2021 (grab from CTD)
#adding ORP for 2017
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


}
# Merge with DCM BVR data and keep only relevant rows


#add a column with PAR calculated from the light availability fraction (will come back to this later)
#add a column with light availability fraction using the interpolated PAR

#Adding nutrients
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

# Merge with DCM BVR data and keep only relevant rows
final_data0 <- final_data %>%
  left_join(interpolated_data, by = c("Date", "Depth_m"), relationship = "many-to-many") %>%
  select(-TN_ugL,-TP_ugL, -NH4_ugL, -NO3NO2_ugL,-SRP_ugL, -DOC_mgL, -DIC_mgL, -DC_mgL) %>% # Remove unnecessary columns
  filter(Depth_m %in% DCM_BVRdata$Depth_m)|> # Keep only rows with depths present in DCMdata
  mutate(Reservoir = "BVR")|>
  mutate(Site = 50)|>
  relocate(Reservoir, .before = 1)|>
  relocate(Site, .before = 1)

#calculating np ratio
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
#added np ratio to dataframe
final_data0 <- final_data0 %>%
  mutate(np_ratio = calculate_np_ratio(interp_TN_ugL,interp_TP_ugL))|>
  relocate(np_ratio, .before = interp_TN_ugL)
}

#visualizing metdata

metdata0 <- metdata|>
  mutate(Date = as_date(DateTime))|>
  mutate(DOY = yday(Date))|>
  relocate(DOY, .before = DateTime)|>
  relocate(Date, .before = Date)|>
  mutate(DateTime = ymd_hms(DateTime))

#function for plotting meterological variables
metplots <- function(yearz, variable, maxx = NULL){
  
  metviz <- metdata0|>
    filter(year(DateTime) == yearz) #filtering for the year specified
  
  
  ggplot(metviz, aes(x = DateTime, y = {{variable}}))+
    geom_path()+
    ggtitle(paste(deparse(substitute(variable)), yearz))+
    theme_minimal()+
    scale_y_continuous(limits = c(0, maxx))  # setting consistent y-axis limits
  
}


#Rain_Total_mm, Total rainfall over measurement interval

#precipitation make new dataframe with daily averages:

metdata0 <- metdata0 |> 
  group_by(DOY, year(DateTime)) |> 
  mutate(precip_daily = sum(Rain_Total_mm, na.rm = TRUE))|>
  ungroup()|>
  relocate(Date, .before = DateTime)|>
  relocate(precip_daily, .before = Rain_Total_mm)



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
print(b1)

print(precips)

#AirTemp_C_Average, Air temperature averaged over measurement interval. See metadata for QA/QC applied to this variable.	

b1 <- metplots(2015, AirTemp_C_Average, maxx = 50)
b2 <- metplots(2016, AirTemp_C_Average, maxx = 50)
b3 <- metplots(2017, AirTemp_C_Average, maxx = 50)
b4 <- metplots(2018, AirTemp_C_Average, maxx = 50)
b5 <- metplots(2019, AirTemp_C_Average, maxx = 50)
b6 <- metplots(2020, AirTemp_C_Average, maxx = 50)
b7 <- metplots(2021, AirTemp_C_Average, maxx = 50)
b8 <- metplots(2022, AirTemp_C_Average, maxx = 50)
b9 <- metplots(2023, AirTemp_C_Average, maxx = 50)


temps<- plot_grid(
  b1, b2, b3,
  b4, b5, b6, 
  b7, b8, b9,
  ncol = 3
)

print(temps)

#dailyaverage and dailymax for temps
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

#adding precip_daily, daily_airtempavg, maxdaily_airtemp, mindaily_airtemp to final_data0

metdata_join <- metdata0|>
  select(Date, precip_daily, daily_airtempavg, maxdaily_airtemp, mindaily_airtemp)|>
  distinct()

final_data0 <- final_data0|>
  left_join(metdata_join, by = c("Date"), relationship = "many-to-many")|>
  filter(Date %in% final_data0$Date)



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

#waterlevels

#using rLakeAnalyzer to calculate buoyancy freq, thermocline

BVRbath <- bath|>
  filter(Reservoir == "BVR")

final_data0 <- final_data0|>
  group_by(CastID)|>
  mutate(buoyancy_freq = c(buoyancy.freq(Temp_C, Depth_m), NA)) #added for padding for the last value
#need to make sure this makes sense
  
#dataframe with thermocline
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

final_datatest <- final_datatest|>
  group_by(CastID)|>
  fill(thermocline_depth, .direction = "updown")|>
  ungroup()
  


# Get system data file paths
wtr.path <- system.file('extdata', 'Sparkling.wtr', package="rLakeAnalyzer")
bth.path <- system.file('extdata', 'Sparkling.bth', package="rLakeAnalyzer")
# Load data for example lake, Sparkilng Lake, Wisconsin.
wtr = load.ts(wtr.path)
bth = load.bathy(bth.path)
## Not run:
# Generate default plot
schmidt.plot(wtr,bth)
## End(Not run)




























#preparing separate data frame so I can add a line where the chlorophyll maxima is in the heatmaps (for bluegreens currently, if I want to look at total change in the flora_heatmap function)
chlorophyll_data <- final_data0 |>
  filter(DCM_totalconc > 20)|> #choosing geater than 20 as the bloom
  filter(DCM == TRUE)|>
  filter(!(month(DateTime) == 8 & year(DateTime) == 2017 & Bluegreens_ugL < 35))|> #filter out weird drop in 2017
  filter(!(CastID == 395)) #filter out weird drop in 2016

#might be interesting to add lines for other phytos and see how they compare

flora_heatmap <- function(fp_data, reservoir, year, site, z, unitz, chlorophyll_data = NA, max_legend_value = NA)
{
  
  #subset to relevant data
  fp <- fp_data %>%
    filter(Reservoir == reservoir & year(DateTime) == year & Site == site) %>%
    select(CastID, DateTime, Depth_m, {{z}}) 
  
  #slice by depth for each reservoir
  if (reservoir == "FCR"){
    
    depths = seq(0.1, 9.3, by = 0.3)
    df.final<-data.frame()
    
    for (i in 1:length(depths)){
      
      fp_layer <- fp %>% 
        group_by(CastID) %>% 
        slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
      
      # Bind each of the data layers together.
      df.final = bind_rows(df.final, fp_layer)
      
    }
    
    
  } else if (reservoir == "BVR"){
    
    depths = seq(0.1, 10, by = 0.3)
    df.final<-data.frame()
    
    for (i in 1:length(depths)){
      
      fp_layer<-fp %>% group_by(CastID) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
      
      # Bind each of the data layers together.
      df.final = bind_rows(df.final, fp_layer)
      
    }
    
  } else if(reservoir == "CCR"){
    
    depths = seq(0.1, 20, by = 0.3)
    df.final<-data.frame()
    
    for (i in 1:length(depths)){
      
      fp_layer<-fp %>% group_by(CastID) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
      
      # Bind each of the data layers together.
      df.final = bind_rows(df.final, fp_layer)
      
    }
  } else if(reservoir == "GWR"){
    
    depths = seq(0.1, 12, by = 0.3)
    df.final<-data.frame()
    
    for (i in 1:length(depths)){
      
      fp_layer<-fp %>% group_by(CastID) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
      
      # Bind each of the data layers together.
      df.final = bind_rows(df.final, fp_layer)
      
    }
  } else if(reservoir == "SHR"){
    
    depths = seq(0.1, 30, by = 0.3)
    df.final<-data.frame()
    
    for (i in 1:length(depths)){
      
      fp_layer<-fp %>% group_by(CastID) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
      
      # Bind each of the data layers together.
      df.final = bind_rows(df.final, fp_layer)
      
    } 
    
  }
  
  #wrangle final dataframe for plotting
  # Re-arrange the data frame by date
  fp_new <- arrange(df.final, DateTime)
  
  # Round each extracted depth to the nearest 10th. 
  fp_new$Depth_m <- round(as.numeric(fp_new$Depth_m), digits = 0.5)
  
  # Convert to DOY
  fp_new$DOY <- yday(fp_new$DateTime)
  
  #trying to address error in missing values and Infs here!!!!!
  fp_new <- fp_new|>
    filter(!is.na(DOY) & !is.na(Depth_m) & !is.na(fp_new[[z]]) &
                     !is.infinite(DOY) & !is.infinite(Depth_m) & !is.infinite(fp_new[[z]]))
  
  
  fig_title <- paste(reservoir, year, "Site", site, z, sep = " ")
  
  interp <- interp(x=fp_new$DOY, y = fp_new$Depth_m, z = unlist(fp_new[z]),
                   xo = seq(min(fp_new$DOY), max(fp_new$DOY), by = .1), 
                   yo = seq(min(fp_new$Depth_m), max(fp_new$Depth_m), by = 0.01),
                   extrap = T, linear = T, duplicate = "strip")
  interp <- interp2xyz(interp, data.frame=T)
  
  # Prepare chlorophyll maxima data for line
  chlorophyll_data <- chlorophyll_data %>%
    filter(Reservoir == reservoir & year(DateTime) == year & Site == site) %>%
    mutate(DOY = yday(DateTime))|>
    filter(DOY <= max(fp_new$DOY) & DOY >= min(fp_new$DOY))
  
  p1 <- ggplot(interp, aes(x=x, y=y))+
    geom_raster(aes(fill=z))+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0), breaks = seq(1, 366, by = 30), 
                       labels = function(x) format(as.Date(x - 1, origin = paste0(year, "-01-01")), "%b")) +
    scale_fill_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +
    geom_path(data = chlorophyll_data, aes(x = DOY, y = DCM_depth, color = Bluegreens_ugL), size = 1.2) + # Color line by Bluegreens_ugL
    scale_color_gradient(low = "blue", high = "red") + # Adjust color scale as needed
    labs(x = "Day of year", y = "Depth (m)", title = fig_title,fill= unitz, color = "Bluegreens (µg/L)")+
    theme_bw()+
    theme(
      legend.text = element_text(size = 8), # Adjust text size in legend
      legend.title = element_text(size = 10), # Adjust title size in legend
      legend.key.size = unit(0.5, "cm") # Adjust the size of legend keys
    )
  
  print(p1)
  
}

#heatmaps for flora
{
  
  b1 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2014, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b2 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2015, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b3 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2016, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b4 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2017, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b5 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2018, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b6 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2019, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b7 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2020, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b8 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2021, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b9 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2022, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b10 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2023, site = 50, z = "Bluegreens_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  
  bluegreens <- plot_grid(
    b1, b2, b3,
    b4, b5,b6,
    b7, b8, b9,
    ncol = 3
  )
  
  print(bluegreens)
  
  p1 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2022, site = 50, z = "TotalConc_ugL")
  p2 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2022, site = 50, z = "BrownAlgae_ugL")
  p3 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2022, site = 50, z = "GreenAlgae_ugL")
  p4 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2022, site = 50, z = "Bluegreens_ugL")
  p5 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2022, site = 50, z = "MixedAlgae_ugL")
  p6 <- flora_heatmap(fp_data = current_df, reservoir = "BVR", year = 2022, site = 50, z = "YellowSubstances_ugL")
  
  final_plot <- plot_grid(
    p1, p2, p3,
    p4, p5, p6,
    ncol = 3  # Specify the number of columns
  )
  
}

#heatmaps for SFe_mgL
{

dataforheatmap <- final_data0 |>
  filter(!is.na(interp_SFe_mgL))  # Remove rows with NA in interp_SFe_mgL


b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_SFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)

soluble_iron <- plot_grid(
  b1, b2, b3,
  b4, b5, b6,
  b7, b8, b9,
  ncol = 3  # Specify the number of columns
)

print(soluble_iron)
}

#heatmaps for interp_TFe_mgL
{
  
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_TFe_mgL))  # Remove rows with NA in interp_SFe_mgL
  
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_TFe_mgL", unitz = "mgL", chlorophyll_data, max_legend_value = 35)
  
  total_iron <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7, b8, b9,
    ncol = 3  # Specify the number of columns
  )
  
  print(total_iron)
}
#heatmaps for SMn_mgL
{
  
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_SMn_mgL))  # Remove rows with NA in interp_SMn_mgL
  
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_SMn_mgL", unitz = "mgL", chlorophyll_data)
  
  SMn_mgL_plot <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7, b8, b9,
    ncol = 3  # Specify the number of columns
  )
  
  print(SMn_mgL_plot)
}
#heatmaps for interp_SCa_mgL
#error
{
  
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_SCa_mgL))  # Remove rows with NA in interp_SMn_mgL
  
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_SCa_mgL", unitz = "mgL", chlorophyll_data)
  
  interp_SCa_mgL <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7, b8, b9,
    ncol = 3  # Specify the number of columns
  )
  
  print(interp_SCa_mgL)
}
#more metals to add

#heatmaps for interp_CO2_umolL
#no data for 2014
#all data colinear for b2 (will look into this)
{
dataforheatmap <- final_data0 |>
  filter(!is.na(interp_CO2_umolL))  # Remove rows with NA in interp_SFe_mgL

#b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_CO2_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)
b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_CO2_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)
b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_CO2_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)
b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_CO2_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)
b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_CO2_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)
b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_CO2_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)
b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_CO2_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)
b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_CO2_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)
b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_CO2_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)

CO2_plots <- plot_grid(
  b3, b4, b5,
  b6, b7, b8,
  b9, b10,
  ncol = 3
)

print(CO2_plots)
}

#heatmaps for interp_CH4_umolL (no data for 2014)
{
dataforheatmap <- final_data0 |>
  filter(!is.na(interp_CH4_umolL))  # Remove rows with NA in interp_SFe_mgL

#not sure why not working b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_CH4_umolL", unitz = "µmol/L", chlorophyll_data, max_legend_value = 700)
b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_CH4_umolL", unitz = "µmol/L", chlorophyll_data)
b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_CH4_umolL", unitz = "µmol/L", chlorophyll_data)
b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_CH4_umolL", unitz = "µmol/L", chlorophyll_data)
b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_CH4_umolL", unitz = "µmol/L", chlorophyll_data)
b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_CH4_umolL", unitz = "µmol/L", chlorophyll_data)
b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_CH4_umolL", unitz = "µmol/L", chlorophyll_data)
b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_CH4_umolL", unitz = "µmol/L", chlorophyll_data)
b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_CH4_umolL", unitz = "µmol/L", chlorophyll_data)

methane_plots <- plot_grid(
  b3, b4, b5,
  b6, b7, b8,
  b9, b10,
  ncol = 3
)

print(methane_plots)
}

#heatmaps for interp_PAR_umolm2s
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_PAR_umolm2s))|>  # Remove rows with NA in interp_PAR_umolm2s
    filter(month(DateTime) < 11)

 
#only one day b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data, max_legend_value = 50)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data, max_legend_value = 50)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data, max_legend_value = 50)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data, max_legend_value = 50)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data, max_legend_value = 50)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data, max_legend_value = 50)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data, max_legend_value = 50)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data, max_legend_value = 50)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_PAR_umolm2s", unitz = "µmol/m2s", chlorophyll_data, max_legend_value = 50)
  
  interp_PAR_plots <- plot_grid(
    b2, b3, b4,
    b5, b6, b7,
    b8, b9, b10,
    ncol = 3
  )
  
  print(interp_PAR_plots)
}

#heatmaps for interp_DO_mgL
#b1 and b10 not working

{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_DO_mgL))|>  # Remove rows with NA in interp_DOC_mgL
    filter(month(DateTime) < 11)
  
  #b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_DO_mgL", unitz = "mg/L", chlorophyll_data)
  
  interp_DO_plots <- plot_grid(
    b2, b3, b4,
    b5,b6, b7,
    b8, b9,
    ncol = 3
  )
  
  print(interp_DO_plots)
}

#heatmaps for interp_Cond_uScm

{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_Cond_uScm))|>  # Remove rows with NA in interp_Cond_uScm
    filter(month(DateTime) < 11)
  
#no data b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)
#colinear  #b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_Cond_uScm", unitz = "uScm", chlorophyll_data)
  
  interp_Cond_uScm <- plot_grid(
    b2, b3, b5,
    b6, b7, b8,
    b9, b10,
    ncol = 3
  )
  
  print(interp_Cond_uScm)
}

#visualization for interp_ORP_mvV

{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_ORP_mV))|> 
    filter(month(DateTime) < 11) 
  
  looking <- dataforheatmap|>
    filter(year(DateTime) == 2021)|>
    select(DateTime, Depth_m, interp_ORP_mV)
  
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_ORP_mV", unitz = "mV", chlorophyll_data, max_legend_value = 400)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_ORP_mV", unitz = "mV", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_ORP_mV", unitz = "mV", chlorophyll_data, max_legend_value = 400)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_ORP_mV", unitz = "mV", chlorophyll_data, max_legend_value = 400)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_ORP_mV", unitz = "mV", chlorophyll_data, max_legend_value = 400)
#there just is no data b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_ORP_mV", unitz = "mV", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_ORP_mV", unitz = "mV", chlorophyll_data, max_legend_value = 400)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_ORP_mV", unitz = "mV", chlorophyll_data, max_legend_value = 400)
  
  ORP <- plot_grid(
    b3, b4, b5,
    b6, b7, b9,
    b10,
    ncol = 3
  )
  
  print(ORP)
}

#heatmaps for Temp_C
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(Temp_C))|>  # Remove rows with NA in Temp_C
    filter(month(DateTime) < 11)
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "Temp_C", unitz = "C", chlorophyll_data)
  
  temp_plots <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7,b8, b9,
    ncol = 3
  )
  
  print(temp_plots)
}

#interp_pH


{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_pH))|>  # Remove rows with NA in interp_SFe_mgL
    filter(month(DateTime) < 11) #doing before November because an NP ratio greater than 6000 is crazy, so I am excluding it

  looking <- dataforheatmap|>
    filter(year(Date) == 2017, !is.na(interp_pH))|>
    select(Date, Depth_m, interp_pH)
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_pH", unitz = "pH", chlorophyll_data, max_legend_value = 8)
  
  pH <- plot_grid(
    b1, b2, b3, 
    b4, b5, b6, 
    b7, b8, b9,
    ncol = 3
  )
  
  print(pH)
}

#heatmaps for np_ratio
{
dataforheatmap <- final_data0 |>
  filter(!is.na(np_ratio))|>  # Remove rows with NA in interp_SFe_mgL
    filter(month(DateTime) < 11) #doing before November because an NP ratio greater than 6000 is crazy, so I am excluding it
  
b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data)
b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data)
b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data)
b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data)
b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data)  
b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data, max_legend_value = 200)
b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data)
b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data)
b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data)
b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "np_ratio", unitz = "NP ratio", chlorophyll_data)

np_ratio_plots <- plot_grid(
  b1, b2, b3, 
  b4, b5, b6,
  b7, b8, b9,
  ncol = 3
)

print(np_ratio_plots)
}

#heatmaps for interp_TN_ugL
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_TN_ugL))|>  # Remove rows with NA in interp_SFe_mgL
    filter(month(DateTime) < 11) #doing before November because an NP ratio greater than 6000 is crazy, so I am excluding it
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_TN_ugL", unitz = "µg/L", chlorophyll_data)
  
  interp_TN_plots <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7, b8, b9,
    ncol = 3
  )
  
  print(interp_TN_plots)
  }

#heatmaps for interp_NH4_ugL
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_NH4_ugL))|>  # Remove rows with NA in interp_SFe_mgL
    filter(month(DateTime) < 11) #doing before November because an NP ratio greater than 6000 is crazy, so I am excluding it
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_NH4_ugL", unitz = "µg/L", chlorophyll_data)
  
  interp_NH4_plots <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7, b8, b9,
    ncol = 3
  )
  
  print(interp_NH4_plots)
}

#heatmaps for interp_NO3NO2_ugL
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_NO3NO2_ugL))|>  # Remove rows with NA in interp_SFe_mgL
    filter(month(DateTime) < 11 & month(DateTime) > 4) 
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_NO3NO2_ugL", unitz = "µg/L", chlorophyll_data, max_legend_value = 15)
  
  interp_NO3NO2_plots <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7, b8, b9,
    ncol = 3
  )
  
  print(interp_NO3NO2_plots)
}

#heatmaps for interp_SRP_ugL
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_SRP_ugL))|>  # Remove rows with NA in interp_SFe_mgL
    filter(month(DateTime) < 11) #doing before November because an NP ratio greater than 6000 is crazy, so I am excluding it
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_SRP_ugL", unitz = "µg/L", chlorophyll_data)
  
  interp_SRP_plots <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7, b8, b9,
    ncol = 3
  )
  
  print(interp_SRP_plots)
}

#heatmaps for interp_TP_ugL
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_TP_ugL))|>  # Remove rows with NA in interp_SFe_mgL
    filter(month(DateTime) < 11) #doing before November because an NP ratio greater than 6000 is crazy, so I am excluding it
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_TP_ugL", unitz = "µg/L", chlorophyll_data)
  
  interp_TP_plots <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7, b8, b9,
    ncol = 3
  )
  
  print(interp_TP_plots)
}

#heatmaps for interp_DOC_mgL
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_DOC_mgL))|>  # Remove rows with NA in interp_DOC_mgL
    filter(month(DateTime) < 11)
  
  b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)
  b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)
  b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)
  b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_DOC_mgL", unitz = "mg/L", chlorophyll_data)
  
  interp_DOC_plots <- plot_grid(
    b1, b2, b3,
    b4, b5, b6,
    b7, b8, b9,
    ncol = 3
  )
  
  print(interp_DOC_plots)
}

#heatmaps for interp_DIC_mgL
#ERROR HEREEEEE
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_DIC_mgL))|>  # Remove rows with NA in interp_DOC_mgL
    filter(month(DateTime) < 11)
  
  #b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)
  #b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)
  #b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)
  #b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_DIC_mgL", unitz = "mg/L", chlorophyll_data)
  
  interp_DIC_plots <-  plot_grid(
    b5, b6, b7, 
    b8, b9, b10,
    ncol = 3
  )
  
  print(interp_DIC_plots)
}


#heatmaps for interp_DC_mgL
#ERROR HERE AS WELL
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_DC_mgL))|>  # Remove rows with NA in interp_DOC_mgL
    filter(month(DateTime) < 11)
  
  #b1 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2014, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)
  #b2 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2015, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)
  #b3 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2016, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)
  #b4 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2017, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)
  b5 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2018, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)  
  b6 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2019, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)
  b7 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2020, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)
  b8 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2021, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)
  b9 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2022, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)
  b10 <- flora_heatmap(fp_data = dataforheatmap, reservoir = "BVR", year = 2023, site = 50, z = "interp_DC_mgL", unitz = "mg/L", chlorophyll_data)
  
  interp_DC_plots <- plot_grid(
    b5, b6, b7, 
    b8, b9, b10,
    ncol = 3
  )
  
  print(interp_DC_plots)
}


#looking at weird drop really quick will delete later
#drop in 2015 mid-may
#drop in July from about 6 ft to 9.5ish feet


#drop in 2017 a couple of days into July from above 5.0 ft down to the bottom
seventeendrop <- current_df|>
  filter(year(DateTime) == 2017 & Reservoir == "BVR" & month(DateTime)== 8 & day(DateTime) == 10)
#2017-08-10 it's only 30ish ug throughout at 9.90. i think i should remove this
#it is flagged 3 throughout for RFU 590nm (will look into this)

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



