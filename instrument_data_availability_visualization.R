#Instrument data availability (raw)

#must have all the dataframes loaded in already from DCMBVRdatawrangling.R


#flora data
flora_check <- current_df|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))
variables <- c("Totals_ugL")


#metals
metals_check <- metalsdf|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))
variables <- c("SFe_mgL", "TFe_mgL", "SMn_mgL", "SCa_mgL",
               "TCa_mgL", "TCu_mgL", "SCu_mgL", "SBa_mgL", "TBa_mgL")

#ghgs
ghgs_check <- ghgs|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))
variables <- c("CO2_umolL", "CH4_umolL")



#secchi
secchi_check <- secchiframe |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))
variables <- c("Secchi_m")


#ysi
ysi_check <- ysi_profiles|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))
variables <- c("DO_mgL", "PAR_umolm2s", "DOsat_percent", "Cond_uScm", "ORP_mV", "pH", "Temp_C")

#CTD
CTD_check <- CTD|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))
variables <- c("Temp_C", "DO_mgL", "DOsat_percent", 
               "Cond_uScm", "SpCond_uScm", "Chla_ugL", "Turbidity_NTU", 
               "pH", "ORP_mV", "PAR_umolm2s", "CDOM_ugL")

data_availability(CTD_check, variables)


#chem
chem_check <- chemistry|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))
variables <- c("TN_ugL", "TP_ugL", "NH4_ugL", "NO3NO2_ugL", "SRP_ugL", 
               "DOC_mgL", "DIC_mgL", "DC_mgL", "DN_mgL")

#met
metcheck <- metdata|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))

#waterlevel data from water level edi package 
wtrlvl_check <- wtrlvl|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))

#waterlevel data using the pressure sensor (platform data) 
BVR_check <- BVRplatform|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime)) 