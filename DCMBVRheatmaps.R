#Heatmaps for BVR_DCM
#Maria Popescu
#alot of help writing the function from Mary



phytos <- read.csv("./phytos.csv")



#need to change these values to reflect Bluegreens not just chla (which it currently is)
# Chlorophyll data for the line visualizing chl max on heatmaps

chlorophyll_data <- final_data0 %>%
  select(Date, Depth_m, Bluegreens_ugL)|> #change to chl here if you want
  mutate(DayOfYear = yday(Date))|>
  group_by(Date) %>%
  slice(which.max(Bluegreens_ugL)) %>%
  ungroup()|>
  mutate(Reservoir = "BVR")|>
  mutate(DateTime = Date)|> #no actual time
  filter(DayOfYear > 133, DayOfYear < 285, Bluegreens_ugL >20)|>
  filter(!(month(Date) == 8 & year(Date) == 2017 & Bluegreens_ugL < 35)) #weird drop here


####heatmaps####
# Preparing separate data frame so I can add a line where the thermocline is (removed line from graph but can add it back in)
thermocline_df <- final_data0|>
  mutate(Reservoir = "BVR")|>
  mutate(DateTime = Date)|> #this doesn't actually have the time but I need it for the rest of the script
  select(Reservoir, Date, thermocline_depth, Temp_C, Depth_m, Site)

Photic_zone<- final_data0|>
  select(Date, PZ)|>
  group_by(Date)|>
  mutate(PZ = mean(PZ))|>
  mutate(DOY = yday(Date))
  

#might be interesting to add lines for other phytos and see how they compare

flora_heatmap <- function(fp_data, year, site, z, unitz, chlorophyll_data = NA, max_legend_value = NA)
{
  
  #subset to relevant data
  fp <- fp_data %>%
    filter(year(DateTime) == year) %>%
    select(CastID, DateTime, Depth_m, {{z}})
    
  
  #slice by depth for each reservoir
    depths = seq(0.1, 10, by = 0.3)
    df.final<-data.frame()
    
    for (i in 1:length(depths)){
      
      fp_layer<-fp %>% group_by(CastID) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
      
      # Bind each of the data layers together.
      df.final = bind_rows(df.final, fp_layer)
      
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
  
  
  fig_title <- paste("BVR", year, "Site", site, z, sep = " ")
  
  interp <- interp(x=fp_new$DOY, y = fp_new$Depth_m, z = unlist(fp_new[z]),
                   xo = seq(min(fp_new$DOY), max(fp_new$DOY), by = .001), #this is what i'm playing with right now
                   yo = seq(min(fp_new$Depth_m), max(fp_new$Depth_m), by = 0.01),
                   extrap = T, linear = T, duplicate = "strip")
  interp <- interp2xyz(interp, data.frame=T)
  
  # Prepare chlorophyll maxima data for line
  #chlorophyll_data <- chlorophyll_data %>%
  #  filter(year(DateTime) == year & site == site) %>%
  #  mutate(DOY = yday(DateTime))|>
  #  filter(DOY <= max(fp_new$DOY) & DOY >= min(fp_new$DOY))
  
  #Photic_zone <- Photic_zone %>%
  #  filter(year(Date) == year) %>%
  #  mutate(DOY = yday(Date))

  p1 <- ggplot(interp, aes(x=x, y=y))+
    geom_raster(aes(fill=z))+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0), breaks = seq(1, 366, by = 30), 
                       labels = function(x) format(as.Date(x - 1, origin = paste0(year, "-01-01")), "%b")) +
    scale_fill_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +
    #geom_path(data = chlorophyll_data, aes(x = DOY, y = Depth_m, color = Bluegreens_ugL), size = 1.2) + # Color line by Bluegreens_ugL
    #scale_color_gradient(low = "blue", high = "red") + # Adjust color scale as needed
    labs(x = "Day of year", y = "Depth (m)", title = fig_title,fill= unitz, color = "Bluegreens (µg/L)")+
    theme_bw()+
    theme(
      legend.text = element_text(size = 8), # Adjust text size in legend
      legend.title = element_text(size = 10), # Adjust title size in legend
      legend.key.size = unit(0.5, "cm") # Adjust the size of legend keys
    )
  
  print(p1)
  
}

#### flora ####
{
  
  b1 <- flora_heatmap(fp_data = phytos, year = 2014, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  b2 <- flora_heatmap(fp_data = current_df, year = 2015, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  b3 <- flora_heatmap(fp_data = current_df, year = 2016, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  b4 <- flora_heatmap(fp_data = current_df, year = 2017, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  b5 <- flora_heatmap(fp_data = current_df, year = 2018, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  b6 <- flora_heatmap(fp_data = current_df, year = 2019, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  b7 <- flora_heatmap(fp_data = current_df, year = 2020, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  b8 <- flora_heatmap(fp_data = current_df, year = 2021, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  b9 <- flora_heatmap(fp_data = current_df, year = 2022, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  b10 <- flora_heatmap(fp_data = current_df, year = 2023, site = 50, z = "Bluegreens_ugL", unitz = "ug/L")
  
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

#### SFe_mgL  ####
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

#### TFe_mgL  ####
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
#### SMn_mgL  ####
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
#### SCa_mgL  ####
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

#### CO2_umolL  ####
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

#### CH4_umolL ####
#(no data for 2014)
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

#### PAR_umolm2s ####
{
  dataforheatmap <- final_data0 |>
    filter(!is.na(interp_PAR_umolm2s))|>  # Remove rows with NA in interp_PAR_umolm2s
    filter(month(Date) < 11)
  
  
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

#### DO_mgL ####
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

#### Cond_uScm  ####

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

#### ORP_mvV  ####

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

#### Temp_C  ####
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

#### pH  ####
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
