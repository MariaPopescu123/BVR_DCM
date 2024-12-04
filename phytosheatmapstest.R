#phytoheatmaps test
#Heatmaps for BVR_DCM
#Maria Popescu
#alot of help writing the function from Mary


pacman::p_load(tidyverse, lubridate, akima, reshape2, 
               gridExtra, grid, colorRamps, RColorBrewer, rLakeAnalyzer,
               reader, cowplot, dplyr, tidyr, ggplot2, zoo, purrr, beepr, forecast, ggthemes)

phytos_for_heatmaps <- read.csv("./phytos_for_heatmaps")


#if you just want to look at the dates within the timeframe chosen use heatmap_data.csv
#if you want to have the data for all dates use ./final_data_alldates.csv

heatmap_data<- phytos_for_heatmaps|>
  mutate(Date = as.Date(Date, format = "%Y-%m-%d"))  # Adjust format if necessary


#need to change these values to reflect Bluegreens not just chla (which it currently is)
# Chlorophyll data for the line visualizing chl max on heatmaps

chlorophyll_data <- heatmap_data %>%
  select(Date, Depth_m, TotalConc_ugL)|>
  mutate(DayOfYear = yday(Date))|>
  group_by(Date) %>%
  slice(which.max(TotalConc_ugL)) %>%
  ungroup()|>
  mutate(Reservoir = "BVR")|>
  mutate(Date = Date)|> #no actual time
  filter(TotalConc_ugL >20)|>
  filter(!(month(Date) == 8 & year(Date) == 2017 & TotalConc_ugL < 35)) #weird drop here


####heatmaps####

#might be interesting to add lines for other phytos and see how they compare

flora_heatmap <- function(fp_data, year, site, z, unitz, chlorophyll_data = NA, max_legend_value = NA)
{
  
  #subset to relevant data
  fp <- fp_data %>%
    mutate(Date = as.Date(Date,format = "%Y-%m-%d"))|>
    filter(lubridate::year(Date) ==year) %>%
    select(CastID, Date, Depth_m, {{z}}) 
  
  
  depths = seq(0.1, 10, by = 0.3)
  df.final<-data.frame()
  
  for (i in 1:length(depths)){
    
    fp_layer<-fp %>% group_by(CastID) %>% slice(which.min(abs(as.numeric(Depth_m) - depths[i])))
    
    # Bind each of the data layers together.
    df.final = bind_rows(df.final, fp_layer)
    
    
  }
  
  #wrangle final dataframe for plotting
  # Re-arrange the data frame by date
  fp_new <- arrange(df.final, Date)
  
  # Round each extracted depth to the nearest 10th. 
  fp_new$Depth_m <- round(as.numeric(fp_new$Depth_m), digits = 0.5)
  
  # Convert to DOY
  fp_new$DOY <- yday(fp_new$Date)
  
  #trying to address error in missing values and Infs here!!!!!
  fp_new <- fp_new|>
    filter(!is.na(DOY) & !is.na(Depth_m) & !is.na(fp_new[[z]]) &
             !is.infinite(DOY) & !is.infinite(Depth_m) & !is.infinite(fp_new[[z]]))
  
  
  fig_title <- paste("BVR", year, "Site", site, z, sep = " ")
  
  interp <- interp(x=fp_new$DOY, y = fp_new$Depth_m, z = unlist(fp_new[z]),
                   xo = seq(min(fp_new$DOY), max(fp_new$DOY), by = .1), 
                   yo = seq(min(fp_new$Depth_m), max(fp_new$Depth_m), by = 0.01),
                   extrap = T, linear = T, duplicate = "strip")
  interp <- interp2xyz(interp, data.frame=T)
  
  # Prepare chlorophyll maxima data for line
  chlorophyll_data <- chlorophyll_data %>%
    filter(year(Date) == year & site == site) %>%
    mutate(DOY = yday(Date))|>
    filter(DOY <= max(fp_new$DOY) & DOY >= min(fp_new$DOY))
  
  p1 <- ggplot(interp, aes(x=x, y=y))+
    geom_raster(aes(fill=z))+
    scale_y_reverse(expand = c(0,0))+
    scale_x_continuous(expand = c(0, 0), breaks = seq(1, 366, by = 30), 
                       labels = function(x) format(as.Date(x - 1, origin = paste0(year, "-01-01")), "%b")) +
    scale_fill_gradientn(colours = blue2green2red(60), na.value = "gray", limits = c(NA, max_legend_value)) +
    geom_path(data = chlorophyll_data, aes(x = DOY, y = Depth_m, color = TotalConc_ugL), size = 1.2) + # Color line by TotalConc_ugL
    scale_color_gradient(low = "blue", high = "red") + # Adjust color scale as needed
    labs(x = "Day of year", y = "Depth (m)", title = fig_title,fill= unitz, color = "TotalConc (µg/L)")+
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
  
  b1 <- flora_heatmap(fp_data = heatmap_data, year = 2014, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b2 <- flora_heatmap(fp_data = heatmap_data,  year = 2015, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b3 <- flora_heatmap(fp_data = heatmap_data,  year = 2016, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b4 <- flora_heatmap(fp_data = heatmap_data,  year = 2017, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b5 <- flora_heatmap(fp_data = heatmap_data,  year = 2018, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b6 <- flora_heatmap(fp_data = heatmap_data,  year = 2019, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b7 <- flora_heatmap(fp_data = heatmap_data,  year = 2020, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b8 <- flora_heatmap(fp_data = heatmap_data,  year = 2021, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b9 <- flora_heatmap(fp_data = heatmap_data,  year = 2022, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  b10 <- flora_heatmap(fp_data = heatmap_data,  year = 2023, site = 50, z = "TotalConc_ugL", unitz = "ug/L", chlorophyll_data = chlorophyll_data)
  
  totals_heatmaps <- plot_grid(
    b1, b2, b3, b4, b5,
    b6, b7, b8, b9, b10,
    ncol = 5
  )
  
  print(totals_heatmaps)
  
  p1 <- flora_heatmap(fp_data = heatmap_data,  year = 2022, site = 50, z = "TotalConc_ugL")
  p2 <- flora_heatmap(fp_data = heatmap_data,  year = 2022, site = 50, z = "BrownAlgae_ugL")
  p3 <- flora_heatmap(fp_data = heatmap_data,  year = 2022, site = 50, z = "GreenAlgae_ugL")
  p4 <- flora_heatmap(fp_data = heatmap_data,  year = 2022, site = 50, z = "TotalConc_ugL")
  p5 <- flora_heatmap(fp_data = heatmap_data,  year = 2022, site = 50, z = "MixedAlgae_ugL")
  p6 <- flora_heatmap(fp_data = heatmap_data,  year = 2022, site = 50, z = "YellowSubstances_ugL")
  
  final_plot <- plot_grid(
    p1, p2, p3,
    p4, p5, p6,
    ncol = 3  # Specify the number of columns
  )
  
}
