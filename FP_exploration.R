# L1 FP data exploration 
# Author: Mary Lofton
# Date: 08JUL24

library(tidyverse)
library(lubridate)

#flora data raw 2024
dat <- read_csv("https://raw.githubusercontent.com/CareyLabVT/Reservoirs/master/Data/DataNotYetUploadedToEDI/Raw_fluoroprobe/fluoroprobe_L1.csv")

#flora data all years

dat <- read.csv("./current_df.csv")

length(unique(dat$CastID)) 

casts <- dat %>%
  mutate(Date = date(DateTime)) %>%
  select(Reservoir, Date, Site) %>%
  distinct() 

# Look at casts post-servicing

# limit to casts you actually want
#change dates here
make_plot <- function(desired_date) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(lubridate)
  
  # Filter sample data based on the desired date
  sample_dat <- casts %>%
    filter(Reservoir == "BVR" & Site == 50 & Date %in% c(desired_date)) %>%
    mutate(SampleDepth = c(3.6), 
           Date = factor(Date, levels = c(desired_date)),
           Date = as.Date(as.character(Date)))  # Convert factor back to Date
  
  # Process the data for plotting
  plot_dat <- dat %>%
    mutate(Date = date(DateTime)) %>%
    filter(!month(Date) == 1) %>%
    group_by(Reservoir, Date, Site) %>%
    mutate(FacetID = paste(CastID, Reservoir, Site, Date, sep = " ")) %>%
    ungroup() %>%
    right_join(sample_dat, by = c("Reservoir", "Date")) %>%
    select(CastID, FacetID, GreenAlgae_ugL, Bluegreens_ugL, BrownAlgae_ugL, MixedAlgae_ugL, YellowSubstances_ugL,
           TotalConc_ugL, Depth_m, SampleDepth) %>%
    pivot_longer(cols = GreenAlgae_ugL:TotalConc_ugL, names_to = "var", values_to = "ugL")
  
  # Generate the plot
  plot_casts <- ggplot(plot_dat, aes(x = ugL, y = Depth_m, group = var, color = var)) +
    geom_path() +
    scale_y_reverse() +
    theme_bw() +
    facet_wrap(facets = vars(FacetID)) +
    xlab("Micrograms per liter") +
    ylab("Depth (m)") +
   # geom_hline(aes(yintercept = SampleDepth)) +
    scale_color_discrete(name = "Variable") +
    scale_linetype_discrete(na.translate = FALSE)
  
  return(plot_casts)
}

#dates that I want: "2014-07-23" "2015-06-25" "2016-06-30" "2017-08-17"
#"2018-08-16" "2019-05-23" "2020-08-20" "2021-08-09" "2022-08-18" "2023-07-17"#change date to see correlations for the singular day that max was the biggest

b1 <- make_plot("2014-07-23")
b2 <- make_plot("2015-06-25")
b3 <- make_plot("2016-06-30")
b4 <- make_plot("2017-08-17")
b5 <- make_plot("2018-08-16")
b6 <- make_plot("2019-05-23")
b7 <- make_plot("2020-08-20")
b8 <- make_plot("2021-08-09")
b9 <- make_plot("2022-08-18")
b10 <- make_plot("2023-07-17")

print(b10)

allfloraplots <- plot_grid(
  b1, b2, b3, b4, b5,
  b6, b7, b8, b9, b10,
  ncol = 5
)

print(allfloraplots)


#ggsave("./Desktop/FP_exploration/FP_casts_2024.png", device = "png",
#       height = 6, width = 10.5, units = "in")  # you will probably want to change the dimensions

## made edits through here - after that you are on your own :-)

group_biomass <- plot_dat %>%
  group_by(FacetID, var) %>%
  slice(which.min(abs(Depth_m - SampleDepth))) %>%
  mutate(DepthDiff = abs(SampleDepth - Depth_m)) %>%
  filter(CastID %in% c(1167, 1171, 1168, 1179))

group_biomass_plot <- ggplot(data = group_biomass, aes(x = var, y = ugL, color = var, fill = var))+
  geom_bar(stat = "identity")+
  #facet_wrap(facets = vars(FacetID), nrow = 1)+
  theme_bw()+
  scale_color_manual(name = "Variable", values = c("cyan4","chocolate4","darkolivegreen4","gray","black","goldenrod"))+
  scale_fill_manual(name = "Variable", values = c("cyan4","chocolate4","darkolivegreen4","gray","black","goldenrod"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  geom_text(aes(x = "Bluegreens_ugL", y = 65, label = round(DepthDiff,2)), color = "black")+
  ggtitle("Text inside each subplot is distance between FP measurement and phyto sample in meters")+
  xlab("")
group_biomass_plot
#ggsave("./Desktop/FP_exploration/FP_group_biomass_at_sample_depth_2024.png", device = "png",
#       height = 3, width = 12, units = "in") 
















