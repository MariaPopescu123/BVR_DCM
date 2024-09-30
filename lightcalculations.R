#####PAR calculations####
CTD <- read.csv("./CTD.csv")

PAR_profiles <- read.csv("PAR_profiles.csv")


CTDfiltered <- CTD|>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  group_by(Date, Depth_m)|>
  summarise(DO_mgL = mean(DO_mgL, na.rm = TRUE),
            PAR_umolm2s = mean(PAR_umolm2s, na.rm = TRUE), 
            DOsat_percent = mean(DOsat_percent, na.rm = TRUE),
            Cond_uScm = mean(Cond_uScm, na.rm = TRUE),
            ORP_mV = mean(ORP_mV, na.rm = TRUE),
            pH = mean(pH, na.rm = TRUE))


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


plot_dat <- final_PAR %>%
  mutate(Year = year(Date), 
         DayOfYear = yday(Date))|> # Extract year and day of the year
  select(Date, Year, DayOfYear, interp_PAR_umolm2s, Depth_m)

redpoints <- plot_dat|>
  filter(DayOfYear>133, DayOfYear<286)


# Plot: x-axis is DayOfYear, y-axis is Year, with a line and highlighted points
ggplot(plot_dat, aes(x = DayOfYear, y = as.factor(Year), group = Year)) +
  geom_line() +  # Line for each year
  geom_point() +  # Data points
  theme_bw() +
  labs(x = "Day of Year", y = "Year", title = "combined_df2") +
  geom_point(data = redpoints, aes(x = DayOfYear, y = as.factor(Year)), 
             color = "red", size = 3) +  # Highlight max points in red
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +  # Adjust x-axis breaks
  theme(panel.grid.minor = element_blank())  # Optional: remove minor grid lines

