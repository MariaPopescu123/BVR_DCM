
####to see dates that data is available for
data_available <- function(data, title_text, variable){

plot_dat <- data %>%
  mutate(Year = year(Date), 
         DayOfYear = yday(Date))|> # Extract year and day of the year
  select(Date, Year, DayOfYear, {{ variable }})

redpoints <- plot_dat|>
  filter(DayOfYear>133, DayOfYear<286)


# Plot: x-axis is DayOfYear, y-axis is Year, with a line and highlighted points
ggplot(plot_dat, aes(x = DayOfYear, y = as.factor(Year), group = Year)) +
  geom_line() +  # Line for each year
  geom_point() +  # Data points
  theme_bw() +
  labs(x = "Day of Year", y = "Year", title = {{title_text}}) +
  geom_point(data = redpoints, aes(x = DayOfYear, y = as.factor(Year)), 
             color = "red", size = 3) +  # Highlight max points in red
  scale_x_continuous(breaks = seq(1, 365, by = 30)) +  # Adjust x-axis breaks
  theme(panel.grid.minor = element_blank())  # Optional: remove minor grid lines
}




