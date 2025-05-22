#function for summarizing variables 

library(dplyr)
library(lubridate)
library(rlang)
library(purrr)

#example of use 
#Date must be edited to be just Date
#variables <- c("SFe_mgL", "TFe_mgL", "SMn_mgL")
#summary_output <- summarize_data_by_week(metalsdf, variables)
weekly_sum_variables <- function(df, variables) {
  # Step 0: Create Date, Week, and Year columns
  df_prepped <- df |>
    mutate(
      Week = lubridate::week(Date),
      Year = lubridate::year(Date)
    ) |>
    select(Date, Depth_m, Week, Year, all_of(variables))
  
  # Step 1: Identify valid Week-Year combos with depth range > 4 m
  valid_weeks <- df_prepped |>
    group_by(Year, Week, Date) |>
    summarise(
      depth_range = max(Depth_m, na.rm = TRUE) - min(Depth_m, na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(depth_range > 4) |>
    distinct(Year, Week)
  
  # Step 2: Filter for valid Week-Year combinations
  df_filtered <- df_prepped |>
    semi_join(valid_weeks, by = c("Year", "Week")) |>
    group_by(Year, Week)
  
  # Step 3: Summarize each variable separately
  summary_list <- lapply(variables, function(var) {
    var_sym <- sym(var)
    
    df_filtered |>
      summarise(
        !!paste0("depth_", var, "_max") := Depth_m[which.max(!!var_sym)],
        !!paste0("depth_", var, "_min") := Depth_m[which.min(!!var_sym)],
        !!paste0(var, "_max_val") := max(!!var_sym, na.rm = TRUE),
        !!paste0(var, "_min_val") := min(!!var_sym, na.rm = TRUE),
        !!paste0(var, "_range") := max(!!var_sym, na.rm = TRUE) - min(!!var_sym, na.rm = TRUE),
        mean_val = mean(!!var_sym, na.rm = TRUE),
        !!paste0("depth_mean_", var) := Depth_m[which.min(abs(!!var_sym - mean(!!var_sym, na.rm = TRUE)))],
        .groups = "drop"
      ) |>
      select(-mean_val)  # Drop helper column
  })
  
  # Step 4: Combine all summaries
  combined_summary <- reduce(summary_list, left_join, by = c("Year", "Week"))
  
  # Step 5: Add earliest date for each Week-Year
  week_dates <- df_filtered |>
    summarise(Date = min(Date), .groups = "drop")
  
  combined_summary |>
    left_join(week_dates, by = c("Year", "Week")) |>
    select(Date, everything())
}