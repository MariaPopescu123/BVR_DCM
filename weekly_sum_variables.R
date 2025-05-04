#function for summarizing variables 

library(dplyr)
library(lubridate)
library(rlang)
library(purrr)

#example of use 
#variables <- c("SFe_mgL", "TFe_mgL", "SMn_mgL")
#summary_output <- summarize_data_by_week(metalsdf, variables)
weekly_sum_variables <- function(df, variables) {
  df_prepped <- df |>
    mutate(
      Date = if ("DateTime" %in% names(data)) as.Date(DateTime) else Date
    ) |>
    mutate(Week = lubridate::week(Date)) |>
    group_by(Year, Week) |>
    select(Date, Depth_m, all_of(variables)) |>
    mutate(Week = week(Date), Year = year(Date))
  
  # Step 1: Identify valid Week-Year combos with at least one Date having depth range > 4 m
  valid_weeks <- df_prepped |>
    group_by(Year, Week, Date) |>
    summarise(depth_range = max(Depth_m, na.rm = TRUE) - min(Depth_m, na.rm = TRUE), .groups = "drop") |>
    filter(depth_range > 4) |>
    distinct(Year, Week)
  
  # Step 2: Filter original data for only those valid Week-Year combos
  df_filtered <- df_prepped |>
    semi_join(valid_weeks, by = c("Year", "Week")) |>
    group_by(Week, Year)
  
  # Step 3: Loop through each variable to summarize
  summary_list <- list()
  
  for (var in variables) {
    var_sym <- sym(var)
    
    summary_df <- df_filtered |>
      reframe(
        !!paste0("depth_", var, "_max") := Depth_m[which.max(!!var_sym)],
        !!paste0("depth_", var, "_min") := Depth_m[which.min(!!var_sym)],
        !!paste0(var, "_max_val") := if (all(is.na(!!var_sym))) NA_real_ else max(!!var_sym, na.rm = TRUE),
        !!paste0(var, "_min_val") := if (all(is.na(!!var_sym))) NA_real_ else min(!!var_sym, na.rm = TRUE),
        !!paste0(var, "_range") := (!!sym(paste0(var, "_max_val"))) - (!!sym(paste0(var, "_min_val")))
      )
    
    summary_list[[var]] <- summary_df
  }
  
  # Combine all summaries by Week and Year
  combined_summary <- reduce(summary_list, left_join, by = c("Week", "Year"))
  
  # Step 4: Add one Date per Week-Year combo (e.g., earliest)
  week_dates <- df_filtered |>
    group_by(Week, Year) |>
    summarise(Date = min(Date), .groups = "drop")
  
  combined_summary |>
    left_join(week_dates, by = c("Week", "Year")) |>
    select(Date, everything())  # Optional: move Date to the front
}

