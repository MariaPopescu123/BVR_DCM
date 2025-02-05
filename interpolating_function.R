#function for interpolating variables in a data frame
#work in progress

metalsdf_filtered <- metalsdf |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  filter(!if_any(starts_with("Flag"), ~. == 68))

variables <- c("SFe_mgL")

#test
interpolate_variable(metalsdf_filtered, variables, expanded_dates)

interpolate_variable <- function(data, variable_list, expanded_dates) {
  
  if (missing(expanded_dates) || !exists("expanded_dates")) {
    stop("Error: `expanded_dates` is required but not provided.")
  }
  
  data <- data |>
    filter(Reservoir == "BVR", Site == 50) |>
    mutate(Date = as_date(DateTime),
           Week = week(Date),
           Year = year(Date),
           DOY = yday(Date))
  
  interpolated_results <- list()  # Store results for each variable
  
  for (var in variable_list) {
    var_daily_summarise <- data |>
      select(Depth_m, Year, Week, DOY, Date, all_of(var)) |>
      group_by(Year, DOY, Week, Date, Depth_m) |>
      summarise(!!var := mean(.data[[var]], na.rm = TRUE), .groups = "drop") |>
      ungroup()
    
    var_depth_rounded <- var_daily_summarise |>
      group_by(Date) |>
      mutate(Depth_m = round(Depth_m, digits = 1)) |>
      ungroup() |>
      group_by(Date, Depth_m) |>
      summarise(!!var := mean(.data[[var]], na.rm = TRUE), .groups = "drop") |>
      ungroup() |>
      mutate(Year = year(Date),
             Week = week(Date))
    
    var_weekly <- var_depth_rounded |>
      group_by(Year, Week, Depth_m) |>
      summarise(!!var := mean(.data[[var]], na.rm = TRUE), .groups = "drop")
    
    var_interpolated <- expanded_dates %>%
      left_join(var_weekly, by = c("Depth_m", "Week", "Year")) %>%
      
      group_by(Year, Week) %>%
      mutate(
        first_valid_depth = ifelse(all(is.na(.data[[var]])), NA_real_, min(Depth_m[!is.na(.data[[var]])], na.rm = TRUE)),
        last_valid_depth = ifelse(all(is.na(.data[[var]])), NA_real_, max(Depth_m[!is.na(.data[[var]])], na.rm = TRUE)),
        Value_interp_depth = ifelse(
          Depth_m >= first_valid_depth & Depth_m <= last_valid_depth,
          zoo::na.approx(.data[[var]], x = Depth_m, na.rm = FALSE),
          NA_real_
        )
      ) %>%
      ungroup() %>%
      
      group_by(Year, Depth_m) %>%
      mutate(
        first_valid_Week = ifelse(all(is.na(.data[[var]])), NA_real_, min(Week[!is.na(.data[[var]])], na.rm = TRUE)),
        last_valid_Week = ifelse(all(is.na(.data[[var]])), NA_real_, max(Week[!is.na(.data[[var]])], na.rm = TRUE)),
        Value_interp_Week = ifelse(
          Week >= first_valid_Week & Week <= last_valid_Week,
          zoo::na.approx(.data[[var]], x = Week, na.rm = FALSE),
          NA_real_
        )
      ) %>%
      ungroup() %>%
      
      mutate(interp_var = coalesce(Value_interp_depth, Value_interp_Week)) %>%
      
      select(-matches(var), -first_valid_depth, -last_valid_depth, -Value_interp_depth,
             -first_valid_Week, -last_valid_Week, -Value_interp_Week) %>%
      rename(!!var := interp_var)
    
    interpolated_results[[var]] <- var_interpolated
  }
  
  return(bind_rows(interpolated_results))
}

  
looking<- interpolate_variable(metalsdf_filtered, variables, expanded_dates)

looking <- looking|>
  filter(!is.na(SFe_mgL))
