#function for interpolating variables in a data frame
#make sure to clean up data before hand. remove flags 


interpolate_variable <- function(data, variable_list, expanded_dates) {
  
  # data <- metalsdf_filtered
  # variable_list <- c("SFe_mgL", "TFe_mgL")
  # expanded_dates <- expanded_dates
  
  if (missing(expanded_dates) || !exists("expanded_dates")) {
    stop("Error: `expanded_dates` is required but not provided.")#do i need this
  }
  
  data <- data |> #changing this temporarily 
    filter(Reservoir == "BVR", Site == 50) |>
    mutate(Date = as_date(DateTime),
           Week = week(Date),
           Year = year(Date),
           DOY = yday(Date))
  
  interpolated_results <- list()  # Store results for each variable
  

  for (var in variable_list) {
    var_daily_summarise <- data |> #sumarising the data per day per depth for replicates
      select(Depth_m, Year, Week, DOY, Date, all_of(var)) |>
      group_by(Year, DOY, Week, Date, Depth_m) |>
      summarise(!!sym(var) := mean(.data[[var]], na.rm = TRUE), .groups = "drop") |>
      ungroup()
    
    var_depth_rounded <- var_daily_summarise |> #rounding to the 0.1 meter and getting the mean
      group_by(Date) |>
      mutate(Depth_m = round(Depth_m, digits = 1)) |>
      ungroup() |>
      group_by(Date, Depth_m) |>
      summarise(!!sym(var) := mean(.data[[var]], na.rm = TRUE), .groups = "drop") |>
      ungroup() |>
      mutate(Year = year(Date),
             Week = week(Date))
    
    var_weekly <- var_depth_rounded |>
      group_by(Year, Week, Depth_m) |>
      summarise(!!sym(var) := mean(.data[[var]], na.rm = TRUE), .groups = "drop")
    
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
      
        #make sure 
      group_by(Year, Depth_m) %>%
      mutate(
        first_valid_Week = ifelse(all(is.na(Value_interp_depth)), NA_real_, min(Week[!is.na(Value_interp_depth)], na.rm = TRUE)), #setting the range where to interpolate
        last_valid_Week = ifelse(all(is.na(Value_interp_depth)), NA_real_, max(Week[!is.na(Value_interp_depth)], na.rm = TRUE)),
        Value_interp_Week = ifelse(
          Week >= first_valid_Week & Week <= last_valid_Week,
          zoo::na.approx(Value_interp_depth, x = Week, na.rm = FALSE),
          NA_real_
        )
      ) %>%
      ungroup() %>%
      
      mutate(interp_var = coalesce(Value_interp_depth, Value_interp_Week)) %>%
      
      select(-matches(var), -first_valid_depth, -last_valid_depth, -Value_interp_depth,
             -first_valid_Week, -last_valid_Week, -Value_interp_Week) %>%
      rename(!!sym(var) := interp_var)
    
    interpolated_results[[var]] <- var_interpolated  # Store in list dynamically
  }
  
  
  final_result <- plyr::join_all(interpolated_results, by = c("Week", "DOY", "Depth_m", "Year", "Date"))  # Combine all results
  return(final_result)
  
}

