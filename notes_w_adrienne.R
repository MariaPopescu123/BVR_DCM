
metalsdf_filtered <- metalsdf |>
  filter(Reservoir == "BVR", Site == 50)|>
  mutate(Date = as_date(DateTime))|>
  filter(!if_any(starts_with("Flag"), ~. == 68))

variables <- c("SFe_mgL", "TFe_mgL", "TMn_mgL", "SMn_mgL")

#test
test <- interpolate_variable(metalsdf_filtered, variables, expanded_dates)  

looking <- final_result|>
  filter(!is.na(SFe_mgL))

#first clean up dataframes

#pmap purr package
#the argument would be the data frame and in the data
#frame you have a column for each argument in your function that you made

#all_instruments_interpolated <- df_witharguments|>
#purrr::pmap(function)
#then afterwards join all(all_instruments_interpolated)

#look at fcr catwalk in edi already published

