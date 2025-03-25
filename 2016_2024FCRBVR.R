#exploring phytoplankton from 2016 to 2024 in FCR and BVR

#in FCR
FCRcheck <- current_df|>
  mutate(Date  = as_date(DateTime)) |> 
  filter(year(DateTime) == 2016, Reservoir == "FCR")
#max 2016 is 90 ugL on 2016-08-15

FCRcheck <- current_df|>
  mutate(Date  = as_date(DateTime)) |> 
  filter(year(DateTime) == 2024, Reservoir == "FCR")
#max 2024 is 363 ugL on 2024-09-02 

BVRcheck <- current_df|>
  mutate(Date  = as_date(DateTime)) |> 
  filter(year(DateTime) == 2016, Reservoir == "BVR")

2016-06-30

BVRcheck <- current_df|>
  mutate(Date  = as_date(DateTime)) |> 
  filter(year(DateTime) == 2024, Reservoir == "BVR")

2024-08-05 339 ugL