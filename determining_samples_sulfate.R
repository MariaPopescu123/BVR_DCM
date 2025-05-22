#Determining Sulfur sampling for FCR and BVR 

#for FCR

#first let's look at where the phytoplankton blooms are occurring
# for 2024


FCR2024 <- flora_heatmap(fp_data = phytos, year = 2024, site = 50, z = "TotalConc_ugL", unitz = "ug/L", max_legend_value = max(phytos$TotalConc_ugL))

#let's see max concentration for the year

FCRmaxphytos <- phytos|>
  filter(Year == 2024)

#answer is 3.5ish

checking <- chemistry|>
  filter(Reservoir == "FCR")
