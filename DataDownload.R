# Downloading data

beep <- function(){
  system("rundll32 user32.dll,MessageBeep") #just putting this here so I can be notified when something that takes a while to run is done
}  

pacman::p_load(tidyverse, lubridate, akima, reshape2, 
               gridExtra, grid, colorRamps, RColorBrewer, rLakeAnalyzer,
               reader, cowplot, dplyr, tidyr, ggplot2, zoo, purrr, beepr, forecast)

#### Loading Data  ####

#ctd data
CTD <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/200/14/0432a298a90b2b662f26c46071f66b8a")
write.csv(CTD, "./CTD.csv", row.names = FALSE )

#flora data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.272.8
current_df <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/272/8/0359840d24028e6522f8998bd41b544e")
write.csv(current_df, "./current_df.csv", row.names = FALSE )


# metals data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.455.8
metalsdf <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/455/8/9c8c61b003923f4f03ebfe55cea8bbfd")
#removed flags for 68 as per Cece's advice
write.csv(metalsdf, "./metalsdf.csv", row.names = FALSE )

#ghgs data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.551.8
ghgs <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/551/8/454c11035c491710243cae0423efbe7b")
#not sure whether or not to remove those with flags 3 and 4
#3 = The difference between the reps are above the limit of quantification and >30% and <50% different from each other. Both replicates were retained but flagged
#4 = The difference between the reps are above the limit of quantification and >50% different from each other. Both replicates were retained but flagged
write.csv(ghgs, "./ghgs.csv", row.names = FALSE )

#secchi data https://portal.edirepository.org/nis/mapbrowse?packageid=edi.198.11
secchiframe <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/198/11/81f396b3e910d3359907b7264e689052")
write.csv(secchiframe, "./secchiframe.csv", row.names = FALSE )

#ysi https://portal.edirepository.org/nis/mapbrowse?packageid=edi.198.11
ysi_profiles <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/198/11/6e5a0344231de7fcebbe6dc2bed0a1c3")
write.csv(ysi_profiles, "./ysi_profiles.csv", row.names = FALSE)

#data from here https://portal.edirepository.org/nis/mapbrowse?packageid=edi.199.12
chemistry <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/199/12/a33a5283120c56e90ea414e76d5b7ddb")
write.csv(chemistry, "./chemistry.csv", row.names = FALSE)

#meteorological data from FCR https://portal.edirepository.org/nis/mapbrowse?packageid=edi.389.8
options(timeout = 300)
metdata <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/389/8/d4c74bbb3b86ea293e5c52136347fbb0")
write.csv(metdata, "./metdata.csv", row.names = FALSE)

#bathymetry data for BVR https://portal.edirepository.org/nis/metadataviewer?packageid=edi.1254.1
bath <- read.csv("https://portal.edirepository.org/nis/dataviewer?packageid=edi.1254.1&entityid=f7fa2a06e1229ee75ea39eb586577184")
write.csv(bath, "./bath.csv", row.names = FALSE)

#waterlevel data
wtrlvl <- read.csv("https://pasta.lternet.edu/package/data/eml/edi/725/4/43476abff348c81ef37f5803986ee6e1") 
write.csv(wtrlvl, "./wtrlvl.csv", row.names = FALSE)
