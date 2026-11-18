
################################################
# PART 1 - Download entire database ---------------------------------------
################################################

# Downloaded on 2nd November 2026
# PREDICTS citation (Hudson et al., 2017)
# Hudson, L.N., Newbold, T., Contu, S., Hill, S.L., Lysenko, I., De Palma, A., Phillips, H.R., Alhusseini, T.I., Bedford, F.E., Bennett, D.J. and Booth, H., 2017. The database of the PREDICTS (projecting responses of ecological diversity in changing terrestrial systems) project. Ecology and evolution, 7(1), pp.145-188.

# 2022 release citation (Contu et al., 2022)
# Sara Contu; Adriana De Palma; Rachel Bates; Jessica Borer et al. (2022). Release of data added to the PREDICTS database (November 2022) [Data set]. Natural History Museum. https://doi.org/10.5519/jg7i52dg

install.packages("predictsr")
citation(package = 'predictsr')

library(dplyr)
library(ggplot2)
library(here)
library(predictsr)

# Only run this section once:

# dir.create(here("Data"), showWarnings = FALSE)

# Database 
PREDICTS <- GetPredictsData(extract = c(2016, 2022))

# Description of data extract columns
col_descriptions <- GetColumnDescriptions()

# Site-level summaries
site_summaries <- GetSitelevelSummaries(extract = c(2016, 2022))

# Save locally so you don't need to re-download
write.csv(PREDICTS, here("Data", "PREDICTS_database_4Nov26.csv"),   row.names = FALSE)
write.csv(col_descriptions, here("Data", "PREDICTS_descriptions.csv"), row.names = FALSE)
write.csv(site_summaries, here("Data", "PREDICTS_site_summaries_4Nov26.csv"), row.names = FALSE)


library(rnaturalearth)
library(tidyterra)
library(sf)

studies <- predicts %>% group_by(SS) %>% 
  summarise(n_spp = n_distinct(Taxon_name_entered),
            Longitude = mean(Longitude),
            Latitude = mean(Latitude))

ne_countries(scale = 50,returnclass = 'sv') -> world

p <- ggplot() + 
  geom_spatvector(data = world,fill="#BF2C2333",col="#00000000") + 
  geom_point(data = studies,
             mapping = aes(x = Longitude,y = Latitude,size=n_spp),
             col = "#1065AB66") +
  scale_x_continuous(name = "") + 
  scale_y_continuous(name = "") + 
  theme_classic()

p

###########################################
# Data pre-processing -----------------------------------------------------
###########################################

# First, for the 0.8% of records that are of an effort-sensitive metric, and where sampling effort differs among sampled sites within the study, we need to correct sampling effort. This is done by rescaling sampling effort within each study, such that the most-sampled site gets a relative effort value of 1, and the sampling effort for all other sites is scaled linearly relative to this maximum value. The CorrectSamplingEffort function in the predictsFunctions package carries out these operations:



