library(haven)
popksq <- read_sas("Z:/Sante/Cancers/TabK82/TABLES/d01_20191202_popksq1q11.sas7bdat")

library(tidyverse)
dat <- popksq %>% filter(ksein == 1 | ksein1 == 1 | ksein2 == 1 | ksein3 == 1 | ksein4 == 1 )

write_sas(dat, path = "subset_ksein.sas7bdat")
saveRDS(dat, file = "subset_ksein.rds")


# EPIC IDs
popepic <- read_sas("//31.10.11.235/Donnees-BRUT/E3N/DEPIC/Population/pop_epic.sas7bdat")
saveRDS(popepic, file = "pop_epic.rds")

idepic <- read_sas("//31.10.11.235/Donnees-BRUT/E3N/Q3/Q3Orange/TABLES/idepic7.sas7bdat")
saveRDS(idepic, file = "idepic7.rds")
