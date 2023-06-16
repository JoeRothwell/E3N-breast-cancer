# Subject selection for TNBC. See email from Amandine 17 May 2021
# 8238 instances of BC downloaded from server
library(haven)
library(tidyverse)
library(readxl)

popksq <- read_sas("Z:/Sante/Cancers/TabK82/TABLES/d01_20191202_popksq1q11.sas7bdat")

library(tidyverse)
dat <- popksq %>% filter(ksein == 1 | ksein1 == 1 | ksein2 == 1 | ksein3 == 1 | ksein4 == 1 )

#write_sas(dat, path = "subset_ksein.sas7bdat")
#saveRDS(dat, file = "subset_ksein.rds")

#bc1 <- read_sas("subset_ksein.sas7bdat")
bc <- readRDS("subset_ksein.rds")

# Get confirmed cases----
bcc <- 
  bc %>% filter(investstop1 == 1 | investstop2 == 1 | investstop3 == 1 | investstop4 == 1)
bcc <- bc

# Get triple negative cases
bcer <- bcc %>% filter(ro1 == "Négatif" | ro2 == "Négatif" | ro3 == "Négatif" | ro4 == "Négatif")
bcpr <- bcer %>% filter(rp1 == "Négatif" | rp2 == "Négatif" | rp3 == "Négatif" | rp4 == "Négatif")
bctn <- bcpr %>% filter(cerb2_1 == "Pas de surexpression" | cerb2_2 == "Pas de surexpression" | 
                          cerb2_3 == "Pas de surexpression" | cerb2_4 == "Pas de surexpression")

table(bctn$pre_inc1)
table(bctn$pre_inc2)
# No prevalent cancers in this subset
# Variables pre_inc3 and 4 are empty

table(bctn$conco) # only 12 of 256 are concomitant
min(bctn$ddiag) # First diagnosis is 18.12.1990
plot(bctn$ddiag) # Only 4 before 1995

# Warning: if repeating the analysis, add code to exclude the few concomitant
# tumours that are mistakenly classified as TN (almalgamation of status).

# Get ident
unique(bctn$ident) %>% length

# Select important variables only (tumour IDs also added)
dat <- bctn %>% select(ident, conco, pre_inc1, pre_inc2, ddiag, ro1, ro2, rp1, rp2, 
                       cerb2_1, cerb2_2, id_tum1, id_tum2)

#dat <- bctn %>% select(IDENT = ident, conco, pre_inc1, pre_inc2, ddiag, 
#                       ro1, ro2, rp1, rp2, cerb2_1, cerb2_2)

#write.csv(dat, "Ident_triple_negatif.csv")

# Biobank data from Rahime and join to original data
bb <- read_xlsx("Ident_Triple_Negatif.xlsx")
tn <- dat %>% inner_join(bb, by = "ident") %>% filter(disponibilite == T)
write.csv(tn, "TNBC_biobank_blood.csv")

plot(tn$ddiag)

# Compare Ident to case control: run prep_data down to meta
# Can omit and go straight to tn.blood below
idents <- read_xls("E3N_cancer du sein_21072014.xls") %>% select(1:2)
meta <- read_csv("metadata.csv", na = "9999") %>% mutate(CODBMB = as.character(CODBMB)) %>%
  left_join(idents, by = "CODBMB")
meta.tn <- meta %>% inner_join(dat, by = "IDENT") #28 TN cases

# Get samples to look for from Laure's file
tn.blood <- read_xlsx("E3N TNBC blood.xlsx") %>% filter(INCLUSION %in% c(3,6)) %>%
  select(ident, INCLUSION)
tn.blood$ident <- as.character(tn.blood$ident)
meta.tn <- tn.blood %>% inner_join(dat, by = "ident") 
write.csv(meta.tn, "recherche_tumeurs_TN.csv")

# Export to send to Rafika to look for tumours


# From previously extracted data, for comparison
sum(meta$CERB2 == 0 & meta$ER == 0 & meta$PR == 0, na.rm = T) # 22 TN cases


### Get positive samples (ER and PR positive)
bcer <- bcc %>% filter(ro1 == "Positif" | ro2 == "Positif" | ro3 == "Positif" | ro4 == "Positif")
bcdp <- bcer %>% filter(rp1 == "Positif" | rp2 == "Positif" | rp3 == "Positif" | rp4 == "Positif") %>%
  # Filter by date according to instructions from Laure
  filter(ddiag > "2005-07-31")

dat <- bcdp %>% select(ident, conco, pre_inc1, pre_inc2, ddiag, 
                       ro1, ro2, rp1, rp2, cerb2_1, cerb2_2)
write.csv(dat, "Ident_erpos_prpos.csv")


# EPIC IDs----
popepic <- read_sas("//31.10.11.235/Donnees-BRUT/E3N/DEPIC/Population/pop_epic.sas7bdat")
saveRDS(popepic, file = "pop_epic.rds")

idepic <- read_sas("//31.10.11.235/Donnees-BRUT/E3N/Q3/Q3Orange/TABLES/idepic7.sas7bdat")
saveRDS(idepic, file = "idepic7.rds")
