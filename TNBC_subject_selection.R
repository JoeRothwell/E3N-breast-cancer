# Subject selection for TNBC. See email from Amandine 17 May 2021
# 8238 instances of BC downloaded from server
library(haven)
library(tidyverse)
library(readxl)

#bc1 <- read_sas("subset_ksein.sas7bdat")
bc <- readRDS("subset_ksein.rds")

# Get confirmed cases
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

# Get ident
unique(bctn$ident) %>% length

# Select important variables only
dat <- bctn %>% select(IDENT = ident, conco, pre_inc1, pre_inc2, ddiag, 
                       ro1, ro2, rp1, rp2, cerb2_1, cerb2_2)

#write.csv(dat, "Ident_triple_negatif.csv")

# Compare Ident to case control: run prep_data down to meta
idents <- read_xls("E3N_cancer du sein_21072014.xls") %>% select(1:2)
meta <- read_csv("metadata.csv", na = "9999") %>% mutate(CODBMB = as.character(CODBMB)) %>%
  left_join(idents, by = "CODBMB")
meta.tn <- meta %>% inner_join(dat, by = "IDENT") #28 TN cases

# From previously extracted data, for comparison
sum(meta$CERB2 == 0 & meta$ER == 0 & meta$PR == 0, na.rm = T) # 22 TN cases
