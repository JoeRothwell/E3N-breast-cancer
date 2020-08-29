library(haven)
library(tidyverse)
library(readxl)
library(broom)

# Datasets-----------------------------------------------------------------------

# Food data
alim <- read_sas("frjour.sas7bdat")

# Not used
#dat1 <- read_sas("nutfra2_all_nut.sas7bdat")
#dat2 <- read_sas("nutfra2_polyph.sas7bdat")
#dat3 <- read_sas("nutrfra2_index.sas7bdat")

# Correspondence ident-COBBMB
id <- read_xls("E3N_cancer du sein_21072014.xls") %>% mutate(ident = IDENT) %>% select(c("CODBMB", "ident"))

# Metadata
meta <- read_csv("metadata.csv", na = "9999")
meta$CODBMB <- as.character(meta$CODBMB)

# Extract controls from metadata and add food data
meta0 <- meta[meta$CT == 0, ]
meta.id <- meta0 %>% left_join(id, by = "CODBMB")
alim <- meta.id %>% left_join(alim, by = "ident")

# Get metabolomics data (unscaled)
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
ints.ctrl <- ints0[meta$CT == 0, ]

# Scale to unit variance
metabolo <- scale(ints.ctrl)

# We now have 791 observations of metadata and foods and 791 observations of metabolite data


# Correlations---------------------------------------------------------------------

# Wine
# Replace intake = NA with median value
alim$VIN[is.na(alim$VIN)] <- median(alim$VIN, na.rm = T)

# Simple correlation for wine
# Exclude zero intakes with alim$VIN > 0

simplecor <- function(x) cor.test(alim$VIN[alim$VIN > 0], x, method = "spearman")
corlist <- apply(metabolo[alim$VIN > 0, ], 2, simplecor)

# Convert to data frame and add compound names, order by correlation
library(broom)
cordat <- map_dfr(corlist, tidy) %>% bind_cols(compound = colnames(metabolo)) %>% arrange(-estimate)


# Partial correlation controlling for BMI and smoking

partialcor <- function(x) {
  
  # Linear model of food intake and confounders
  mod1 <- lm(VIN ~ BMI + SMK, data = alim[alim$VIN > 0, ])
  
  # Linear model of metabolites and confounders
  mod2 <- lm(x ~ BMI + SMK, data = alim[alim$VIN > 0, ])
  
  # Correlate the two sets of residuals              
  cor.test(residuals(mod1), residuals(mod2), type = "spearman")
  
}
  
# Convert to data frame and add compound names, order by correlation
pcorlist <- apply(metabolo[alim$VIN > 0, ], 2, partialcor)  
pcordat <- map_dfr(pcorlist, tidy) %>% bind_cols(compound = colnames(metabolo)) %>% arrange(-estimate)


