# Preparation of data for BC study (after exploration of data)
# 2 datasets for continuous and categorical analysis
library(haven)
library(tidyverse)
library(readxl)
# Get CODBMB/IDENT correspondance
#idents <- read_xls("E3N_cancer du sein_21072014.xls") %>% select(1:2)

# Get supplementary reproductive variables
#rv <- read_sas("d_grossesse_20190107_corrections.sas7bdat")
#men <- read_sas("d01_menopauseq1.sas7bdat")
# Bind together and add IDENT
#rv <- men %>% bind_cols(rv) %>% left_join(idents, by = "IDENT")

# Add ident and duration of OC use
idents <- read_xls("E3N_cancer du sein_21072014.xls") %>% select(1:2)
oc <- read_sas("oc_joe.sas7bdat") %>% select(IDENT, totduroc, ppille)

# First get metadata. Set control variables to cases, calculate age at diagnosis, add reproductive variables
meta <- read_csv("metadata.csv", na = "9999") %>%
  fill(c(DIAGSAMPLING, ER), .direction = "downup") %>% 
  mutate(DURTHSBMBCat = ifelse(DURTHSBMB > 0, 1, 0), AGEdiag = AGE + DIAGSAMPLING,
         CODBMB = as.character(CODBMB)) %>%
  #mutate(CODBMB, as.character) %>%
  left_join(idents, by = "CODBMB") %>%
  left_join(oc, by = "IDENT") %>%
  #left_join(rv, by = "IDENT") %>%
  mutate(durOC = if_else(totduroc > 120, "high", "low")) %>%
  mutate_at(vars(SMK, DIABETE, durOC), as.factor)
  
meta$DIAGSAMPLINGQ4 <- cut_number(meta$DIAGSAMPLING, 4, labels = F)

# For removal of case-control pairs not matched by menopausal status (univariate)
unmatch.pairs <- meta %>% group_by(MATCH) %>% summarise(sum.men = sum(MENOPAUSE)) %>% 
  filter(sum.men == 1) %>% select(MATCH) %>% pull() %>% as.numeric

log.vec <- !(meta$MATCH %in% unmatch.pairs)
meta <- meta[log.vec, ]

#unmatch.pairs1 <- tapply(meta$MENOPAUSE, meta$MATCH, sum)
#log.vec1 <- !(meta$MATCH %in% unmatch.pairs1)
#meta1 <- meta[which(unmatch.pairs1) == T, ]

# For subsetting
pre <- meta$MENOPAUSE == 0
post <- meta$MENOPAUSE == 1 
#agehi <- meta$Age1 == 1
#agelo <- meta$Age1 == 0
fast <- meta$FASTING == 1
fast1 <- meta$FASTING == 1 & meta$MENOPAUSE == 1
fast0 <- meta$FASTING == 1 & meta$MENOPAUSE == 0

# Oestrogen receptor positive
pos <- meta$ER == 1
neg <- meta$ER == 0

# follow-up
fuQ1 <- meta$DIAGSAMPLINGQ4 == 1
fuQ2 <- meta$DIAGSAMPLINGQ4 == 2
fuQ3 <- meta$DIAGSAMPLINGQ4 == 3
fuQ4 <- meta$DIAGSAMPLINGQ4 == 4

# For sensitivity analyses: less than 55 at diagnosis
preS2 <- meta$MENOPAUSE == 0 & meta$AGEdiag < 55
pre0 <- meta$MENOPAUSE == 0 & meta$DIAGSAMPLING > 2

# Menopausal status at diagnosis: compare age at menopause and age at diagnosis
#meta$menodiag <- ifelse(meta$AGEdiag < as.numeric(meta$MENOAGE), 0, 1)
#ages <- select(meta, MENOPAUSE, menodiag) %>% filter(MENOPAUSE == 0)
# Only about 8 pre-menopausal at blood collection were not at diagnosis

# Compound data for forest plots
cmpd.meta <- read.csv("NMR_cmpd_metadata_new.csv")
cmpds.ordered <- cmpd.meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description)-1))
rowvec <- cmpds.ordered$row

# Unscaled data
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")

#ints.all <- read.delim("1507_XMetabolite_std_cpmg_E3N.txt")
#samples <- ints.all$CODBMB %in% meta$CODBMB
#ints <- ints.all[samples, -1]

# Remove problem samples
#filt <- !(meta$RACK %in% c(10, 29, 33, 34))
ints0 <- ints0[log.vec, ]

# Remove top and bottom 1% of each compound and replace with NA
outliers <- function(x) x > quantile(x, probs = 0.99) | x < quantile(x, probs = 0.01)
logicalmat <- apply(ints0, 2, outliers)
ints0[logicalmat] <- NA

# Scale to unit variance
ints <- scale(ints0)

# Revisions for reviewers 13 December 2020
table(fast = meta$FASTING, meno = meta$MENOPAUSE)
# 65.0% pre-menopasual non-fasting, 63.4% post-menopausal non-fasting

# Other preparation steps (not used)

# Replace negative values with half the minimum positive value (not used)
#rm.neg.values <- function(x) ifelse(x < 0, min(x[x > 0])/2, x)
#ints0 <- apply(ints0, 2, rm.neg.values)
# Check
#which(apply(as.matrix(ints0), 2, min) < 0)

# For removal of problem racks for Ethanol
#meta1 <- meta %>% filter(!RACK %in% c(10, 29, 33, 34))

# Get a subset of the most discriminating compounds
# Histidine, NAC, glycerol, ornithine, ethanol, pyruvate, albumin, glutamate, glutamine, 3 fatty acids
#ss <- ints[, c(5,14,15,16,17,19,20,27,28,33,41,42)]






