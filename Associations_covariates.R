# Association BC risk factors - metabolites
library(tidyverse)
meta <- read_csv("metadata.csv", na = "9999")
ctrl <- meta$CT == 0 & meta$MENOPAUSE == 1

# Get control subjects and remove outliers from metabolite data
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")[ctrl, ]
outliers <- function(x) x > quantile(x, probs = 0.99) | x < quantile(x, probs = 0.01)
logicalmat <- apply(ints0, 2, outliers)
ints0[logicalmat] <- NA

# Scale to unit variance
ints <- scale(ints0)
meta1 <- meta %>% filter(CT == 0 & MENOPAUSE == 1) %>% select(FASTING, SMK, ALCOHOL, Life_Alcohol_Pattern_1, 
      BMI, AGE, BP, TTAILLE, Trait_Horm, MENOPAUSE, DIABETE, CO, DURTHSBMB, CENTTIME, STOCKTIME,
      WEEKS)

# Risk factors to be tested
#FASTING, SMK, BMI, AGE, BP, TTAILLE, Trait_Horm, MENOPAUSE, DIABETE, CO
library(broom)
lm.fasting <- function(x) lm(FASTING ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.alc     <- function(x) lm(ALCOHOL ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.lifealc <- function(x) lm(Life_Alcohol_Pattern_1 ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.smk     <- function(x) lm(SMK ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.bmi     <- function(x) lm(BMI ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.age     <- function(x) lm(AGE ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.bp      <- function(x) lm(BP ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.ttaille <- function(x) lm(TTAILLE ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.horm    <- function(x) lm(Trait_Horm ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.diabet  <- function(x) lm(DIABETE ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.co      <- function(x) lm(CO ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)

fits1 <- apply(ints, 2, lm.fasting) %>% map_df(tidy) 
fits1a <- apply(ints, 2, lm.alc) %>% map_df(tidy) 
fits1b <- apply(ints, 2, lm.lifealc) %>% map_df(tidy) 
fits2 <- apply(ints, 2, lm.smk) %>% map_df(tidy) 
fits3 <- apply(ints, 2, lm.bmi) %>% map_df(tidy)
fits4 <- apply(ints, 2, lm.age) %>% map_df(tidy) 
fits5 <- apply(ints, 2, lm.bp) %>% map_df(tidy) 
fits6 <- apply(ints, 2, lm.ttaille) %>% map_df(tidy) 
fits7 <- apply(ints, 2, lm.horm) %>% map_df(tidy) 
fits9 <- apply(ints, 2, lm.diabet) %>% map_df(tidy) 
fits10 <- apply(ints, 2, lm.co) %>% map_df(tidy)

all <- bind_rows(fits1, fits1a, fits1b, fits2, fits3, fits4, fits5, fits6, fits7,
                 fits9, fits10) %>% filter(str_detect(term, "x")) 

library(RColorBrewer)
# Different colour palettes
set2 <- rep(brewer.pal(8, "Set2"), each = 43, length.out=nrow(all))
#accent <- rep(brewer.pal(8, "Accent") , each = 43, length.out=nrow(all))
#rain <- rep(rainbow(12), each = 43, length.out = nrow(all))


plot(-log10(all$p.value), col = set2, pch = 19, cex = 0.6)
#plot(-log10(all$p.value), col = accent, pch = 19, cex = 0.7)
#plot(-log10(all$p.value), col = rain, pch = 19, cex = 0.7)
