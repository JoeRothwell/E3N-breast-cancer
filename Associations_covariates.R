# Association BC risk factors - metabolites
source("BC_prep_data.R")

meta1 <- meta %>% select(FASTING, SMK, BMI, AGE, BP, TTAILLE, Trait_Horm, MENOPAUSE, DIABETE, CO,
                         DURTHSBMB, CENTTIME, STOCKTIME)

# Risk factors to be tested
#FASTING, SMK, BMI, AGE, BP, TTAILLE, Trait_Horm, MENOPAUSE, DIABETE, CO
library(broom)
lm.fasting <- function(x) lm(FASTING ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)
lm.smk     <- function(x) lm(SMK ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)
lm.bmi     <- function(x) lm(BMI ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)
lm.age     <- function(x) lm(AGE ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)
lm.bp      <- function(x) lm(BP ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)
lm.ttaille <- function(x) lm(TTAILLE ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)
lm.horm    <- function(x) lm(Trait_Horm ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)
lm.meno    <- function(x) lm(MENOPAUSE ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)
lm.diabet  <- function(x) lm(DIABETE ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)
lm.co      <- function(x) lm(CO ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta)

fits1 <- apply(ints, 2, lm.fasting)
fits2 <- apply(ints, 2, lm.smk)
fits3 <- apply(ints, 2, lm.bmi)
fits4 <- apply(ints, 2, lm.age)
fits5 <- apply(ints, 2, lm.bp)
fits6 <- apply(ints, 2, lm.ttaille)
fits7 <- apply(ints, 2, lm.horm)
fits8 <- apply(ints, 2, lm.meno)
fits9 <- apply(ints, 2, lm.diabet)
fits10 <- apply(ints, 2, lm.co)

ll <- list(fits1, fits2, fits3, fits4, fits5, fits6, fits7, fits8, fits9, fits10)
map_df(ll, tidy)

fasting <- map_df(fits1, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))
smk <- map_df(fits2, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))
bmi <- map_df(fits3, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))
age <- map_df(fits4, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))
bp <- map_df(fits5, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))
ttaille <- map_df(fits6, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))
horm <- map_df(fits7, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))
meno <- map_df(fits8, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))
diabet <- map_df(fits9, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))
co <- map_df(fits10, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(colnames(ints))

all <- bind_rows(fasting, bmi, ttaille, smk, age, bp, horm, meno, diabet, co)
library(RColorBrewer)
mycols <- rep(economist_pal(fill=TRUE)(3), each = 43, length.out=nrow(all))
mycols2 <- rep(brewer.pal(8, "Set2") , each = 43, length.out=nrow(all))
plot(-log10(all$p.value), col = mycols2, pch = 19, cex = 0.7)
