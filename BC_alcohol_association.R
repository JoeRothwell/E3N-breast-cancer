# Associations between breast cancer and alcohol by menopausal status

library(tidyverse)
library(readxl)

ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
ints <- scale(ints0)

#cmpd.meta <- read.csv("NMR_cmpd_metadata_new.csv")
#eth <- ints[, 14]

# Get data subsets for post and pre-menopausal women
meta <- read_csv("metadata.csv", na = "9999") %>% mutate_at(vars(SMK, DIABETE, RACK, MATCH, SAMPYEAR), as.factor)
meta.pre <- meta %>% filter(MENOPAUSE == 0)
meta.pos <- meta %>% filter(MENOPAUSE == 1)
metaF <- meta %>% filter(FASTING == 1)
meta.preF <- meta %>% filter(MENOPAUSE == 0 & FASTING == 1)
meta.posF <- meta %>% filter(MENOPAUSE == 1 & FASTING == 1)

# Models for alcohol intake (remove DURTHSDIAG from pre and post)
library(survival)
f0 <- clogit(CT ~ I(ALCOHOL/10) + BMI + SMK + DIABETE + RTH + DURTHSDIAG + strata(MATCH), data = meta)
f1 <- clogit(CT ~ I(ALCOHOL/10) + BMI + SMK + DIABETE + RTH + strata(MATCH), data = meta.pre)
f2 <- clogit(CT ~ I(ALCOHOL/10) + BMI + SMK + DIABETE + RTH + strata(MATCH), data = meta.pos)

# Fasting models
f3 <- clogit(CT ~ I(ALCOHOL/10) + BMI + SMK + DIABETE + RTH + DURTHSDIAG + strata(MATCH), data = metaF)
f4 <- clogit(CT ~ I(ALCOHOL/10) + BMI + SMK + DIABETE + RTH + strata(MATCH), data = meta.preF)
f5 <- clogit(CT ~ I(ALCOHOL/10) + BMI + SMK + DIABETE + RTH + strata(MATCH), data = meta.posF)

#Fasting models

library(broom)
mods <- map_df(list(f0, f1, f2, f3, f4, f5), ~tidy(., exponentiate = T, conf.int = T)) %>% 
        filter(str_detect(term, "ALCOHOL/10")) %>% mutate_if(is.numeric, ~round(., 2)) %>% 
        unite(OR.CI, estimate, conf.low, conf.high, sep = "-")

# Result: Alcohol is borderline associated with risk overall and post-menopausal.



# Models adjusted for ethanol intake (old)
#cmpd.meta <- read.csv("NMR_cmpd_metadata_new.csv")
#eth <- ints[, 14]
#meta <- cbind(meta, eth)

library(survival)
base <- CT ~ eth + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH)
fit1 <- clogit(update(base, . ~ .), data = meta)

# categorical

quartiles <- ints0 %>% mutate_all(funs(cut_number(., n = 4, labels = 1:4))) 
eth_q <- cut_number(pull(ints0[, 14]), n = 4, labels = 1:4)
Q1Q4 <- eth_q %in% c(1, 4)
meta <- cbind(meta, eth_q)
meta2 <- meta[Q1Q4, ]
# Need to drop levels 2 and 3 or else model does not converge
meta2$eth_q <- droplevels(meta2$eth_q)

fit3 <- clogit(CT ~ eth_q + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSDIAG + CENTTIME + STOCKTIME +
         strata(MATCH), data = meta2)

library(broom)
t <- map_df(list(fit1, fit2), tidy) %>% filter(str_detect(term, "eth"))

par(mar=c(5,4,2,2))
library(metafor)
forest(t$estimate, ci.lb = t$conf.low, ci.ub = t$conf.high, refline = 1,
       xlab = "OR per SD increase in concentration", 
       transf = exp, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"))

t1 <- map_df(fit1, tidy) %>% filter(str_detect(term, "eth"))

forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1,
       xlab = "OR per SD increase in concentration", 
       transf = exp, efac = 0.5,
       pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"))