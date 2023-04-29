# Generate R markdown for baseline characteristics table
# To add something about triple negative?
# Source prep data to remove the 10 subjects not matched on menopausal status
source("BC_prep_data.R")
library(haven)
library(tidyverse)

# Get reproductive variables and CODBMB/IDENT correspondence and bind together (same number and order)
rv <- read_sas("d_grossesse_20190107_corrections.sas7bdat") %>% select(-IDENT)
men <- read_sas("d01_menopauseq1.sas7bdat")
rvmen <- men %>% bind_cols(rv)
ident <- read_xls("E3N_cancer du sein_21072014.xls") %>% select(1:2)
oc <- read_sas("oc_joe.sas7bdat") %>% select(IDENT, totduroc, ppille)

# 1582 subjects in table after removal of low quality spectra
meta0 <- read_csv("metadata.csv", na = "9999") %>% #filter(!(MATCH %in% unmatch_pairs)) %>%
  select(CODBMB, CT,             # case-control status
         AGE,                   
         BMICat1, 
         RTHCat1,               # Waist-hip ratio categorical
         MENOPAUSE,             # menopausal status at blood collection
         
         SMK, DIABETE, Life_Alcohol_Pattern_1, ALCOHOL, BP, 
         
         CO,                    # previous oral contraceptive use
         Trait_Horm,            # menopausal treatment therapy taken 24h before blood collection
         #DURTHSDIAG,            # Duration of use of therapy at date of diagnosis
         DURTHSBMB,             # Duration of use of therapy at blood collection
         #CENTTIMECat1,          # time before centrifugation (?) OMIT not of interest
         FASTING, 
         STOCKTIME,             # Storage time (years)
         BEHAVIOUR,             # Tumour behaviour
         SUBTYPE, 
         #CERB2,                # HER2 receptor OMIT mostly unknown
         ER,                    # Estrogen receptor positive
         PR,                    # Progesterone receptor positive
         SBR, GRADE, STADE, DIAGSAMPLING,
         DIAGSAMPLINGCat1) 

meta <- meta0 %>% 
  mutate(CODBMB = as.character(CODBMB), AGEdiag = AGE + DIAGSAMPLING) %>%
  left_join(ident, by = "CODBMB") %>%
  mutate_at(vars(-AGE, -AGEdiag, -STOCKTIME, -ALCOHOL, -DURTHSBMB, -CODBMB), as.factor) %>%
  left_join(rvmen, by = "IDENT")

# For reviewers comments: get duration of OC use
table()
meta1 <- meta %>% left_join(oc, by = "IDENT")

# Check how many pairs have duration of OC data
meta1$durOC <- ifelse(is.na(meta1$totduroc), 0, 1)


library(kableExtra)
library(qwraps2)
options(qwraps2_markup = "markdown")

summ <-
  list(#"Total subjects"                  = list("N"   =   ~ n()),
    
       "Age at blood collection (years)" = list("Mean"     =  ~ mean_sd(AGE, digits = 1)),
       
     "Menopausal status at blood collection" =
         
         list("Pre-menopausal"   =  ~ n_perc0(MENOPAUSE == 0, digits = 1),
              "Post-menopausal"  =  ~ n_perc0(MENOPAUSE == 1, digits = 1)),
     
     "Follow-up time to cancer diagnosis" =
       
       list("5 years or less"   =  ~ n_perc0(DIAGSAMPLINGCat1 == 1, digits = 1),
            "More than five years"  =  ~ n_perc0(DIAGSAMPLINGCat1 == 2, digits = 1)),
       
     "BMI" =    
       
        list("Underweight or normal"   =  ~ n_perc0(BMICat1 == 1, na_rm = T, digits = 1),
             "Overweight"              =  ~ n_perc0(BMICat1 == 2, na_rm = T, digits = 1),
             "Obese"                   =  ~ n_perc0(BMICat1 == 3, na_rm = T, digits = 1),
             "Unknown"                 =  ~ n_perc0(is.na(BMICat1), digits = 1)),
       
     "Waist to hip ratio" =           
       
        list("< 0.8" = ~ n_perc0(RTHCat1 == 0, na_rm = T, digits = 1),
             "> 0.8" = ~ n_perc0(RTHCat1 == 1, na_rm = T, digits = 1),
             "Unknown" = ~ n_perc0(is.na(RTHCat1), digits = 1)),
       
     "Smoking status" = 
       
        list("Yes" = ~ n_perc0(SMK == 1, digits = 1),
              "No"  = ~ n_perc0(SMK == 0, digits = 1)),
       
     "Diabetic status"  = 
       
        list("Yes" = ~ n_perc0(DIABETE == 1, digits = 1),
              "No"  = ~ n_perc0(DIABETE == 0, digits = 1)),
       
     "Lifetime alcohol drinking pattern" =
       
        list("Non-consumers (0 g/day)"      =  ~ n_perc0(Life_Alcohol_Pattern_1 == 0, na_rm = T, digits = 1),
              "Light consumers (1-10 g/day)" =  ~ n_perc0(Life_Alcohol_Pattern_1 == 1, na_rm = T, digits = 1),
              "Drinkers (>10 g/day)"         =  ~ n_perc0(Life_Alcohol_Pattern_1 == 2, na_rm = T, digits = 1),
              "Unknown"                      =  ~ n_perc0(is.na(Life_Alcohol_Pattern_1), digits = 1)),
     
     "Alcohol intake (g/day)" =
       
       list("Mean (SD)" = ~ mean_sd(ALCOHOL, digits = 1)),
       
     "Blood pressure" =
       
      list("Normal tension"  =  ~ n_perc0(BP == 0, na_rm = T, digits = 1),
              "Hypertension"    =  ~ n_perc0(BP == 1, na_rm = T, digits = 1),
              "Unknown"  =  ~ n_perc0(is.na(BP), digits = 1)),
     
     "Previous breastfeeding" =
       
       list("Yes" =   ~ n_perc0(allaitement == 1, na_rm = T, digits = 1),
            "No" =   ~ n_perc0(allaitement == 0, na_rm = T, digits = 0),
            "Unknown"  =  ~ n_perc0(is.na(allaitement), digits = 1)),
       
     "Previous oral contraceptive use" =
       
        list("Yes" = ~ n_perc0(CO == 1, digits = 1),
              "No" = ~ n_perc0(CO == 0, digits = 1)),
       
     "Menopause hormone therapy at blood collection" =
       
        list("Yes" =  ~ n_perc0(Trait_Horm == 1 & MENOPAUSE == 1, digits = 1),
             "No"  =  ~ n_perc0(Trait_Horm == 0 & MENOPAUSE == 1, digits = 1),
             "Unknown"  =  ~ n_perc0(is.na(Trait_Horm) & MENOPAUSE == 1, digits = 1)),
       
     "Duration of use of menopause hormonal treatment at baseline" =
       
        list("Mean (SD)" =  ~ mean_sd(DURTHSBMB[MENOPAUSE == 1], digits = 1)),
       
     "Fasting status" =
       
        list("Fasting"     = ~ n_perc0(FASTING == 1, digits = 1),
              "Non-fasting" = ~ n_perc0(FASTING == 0, digits = 1)),
       
     "Biobank storage time (years)" =
       
        list("Mean (SD)" = ~ mean_sd(STOCKTIME, digits = 1))

  )

# By case-control status
st <- meta %>% group_by(CT) %>% summary_table(summ)
print(st, cnames = c("Controls (N=791)", "Cases (N=791)"))

# By case-control and menopausal status (for reviewers)
st <- meta %>% group_by(CT, MENOPAUSE) %>% summary_table(summ)
print(st, cnames = c("Pre, controls (N=180)", "Pre, cases (N=179)", 
                     "Post, controls (N=611)", "Post, cases (N=612)"))
# Copy and paste output into table1_bc_ms.Rmd and render to word/pdf etc


# Automatic (can't make work properly)
#summ <- meta0 %>% qsummary(., numeric_summaries = list("Mean (SD)" = "~ mean_sd(%s)"),
#                                   n_perc_args = list(digits = 1, show_symbol = T))


# Table of tumour characteristics by menopausal status
summ1 <-
  list("Age at diagnosis" =
      list("Mean (SD)" =  ~ mean_sd(AGEdiag, na_rm = T, digits = 1)),
    
    "Time between sampling and diagnosis" =
      
      list("5 years or less"   =  ~ n_perc0(DIAGSAMPLINGCat1 == 1, na_rm = T, digits = 1),
           "More than 5 years" =  ~ n_perc0(DIAGSAMPLINGCat1 == 2, na_rm = T, digits = 1)),
    
    "Tumor behavior" =
      
      list("In situ"  = ~ n_perc0(BEHAVIOUR == 2, na_rm = T, digits = 1),
           "Invasive" = ~ n_perc0(BEHAVIOUR == 3, na_rm = T, digits = 1),
           "Unknown"  = ~ n_perc0(is.na(BEHAVIOUR), digits = 1)),
    
    "Subtype" =
      
      list("Lobular"  =  ~ n_perc0(SUBTYPE == 1, na_rm = T, digits = 1),
           "Ductal"   =  ~ n_perc0(SUBTYPE == 2, na_rm = T, digits = 1),
           "Tubular"  =  ~ n_perc0(SUBTYPE == 3, na_rm = T, digits = 1),
           "Mixed"    =  ~ n_perc0(SUBTYPE == 4, na_rm = T, digits = 1),
           "Others"   =  ~ n_perc0(SUBTYPE == 5, na_rm = T, digits = 1),
           "Unknown"  =  ~ n_perc0(is.na(SUBTYPE), digits = 1)),
    
    "Estrogen receptor" =
      
      list("Negative"   =  ~ n_perc0(ER == 0, na_rm = T, digits = 1),
           "Positive"  =  ~ n_perc0(ER == 1, na_rm = T, digits = 1),
           "Unknown"   =  ~ n_perc0(is.na(ER), digits = 1)),
    
    "Progesterone receptor" =
      
      list("Negative"  =  ~ n_perc0(PR == 0, na_rm = T, digits = 1),
           "Positive" =  ~ n_perc0(PR == 1, na_rm = T, digits = 1),
           "Unknown"  =  ~ n_perc0(is.na(PR), digits = 1)),
    
    "Grade" =
      
      list("I"      =  ~ n_perc0(SBR == 1, na_rm = T, digits = 1),
           "II"  =  ~ n_perc0(SBR == 2, na_rm = T, digits = 1),
           "III"   =  ~ n_perc0(SBR == 3, na_rm = T, digits = 1),
           "Unknown"          =  ~ n_perc0(is.na(SBR), digits = 1)),
    
    "Stage" =
      
      list("1"    =  ~ n_perc0(STADE == 1, na_rm = T, digits = 1),
           "2"    =  ~ n_perc0(STADE == 2, na_rm = T, digits = 1),
           "3"    =  ~ n_perc0(STADE == 3, na_rm = T, digits = 1),
           "4"    =  ~ n_perc0(STADE == 4, na_rm = T, digits = 1),
           "Unknown"  =  ~ n_perc0(is.na(STADE), digits = 1))
    
  )

st1 <- meta %>% #filter(CT == 1) %>% group_by(MENOPAUSE) %>% 
  summary_table(summ1)
print(st1)
print(st1, cnames = c("Pre-menopausal", "Post-menopausal"))

# ---------------------------------------------------------------

# Models for p-values for table
# Chi-sq and wilcox test
chisq.test(meta$CT, meta$RTHCat1)$p.value
chisq.test(meta$CT, meta$BMICat1)$p.value
chisq.test(meta$CT, meta$SMK)$p.value
chisq.test(meta$CT, meta$DIABETE)$p.value
chisq.test(meta$CT, meta$BP)$p.value
chisq.test(meta$CT, meta$Life_Alcohol_Pattern_1)$p.value
chisq.test(meta$CT, meta$Trait_Horm)$p.value
chisq.test(meta$CT, meta$allaitement)$p.value
chisq.test(meta$CT, meta$CO)$p.value
wilcox.test(meta$ALCOHOL ~ meta$CT, paired = T)$p.value
wilcox.test(meta$DURTHSBMB ~ meta$CT, paired = T)$p.value
wilcox.test(meta$STOCKTIME ~ meta$CT, paired = T)$p.value

# McNemar test for binaries
mcnemar.test(meta$CT, meta$SMK)
mcnemar.test(meta$CT, meta$DIABETE)
mcnemar.test(meta$CT, meta$CO)
mcnemar.test(meta$CT, meta$allaitement)

# MHT use only in post-menopausal
meta.post <- meta[ meta$MENOPAUSE == 1, ]
chisq.test(meta.post$CT, meta.post$Trait_Horm)$p.value
wilcox.test(meta.post$DURTHSBMB ~ meta.post$CT)$p.value

# Pre/post menopausal models
chisq.test(meta$MENOPAUSE, meta$DIAGSAMPLINGCat1)$p.value
chisq.test(meta$MENOPAUSE, meta$BEHAVIOUR)$p.value
fisher.test(meta$MENOPAUSE, meta$SUBTYPE)$p.value
chisq.test(meta$MENOPAUSE, meta$ER)$p.value
chisq.test(meta$MENOPAUSE, meta$PR)$p.value
chisq.test(meta$MENOPAUSE, meta$GRADE)$p.value
fisher.test(meta$MENOPAUSE, meta$STADE)$p.value


# Get individual receptor status for cases
table(meta$ER[meta$CT == 1], useNA = "always") # No missings
table(meta$CERB2[meta$CT == 1], useNA = "always") # 568 missings
table(meta$PR[meta$CT == 1], useNA = "always") # 190 missings


# Identification of triple negative cases
sum(meta$CERB2 == 0 & meta$ER == 0 & meta$PR == 0, na.rm = T)
sum(meta$ER == 0 & meta$PR == 0, na.rm = T)

