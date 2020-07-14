# Breast cancer study multivariate models from unbinned NMR raw data, data from Elodie Jobard 27-6-2019
# Read in data with outliers and QCs (n=1739) or with QCs only
# For multivariate analysis, skip straight to adjmat (line 80) or ROC workspace (line 150)
library(tidyverse)
library(readxl)
library(janitor)
library(MetabolAnalyze)
rawints <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt") %>% select(-8501)

# Metadata with 112 QCs. Assess multivariate quality
metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")

# PCA scores plots (THAW_DATE removes QCs/blanks)
dat1 <- rawints %>% remove_constant() %>% as.matrix
pca <- prcomp(scaling(dat1, type = "pareto"), scale. = F, rank. = 10)
pca2d(pca, group = metaQC$WEEKS)
box(which = "plot", lty = "solid")

library(ggplot2)
ggplot(data.frame(pca$x), aes(PC1, PC2, colour = as.factor(metaQC$TYPE_ECH), 
                              shape = as.factor(metaQC$TYPE_ECH))) + 
  geom_point() + theme_bw() +
  xlab("Score on PC1 (37.6% variance explained)") + ylab("Score on PC2 (25.7% variance explained)") +
  scale_color_discrete(labels = c("Quality controls", "Experimental samples")) +
  scale_shape_discrete(labels = c("Quality controls", "Experimental samples")) +
  theme(legend.position = c(0.83,0.1), legend.title = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggtitle("A")

# Data prep for multivariate models
unscale <- rawints %>% filter(!is.na(metaQC$THAW_DATE)) %>% remove_constant() %>% as.matrix
concs <- scaling(unscale, type = "pareto")

# Calculation of RSDs for chemical shift regions
unscaleQC <- rawints %>% filter(is.na(metaQC$THAW_DATE)) %>% remove_constant() %>% as.matrix
rsd <- apply(unscaleQC, 2, function(x) abs(mean(x))/sd(x))
hist(rsd, breaks = 30)
mean(rsd); median(rsd); IQR(rsd); quantile(rsd)


ggplot(tibble(rsd), aes(x = rsd)) + geom_histogram(colour = "black", fill = "grey80", binwidth = 2) + 
  theme_bw() + xlab("Relative Standard Deviation (%)") + ylab("Number of chemical shift regions") +
  annotate("text", x = 80, y = 1200, label = "Mean RSD = 12.3%") +
  annotate("text", x = 80, y = 1100, label = "Median RSD = 6.9%, IQR = 1.3-18.7%") + ggtitle("B")

# DIAGSAMPLINGCat1, 1: <5y, 2: >5y. DIAGSAMPLINGCat4, same for 2y
meta <- metaQC %>% filter(!is.na(THAW_DATE)) %>%
  select(CT, MATCH, WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, CENTTIME, SAMPYEAR, 
         DIAGSAMPLING, DIAGSAMPLINGCat1, DIAGSAMPLINGCat4, STOCKTIME, DURTHSBMB) %>%
  mutate(DURTHSBMBCat = ifelse(DURTHSBMB > 0, 1, 0), AGEdiag = AGE + DIAGSAMPLING) %>%
  mutate_at(vars(-AGE, -BMI, -CENTTIME, -DIAGSAMPLING, -DIAGSAMPLINGCat1, -DURTHSBMB, 
                 -DIAGSAMPLINGCat4, -AGEdiag), as.factor)


#----

# PC-PR2 to find main sources of variability and remove using residuals method
# (for manuscript; can skip straight to calculation of adjusted matrix)
Z_Meta <- meta %>% select(-MATCH, -CT, -CENTTIME, -DIAGSAMPLINGCat1)

library(pcpr2)
props.raw <- runPCPR2(concs, Z_Meta)
plot(props.raw)
# Greatest sources of variability are BMI > DIABETE > PLACE

# ----

# Adjust for fixed effects only. Random effects model with lme4 did not work, boundary fit or didn't converge
# Note: NAs in CENTTIME, so adjusted matrix loses these 10 observations
adj <- function(x) residuals(lm(x ~ PLACE + WEEKS + #DIABETE + 
                                  FASTING + CENTTIMECat1 + STOCKTIME, data = meta))
adjmat <- apply(concs, 2, adj)
# Final data matrix used for models below
#saveRDS(adjmat, "adjusted_NMR_features.rds")

#----

# PCPR2 of transformed matrix
props.adj <- runPCPR2(adjmat, Z_Meta)

par(mfrow = c(1, 2))
plot(props.raw, main = "Raw feature intensities", font.main = 1)
plot(props.adj, main = "Transformed to residuals of linear model of intensity on confounders*", font.main = 1)

# Can also check PCA 
library(pca3d)
scores.adj <- prcomp(adjmat, scale. = F)
pca2d(scores.adj, group = scores$MENOPAUSE)

#----

# Final data matrix is adjmat, n = 1572

# Multivariate analysis. Subsets to be made:
# 1. All samples; 2. Post-menopausal; 3. Pre-menopausal only
# 4. >55 years at diagnosis; 5. <55 years at diagnosis (AGEdiag)
# 6. <2 year follow-up; 7. >5 year follow-up (tdiag)
adjmat <- readRDS("adjusted_NMR_features_no_age.rds")

# Give controls the time to diagnosis category and age at diagnosis of the corresponding case 
meta <- meta %>% group_by(MATCH) %>% 
  fill(c(DIAGSAMPLINGCat1, DIAGSAMPLINGCat4, AGEdiag), .direction = "downup") %>% filter(!is.na(CENTTIME))

all <- data.frame(class = as.factor(meta$CT), adjmat)

# Logical vectors (need to group for diagnosed > and < 5 years)
post <- meta$MENOPAUSE == 1 
pre  <- meta$MENOPAUSE == 0
diagover55 <- meta$AGEdiag >= 55
diagunder55 <- meta$AGEdiag < 55
followup2 <- meta$DIAGSAMPLINGCat4 == 1
followup5 <- meta$DIAGSAMPLINGCat1 == 2

# Get median follow up times
median(meta$DIAGSAMPLING, na.rm = T)
median(meta$DIAGSAMPLING[post], na.rm = T)
median(meta$DIAGSAMPLING[pre], na.rm = T)
median(meta$DIAGSAMPLING[diagover55], na.rm = T)
median(meta$DIAGSAMPLING[diagunder55], na.rm = T)
median(meta$DIAGSAMPLING[followup2], na.rm = T)
median(meta$DIAGSAMPLING[followup5], na.rm = T)

library(caret)
library(pROC)
# Function to train and fit the model and do a ROC analysis
bc.roc <- function(dat, get.roc = T, ...) {
  
  # Split into training and test sets on class, 75% to training set
  inTrain <- createDataPartition(y = dat$class, p = 0.75, list = F)
  training <- dat[inTrain, ]
  testing <- dat[-inTrain, ]
  
  # Cross validation. If the sample size is large can use a large number of folds (10)
  set.seed(111)
  folds <- createMultiFolds(y = training$class, ...)
  control <- trainControl("repeatedcv", index = folds, selectionFunction = "oneSE")
  print(sapply(folds, length))
  
  # Train PLS model
  mod0 <- train(class ~ ., data = training, method = "pls", metric = "Accuracy", 
                trControl = control, tuneLength = 20) 
  
  plot(mod0)
  print(confusionMatrix(mod0))
  
  # Predict test set
  predictions0 <- predict(mod0, newdata = testing)
  cM <- confusionMatrix(predictions0, reference = testing$class)
  if(get.roc == F) return(cM)
  # Get AUC
  predictions1 <- predict(mod0, newdata = testing, type = "prob")
  result0 <- roc(testing$class, predictions1$`0`, ci = T)
  
}

p0 <- bc.roc(all, k = 10)
p1 <- bc.roc(all[post, ], k = 10)
p2 <- bc.roc(all[pre, ], k = 5, times = 5)
p3 <- bc.roc(all[diagover55, ], k = 10)
p4 <- bc.roc(all[diagunder55, ], k = 5, times = 5)
p5 <- bc.roc(all[followup2, ], k = 5, times = 5)
p6 <- bc.roc(all[followup5, ], k = 10)

# List ROC objects and extract AUCs to a data frame
ll <- list(p0, p1, p2, p3, p4, p5, p6)
aucs <- sapply(ll, "[", 16)
auc.df <- do.call(rbind, auc.ci)

# Get accuracy only (no ROC)
a0 <- bc.roc(all, k = 10, get.roc = F)
a1 <- bc.roc(all[post, ], k = 10, get.roc = F)
a2 <- bc.roc(all[pre, ], k = 5, times = 5, get.roc = F)
a3 <- bc.roc(all[agehi, ], k = 10, get.roc = F)
a4 <- bc.roc(all[agelo, ], k = 10, get.roc = F)
a5 <- bc.roc(all[early, ], k = 10, get.roc = F)

#save.image("ROC_workspace.RData")
load("ROC_workspace.RData")

par(mfrow=c(2,1))
plot.roc(p0, ci = T, grid = T, print.auc = T)
plot.roc(p1, grid = T, print.auc = T)
plot.roc(p2, grid = T, print.auc = T)
plot.roc(p3, grid = T, print.auc = T)
plot.roc(p4, grid = T, print.auc = T)
plot.roc(p5, grid = T, print.auc = T)
plot.roc(p6, grid = T, print.auc = T)

#save.image("ROC_workspace.RData")

library(ggROC)
ggroc(p4, colour = "darkblue", size = 1) + theme_bw() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  ggsave("roc_late.png")

# Compare models
models <- resamples(list("All" = p0, "Post-menopausal" = p1,
                         "Early diagnosis" = p3, "Late diagnosis" = p4))
bwplot(models, metric = "Accuracy")

# Make data frame out of AUCS[CI]
aucs <- t(matrix(c(p0$ci, p1$ci, p2$ci, p3$ci, p4$ci, p5$ci, p6$ci), nrow = 3))
library(dplyr)
df <- data.frame(aucs[, c(2,1,3)])
library(metafor)

par(mar=c(5,4,1,2))
forest(x = df$X1, ci.lb = df$X2, ci.ub = df$X3, slab = c(
  "All study subjects", "Post-menopausal", "Pre-menopausal", "Diagnosis < 5 y",
  "Diagnosis > 5 y", "Pre + Post HT users", "Post non-HT users"),
  refline = 0.5, xlab = "ROC AUC", cex = 1, pch = 22
)
hh <- par("usr")
text()
text(hh[1], 9, "Study group", pos = 4, cex = 1)
text(hh[2], 9, "AUC [95% CI]", pos = 2, cex = 1)



# Description of files

# 1694 obs. of 8501 NMR variables (outliers removed, one NA variable in 8501st col)
#raw <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt") %>% select(-8501)

# 1739 obs. of 8500 NMR variables (all samples, QCs, blanks)
#raw1 <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";")

# Metadata are stored in the following Excel file. The third sheet has the metadata w/o outliers
# (hope it's in the same order as the NMR data!)
#metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")






