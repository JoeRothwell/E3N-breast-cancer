# Breast cancer study multivariate models from unbinned NMR raw data, data from Elodie Jobard 27-6-2019
# Read in data with outliers and QCs (n=1739) or with QCs only

library(tidyverse)
library(readxl)
rawints <- read_tsv("1510_XAlignedE3NcpmgssCitPEGfinal.txt") %>% select(-8501)
bc.meta <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")

# Dataset containing all samples, QCs and blanks
#mat <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";", n_max = 1738) 

# Subset samples only and plot PCA with pareto scaling
prep.data <- function(dat, meta, drop.qc = T, get.pca = T) {
    
  # Subset samples only (removing 112 QCs) from data and metadata
  samp <- !is.na(meta$THAW_DATE)
  mat  <- if(drop.qc == T) dat[samp, ] %>% as.matrix else as.matrix(dat)
  meta <- if(drop.qc == T) meta[samp, ]

  #anyNA(mat)
  #sum(apply(mat, 2, anyNA))
  
  # remove zero variance columns
  zerovar <- sum(apply(mat, 2, var) == 0)
  print(paste("There are", zerovar, "zero variance columns"))
  
  # There are 1116. Place in logical vector
  nonzerovar <- apply(mat, 2, var) != 0
  which(nonzerovar)
  
  mat0 <- mat[ , nonzerovar]
  print(paste("Dimensions", dim(mat0)))

  # Scale and run PCA  
  library(MetabolAnalyze)
  scalemat <- scaling(mat0, type = "pareto")
  
  if(get.pca == F) return(list(scalemat, meta))
  pca <- prcomp(scalemat, scale. = F, center = T, rank. = 10)
  output <- data.frame(pca$x) %>% bind_cols(meta)
}

scores <- prep.data(rawints, bc.meta, get.pca = T)

# Output PCA scores plot for manuscript
library(ggplot2)
ggplot(scores, aes(PC1, PC2, colour = as.factor(TYPE_ECH))) + geom_point() + theme_bw() +
  xlab("Score on PC1") + ylab("Score on PC2") +
  scale_color_discrete(labels = c("Experimental samples", "QCs")) +
  theme(legend.position = "bottom", #legend.justification = c(0, 0), 
        legend.title = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")

library(pca3d)
pca2d(as.matrix(scores[, 1:10]), group = scores$WEEKS)
box(which = "plot", lty = "solid")

# Output Pareto-scaled intensities and corresponding metadata for multivariate analysis
dat <- prep.data(rawints, bc.meta, get.pca = F)

concs <- dat[[1]]
meta <- dat[[2]] %>%
  select(CT, MATCH, WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, CENTTIME, SAMPYEAR, 
         DIAGSAMPLINGCat1, STOCKTIME, DURTHSBMB) %>%
  mutate_at(vars(-AGE, -BMI, -CENTTIME, -DIAGSAMPLINGCat1, -DURTHSBMB), as.factor) %>%
  mutate(DURTHSBMBCat = ifelse(DURTHSBMB > 0, 1, 0))

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
adj <- function(x) residuals(lm(x ~ PLACE + WEEKS + DIABETE + FASTING + CENTTIMECat1 + 
                                  STOCKTIME, data = meta))
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
# 1. All samples; 2. Pre-menopausal only; 3. Post-menopausal only; 4. Diagnosed < 5 years only; 
# 5. Diagnosed > 5 years only; 6. Pre-menopausal or no HT; 7. Post-menopausal and HT.
#adjmat <- readRDS("adjusted_NMR_features.rds")
#adjmat <- readRDS("adjusted_NMR_features_no_age.rds")

# First give each control the time to diagnosis time of the corresponding case
meta <- meta %>% group_by(MATCH) %>% mutate(tdiag = max(as.numeric(DIAGSAMPLINGCat1), na.rm = T)) %>%
  filter(!is.na(CENTTIME))

all <- data.frame(class = as.factor(meta$CT), adjmat)

# Logical vectors (need to group for diagnosed > and < 5 years)
pre   <- meta$MENOPAUSE == 0
post  <- meta$MENOPAUSE == 1 
early <- meta$tdiag == 1
late  <- meta$tdiag == 2
young <- meta$AGE < 55
old   <- meta$AGE >= 55
#noHT  <- meta$DURTHSBMBCat == 0
#HT    <- meta$DURTHSBMBCat == 1

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
  
  #plot(mod0)
  confusionMatrix(mod0)
  
  # Predict test set
  predictions0 <- predict(mod0, newdata = testing)
  
  cM <- confusionMatrix(predictions0, reference = testing$class)
  if(get.roc == F) return(cM)
  
  # Get AUC
  predictions0_1 <- predict(mod0, newdata = testing, type = "prob")
  result0 <- roc(testing$class, predictions0_1$`0`, ci = T)
  
}

p0 <- bc.roc(all, k = 10)
p1 <- bc.roc(all[post, ], k = 10)
p2 <- bc.roc(all[pre, ], k = 5, times = 5)
p3 <- bc.roc(all[early, ], k = 10)


a0 <- bc.roc(all, k = 10, get.roc = F)
a1 <- bc.roc(all[post, ], k = 10, get.roc = F)
a2 <- bc.roc(all[pre, ], k = 5, times = 5, get.roc = F)
a3 <- bc.roc(all[early, ], k = 10, get.roc = F)



p4 <- bc.roc(all[late, ], k = 10)
#p5 <- bc.roc(all[noHT, ], k = 10)
#p6 <- bc.roc(all[HT, ], k = 10)


par(mfrow=c(2,2))
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

# 1739 obs. of 8500 NMR variables (all samples)
#raw1 <- read_delim("X_AlignedE3NData_cpmg_ssCitPEG_1112.txt", delim = ";")

# Metadata are stored in the following Excel file. The third sheet has the metadata w/o outliers
# (hope it's in the same order as the NMR data!)
#metaQC <- read_xlsx("1510_MatriceY_CohorteE3N_Appar.xlsx", sheet = 4, na = ".")





