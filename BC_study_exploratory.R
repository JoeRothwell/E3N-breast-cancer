# Breast cancer NMR metabolomics data exploratory analysis
# Scaled and unscaled data (already done in SIMCA). Updated after receiving scaled data.
library(tidyverse)
library(readxl)
library(gplots)

# Read unscaled data and metadata (each 1882 obs)
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
ints <- scale(ints0)
meta <- read.csv("metadata.csv", na = "9999")


# Some exploratory analyses----
# Follow-up time plot for paper
meta$menopause <- ifelse(meta$MENOPAUSE == 1, "Post-menopausal", "Pre-menopausal")

library(ggplot2)
plot1 <- ggplot(meta, aes(x = DIAGSAMPLING)) + 
  geom_histogram(colour = "black", fill = "grey80")+ theme_minimal() + 
  scale_fill_manual(values = c("white", "grey80")) + 
  xlab("Time from blood collection to diagnosis (years)") + ylab("No. cases") +
  theme(legend.position = c(0.8, 0.8)) +
  facet_grid(.~menopause, scales = "free_y") + ggtitle("")

# Correlation heatmap for paper
library(corrplot)
colnames(ints)[1] <- "3-Hydroxybutyrate"
colnames(ints)[41] <- "N-acetyl glycoproteins"
cormat <- cor(ints, use = "pairwise.complete.obs")

# Corrplot and corrr
colnames(cormat) <- rep("", 43)
#rownames(cormat) <- NULL
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.7,  tl.srt = 30,
         hclust.method = "ward", order = "hclust", type = "full")

library(corrr)
rplot(cormat, shape = 15)

# Dendrogram
library(dendextend)
dend <- hclust(dist(cormat)) %>% as.dendrogram
par(mar=c(1,1,1,8))
dend %>% 
  set("labels_col", value = c("skyblue", "orange", "grey"), k = 3) %>%
  set("branches_k_color", value = c("skyblue", "orange", "grey"), k = 3) %>%
  plot(horiz = T, axes = F)


# Read scaled data and subset samples from meta
ints.all <- read.delim("1507_XMetabolite_std_cpmg_E3N.txt")

which(apply(as.matrix(ints0), 2, min) < 0)
# 3 compounds have values < 0: formate, hypoxanthine, inosine

# subset IDs to get subjects included in CC. Get positions of final CC samples in metadata
samples <- ints.all$CODBMB %in% meta$CODBMB
ints <- ints.all[samples, ]

# Distributions all data
par(mfrow = c(2,1))
hist(as.matrix(ints0), breaks = 50, col = "dodgerblue")
hist(as.matrix(ints),  breaks = 50, col = "dodgerblue")

# Plot metabolites individually
library(RColorBrewer)

#palette(rainbow(5))
palette(brewer.pal(5, "Dark2"))
plot.ts(ints[, 1:10], type = "p", col = meta$STOCKTIME, main = "Hydroxybutyrate - Creatine")
plot.ts(ints[, 11:20], type = "p", col = meta$STOCKTIME, main = "Creatinine - Glutamate")
plot.ts(ints[, 21:30], type = "p", col = meta$STOCKTIME, main = "Glycerol - Isoleucine")
plot.ts(ints[, 31:40], type = "p", col = meta$STOCKTIME, main = "Lactate - Malonate")
plot.ts(ints[, 41:ncol(ints)], type = "p", col = meta$STOCKTIME, main = "Malonate - Succinate")

# Problem compounds only
plot.ts(ints[, c(13, 14, 23, 24, 39, 43)], type = "p", col = meta$SAMPYEAR)

# Or use walk2 to go through plotting all columns
par(mfrow = c(5, 1), mai = c(0.3, 0.5, 0.2, 0.1))
walk2(ints, colnames(ints), ~ plot(.x, main = .y, col = meta$SAMPYEAR))
walk2(ints, colnames(ints), ~ boxplot2(.x ~ meta$MENOPAUSE, main = .y, top = T))

# Plot distributions by compound and menopause
ints1 <- cbind(meno = meta$MENOPAUSE, ints) %>% data.frame
ints.melt <- gather(ints1, compound, value, -meno)
ggplot(ints.melt, aes(as.factor(meno), value)) + geom_boxplot(outlier.shape = NA) + theme_bw() +
  geom_jitter(width=0.3, alpha = 0.3, size = 0.1) + facet_wrap( ~ compound, ncol = 8, scales = "free_y")

# For baseline characteristics table: see BC_baseline_char.R

# Median intensities before scaling
dev.off()
Intensity <- apply(ints0, 2, median)
plot(Intensity, col = "white", main = "Median compound intensities")
text(Intensity, labels = colnames(ints0))
# Fatty acids and glucose strongest signals, then lactate, glycerophosphocholine.

# Visualise case-control differences
ggplot(ints, aes(x= as.factor(meta$CT), y=log(Hypoxanthine))) + 
  geom_line(aes(group = meta$MATCH), alpha = 0.5) + 
  geom_point(aes(group = meta$MATCH), alpha = 0.5) + 
  theme_minimal() + facet_grid( ~ meta$MENOPAUSE)


# Run PCA of all samples
pca0 <- prcomp(ints0, scale. = F, center = F)
biplot(pca0)

pca <- prcomp(ints[, -1], scale.=F)
biplot(pca, cex = 0.7)
scores <- data.frame(pca$x)
ggplot(scores, aes(x = PC1, y = PC2, col = as.factor(meta$MENOPAUSE))) + 
  geom_point() + scale_colour_manual(values = c("black", "grey"))



# Plot 
library(pca3d)
par(mfrow = c(1, 2))
pca2d(pca)
title("A", font.main = 1, adj = 0)
box(which = "plot", lty = "solid")

# Run PCPR2 to compare sources of variability
library(pcpr2)
alldata <- inner_join(meta, ints, by = "CODBMB")
X_DataMatrixScaled <- select(alldata, `3Hydroxybutyrate`:Succinate) %>% as.matrix
Z_Meta <- alldata %>%
  select(WEEKS, PLACE, AGE, BMI, MENOPAUSE, FASTING, SMK, DIABETE, CENTTIMECat1, SAMPYEAR, STOCKTIME) %>%
  mutate_at(vars(-AGE, -BMI, -WEEKS, -STOCKTIME), as.factor)

par(mar=c(6,5,4,2))
props <- runPCPR2(X_DataMatrixScaled, Z_Meta)
plotProp(props, main = "B", font.main = 1, adj = 0)

# PCA of subject IDs only
ints1 <- alldata %>% select(`3Hydroxybutyrate`:Succinate)
pca1 <- prcomp(ints1, scale.=F)

par(mfrow = c(1,2))
plt <- pca2d(pca1, group = alldata$CT)
title("Metabolite profiles of 1582 samples", font.main = 1)
box(which = "plot", lty = "solid")
legend("topleft", legend = plt$groups, col=plt$colors, pch=plt$pch)

# Adjust using residuals method
adj <- function(x) residuals(lm(x ~ BMI + SMK + DIABETE, data = alldata))
adjmat <- apply(ints1, 2, adj)

# Repeat PCA
pca2 <- prcomp(adjmat, scale. = F)
pca2d(pca2, group = alldata$CT)
title("Residuals-adjusted metabolite\nprofiles of 1582 samples", font.main = 1)
box(which = "plot", lty = "solid")
legend("topleft", legend = plt$groups, col=plt$colors, pch=plt$pch)

# Test datasets for scaling -----------------------

library(tidyverse)
library(readxl)

# Original data, 1623 observations of 44 compounds - appears to be Pareto scaled
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")

# Update 11/10/19: Updated with unscaled data (will scale to unit variance). Added "_unscaled" to filename.
# Note: NAC1 and NAC2 are merged and only 1582 observations
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt") #%>% as.matrix
ints1 <- scale(ints0)

# Lifestyle data. Subset variables needed
meta <- read_csv("Lifepath_meta.csv", na = "9999")
samples <- ints$CODBMB %in% meta$CODBMB
ints <- ints[samples, -1]

# Exploratory of different datasets. Find differences between received datasets (formerly BC_compounds_scaling)
library(MetabolAnalyze)

# 1. Original data (1507); 2. New data (no scaling); 3. New data (unit scaled); New data (Pareto scaled)
new     <- ints0 %>% as.matrix
unit    <- scale(ints0)
pareto1 <- scaling(ints0, type = "pareto")

ll <- lapply(list(ints, new, unit, pareto1), function(x) prcomp(x, scale. = F, center = F))
par(mfrow = c(2,2))
lapply(ll, function(x,y) { pca2d(x)
  title(main = paste("Scores plot"))
  box(which = "plot")
})

#how to make hotpink, limegreen, orange, dodgerblue and give different titles?

# Conclusion: the old dataset seems to be a unit variance scaled version of the new dataset
# (although they are not exactly the same; the new datasets has merged NAC1 and NAC2)
# Decision: use unscaled dataset and apply unit scaling, discard old scaled dataset

# Moved from BC_dendrograms.R

library(corrplot)
cormat <- cor(ints1)
colnames(cormat) <- NULL
corrplot(cormat, method = "square", tl.col = "black", tl.cex = 0.8)

dend1 <- cormat %>% dist %>% hclust %>% as.dendrogram

plot(dend0)
plot(dend1)

library(dendextend)
dend1 %>% set("labels_col", "blue") %>% hang.dendrogram() %>% 
  plot(main = "Change label's color")

dend1 %>% set("labels_cex", 1) %>% set("labels_col", value = c(3,4), k=2) %>% 
  plot(main = "Color labels \nper cluster")

dend1 %>% set("branches_k_color", k = 3) %>% plot(main = "Nice defaults")

library(circlize)


rownames(cormat) <- abbreviate(rownames(cormat), named = F)
dend <- cormat %>% dist %>% hclust %>% as.dendrogram %>%
  set("branches_k_color", k=3) %>% set("labels_cex", c(1)) #%>%
#hang.dendrogram(hang_height = 0.3)

circos.par(gap.after = 0) #gap.degree for gaps between sectors

circlize_dendrogram(dend, dend_track_height = 0.85)



# Aggregation of NMR features----
# Formerly NMR_aggregation.R

# Calculation of compound table by aggregation of 107 NMR features corresponding to 43 compounds
# Raw unscaled data received by email 11/10/2019

# 3 compounds have values < 0: formate, hypoxanthine, inosine
# 130 clusters were used to make the 43 compounds
library(readxl)
library(tidyverse)
library(reshape2)

# Get cluster data
dat <- read_xlsx("1510_ClusterAnnotation_E3N.xlsx", sheet = 3, skip = 1)
dat0 <- data.frame(t(dat)) %>% rownames_to_column("cluster")

# Get compound names (with spaces)
cmpds <- read_xlsx("1510_ClusterAnnotation_E3N.xlsx", sheet = 3, col_names = F, n_max = 1) %>% 
  gather() %>% select(value)

# Bind together, drop the sums to leave just the features
feats <- cbind(cmpds, dat0) %>% 
  fill(value) %>%
  filter(!str_detect(cluster, "SUM")) %>%
  select(-cluster)

# Get minimum values
ff <- which(apply(feats[, -1], 1, min) < 0)
feats$value[ff]

# Recalculate sums from groups
cmpds <- feats %>% group_by(value) %>% summarise_all(sum) %>% column_to_rownames("value")
# Get minimum values
which(apply(cmpds, 1, min) < 0)


# Correlation ethanol - alcohol intake----
# Formerly cor-ethanol_alcohol
library(tidyverse)
meta <- read_csv("metadata.csv", na = "9999")
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")

### Categorical analysis of plasma ethanol, no points excluded. Plot of unscaled data
plot(meta$ALCOHOL, ints0$Ethanol)
plot(log(meta$ALCOHOL), ints0$Ethanol)
ethanol.cat <- cut_number(ints0$Ethanol, n = 4, labels = 1:4)

fit <- lm(log(meta$ALCOHOL+0.1) ~ ethanol.cat)
fit1 <- lm(log(meta$ALCOHOL+0.1) ~ ethanol.cat + meta$BMI + meta$MENOPAUSE)
# No difference in alcohol intake between Q1 and Q4 ethanol

# Check for weekday/weekend differences
library(lubridate)
meta$sampdate2 <- dmy(meta$SAMPDATE)
meta$sampday <- weekdays(meta$sampdate2)
table(meta$sampday)

boxplot(ints0$Ethanol ~ sampday, data = meta, ylim = c(0.02, 0.05))
library(ggplot2)
ggplot(meta, aes(y = ints0$Ethanol, x = sampday)) + geom_boxplot() +
  geom_jitter(width = 0.2, size = 0.5) + theme_bw()

# Continuous analysis of plasma ethanol
# Remove top and bottom 1% of each compound and replace with NA
outliers <- function(x) x > quantile(x, probs = 0.99) | x < quantile(x, probs = 0.01)
logicalmat <- apply(ints0, 2, outliers)
ints0[logicalmat] <- NA

# Scale to unit variance
ints <- scale(ints0)

plot(meta$ALCOHOL, ints[, 14])
plot(log(meta$ALCOHOL), ints[, 14], xlab = "log(alcohol, g/day)", ylab = "Ethanol concentration",
     main = "Correlation R-squared = -0.036" )

vec <- !is.na(ints[, 14])
fit <- lm(log(meta$ALCOHOL + 0.001) ~ ints[, 14])
summary(fit)

fit1 <- lm(meta$ALCOHOL ~ ints[, 14])
summary(fit1)

cor(log(meta$ALCOHOL + 0.001), ints[, 14], use = "pairwise.complete.obs")
cor.test(log(meta$ALCOHOL + 0.001), ints[, 14], use = "pairwise.complete.obs")


# Get a subset of the most discriminating compounds and plot PCA
# Histidine, NAC, glycerol, ornithine, ethanol, pyruvate, albumin, glutamate, glutamine, 3 fatty acids
ss <- ints[, c(5,14,15,16,17,19,20,27,28,33,41,42)]
library(zoo)
ss1 <- na.aggregate(ss, function(x) median(x))
pca <- prcomp(ss1, scale. = F)$x

dat <- cbind(meta, pca)
ggplot(dat, aes(DIAGSAMPLING, PC2, colour = as.factor(MENOPAUSE))) + geom_point() + 
  theme_bw() + scale_colour_manual(values = c("black", "grey"))

ggplot(dat, aes(PC1, PC2, colour = as.factor(MENOPAUSE))) + geom_point() +
  theme_bw()

biplot(pca)


# Barplots alc/eth intake----
# Formerly script Barplots_alcohol.R
library(tidyverse)

# Read 1623 observations of 44 intensity variables (appears to be final scaled data) and metadata
#ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")
ints <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt") %>% scale()

# Lifestyle data. Subset variables needed
meta <- read_csv("metadata.csv", na = "9999")
#data <- left_join(meta, ints, by = "CODBMB")
data <- bind_cols(meta, data.frame(ints))

# Barplot of alcohol intake g/day by menopausal status
mat <- data %>% group_by(CT, MENOPAUSE) %>% summarise(alcohol = mean(ALCOHOL)) %>% 
  spread(MENOPAUSE, alcohol) %>% as.matrix

library(gplots)
barplot2(mat[ , -1], beside = T, col = c("grey12", "grey82"), ylim = c(0, 20),
         names.arg = c("Pre-menopausal", "Post-menopausal"), legend = c("Controls", "Cases"),
         ylab = "Alcohol intake g/day")

pval.pre <- wilcox.test(ALCOHOL ~ CT, data = data, subset = MENOPAUSE == 0)$p.value
pval.post <- wilcox.test(ALCOHOL ~ CT, data = data, subset = MENOPAUSE == 1)$p.value

# Boxplot and barplot of alcohol by menopausal status
boxplot(data$Ethanol ~ data$CT + data$MENOPAUSE, varwidth = T, outline = F,
        names = c("Control, pre", "Case, pre", "Control, post", "Case, post"),
        col = "dodgerblue", ylab = "Plasma ethanol conc (scaled)")

mat <- data %>% group_by(CT, MENOPAUSE) %>% summarise(ethanol = mean(Ethanol)) %>% 
  spread(MENOPAUSE, ethanol) %>% as.matrix

barplot2(mat[ , -1], beside = T, col = c("grey12", "grey82"), #plot.ci = T,
         names.arg = c("Pre-menopausal", "Post-menopausal"), legend = c("Controls", "Cases"))


# Correlation ethanol - alcohol intake
plot(data$ALCOHOL, data$Ethanol, xlab = "alcohol, g/day", ylab = "Ethanol concentration", main = "r = -0.036" )

# Subset ethanol data from metabolomics data
#cor.test(meta$ALCOHOL, ints[, 14], use = "pairwise.complete.obs")
cor.test(data$ALCOHOL, data$Ethanol, use = "pairwise.complete.obs", method = "pearson")

data$alcQ4 <- cut_number(data$ALCOHOL, n=4, label = c("Q1", "Q2", "Q3", "Q4"))
boxplot(Ethanol ~ alcQ4, data = data)
fit <- kruskal.test(Ethanol ~ alcQ4, data = data) # 0.7
fit.pre <- kruskal.test(Ethanol ~ alcQ4, data = data[data$MENOPAUSE == 0, ]) # 0.55
fit.pos <- kruskal.test(Ethanol ~ alcQ4, data = data[data$MENOPAUSE == 1, ]) # 0.86

# Plot boxplot of intakes for manuscript panel plot
plotB <- ggplot(data, aes(x = alcQ4, y = Ethanol, fill = as.factor(MENOPAUSE))) + 
  geom_boxplot(outlier.shape = NA, #fill = "grey", 
               position = position_dodge()) +
  scale_fill_manual(values = c("white", "grey")) +
  #geom_jitter(width = 0.2, colour = "black", size = 0.8) + 
  geom_point(position = position_jitterdodge(), size = 0.8) +
  annotate("text", label = "Kruskal-Wallis P = 0.70", x = 2.5, y = 1.5, size = 3) +
  ylim(-1.0, 1.1) + 
  theme_bw(base_size = 10) +
  xlab("Category reported alcohol intake (quartile)") +
  ylab("Plasma ethanol (scaled intensity)") +
  theme(panel.grid.minor = element_blank()) +
  ggtitle("")

# Correlated with outliers but correlation goes away after removing these

# Get a subset of the most discriminating compounds and plot PCA
# Histidine, NAC, glycerol, ornithine, ethanol, pyruvate, albumin, glutamate, glutamine, 3 fatty acids
ss <- ints[, c(5,14,15,16,17,19,20,27,28,33,41,42)]
library(zoo)
ss1 <- na.aggregate(ss, function(x) median(x))
pca <- prcomp(ss1, scale. = F)$x

dat <- cbind(meta, pca)
ggplot(dat, aes(DIAGSAMPLING, PC2, colour = as.factor(MENOPAUSE))) + geom_point() + 
  theme_bw() + scale_colour_manual(values = c("black", "grey"))

ggplot(dat, aes(PC1, PC2, colour = as.factor(MENOPAUSE))) + geom_point() +
  theme_bw()

biplot(pca)

# BC alcohol associations----
# Formerly BC_alcohol_association.R
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



# Association BC risk factors - metabolites
library(tidyverse)
meta <- read_csv("metadata.csv", na = "9999")
ctrl <- meta$CT == 0 #& meta$MENOPAUSE == 1

# Get control subjects and remove outliers from metabolite data
ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")[ctrl, ]
outliers <- function(x) x > quantile(x, probs = 0.99) | x < quantile(x, probs = 0.01)
logicalmat <- apply(ints0, 2, outliers)
ints0[logicalmat] <- NA

# Scale to unit variance
ints <- scale(ints0)
meta1 <- meta %>% filter(CT == 0) %>% select(FASTING, SMK, ALCOHOL, Life_Alcohol_Pattern_1, 
                                             BMI, AGE, BP, TTAILLE, Trait_Horm, MENOPAUSE, DIABETE, CO, DURTHSBMB, CENTTIME, STOCKTIME,
                                             WEEKS)

# Metabolite-covariate associations----
# Formerly script Associations_covariates
# Risk factors to be tested
#FASTING, SMK, BMI, AGE, BP, TTAILLE, Trait_Horm, MENOPAUSE, DIABETE, CO
library(broom)
lm.meno    <- function(x) lm(MENOPAUSE ~ x + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.fasting <- function(x) lm(FASTING ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.alc     <- function(x) lm(ALCOHOL ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.lifealc <- function(x) lm(Life_Alcohol_Pattern_1 ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.smk     <- function(x) lm(SMK ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.bmi     <- function(x) lm(BMI ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.age     <- function(x) lm(AGE ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.bp      <- function(x) lm(BP ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.ttaille <- function(x) lm(TTAILLE ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.horm    <- function(x) lm(Trait_Horm ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.diabet  <- function(x) lm(DIABETE ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)
lm.co      <- function(x) lm(CO ~ x + MENOPAUSE + DURTHSBMB + CENTTIME + STOCKTIME, data = meta1)

fits0 <- apply(ints, 2, lm.meno) %>% map_df(tidy)
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

all <- bind_rows(fits0, fits1, fits1a, fits1b, fits2, fits3, fits4, fits5, fits6, fits7,
                 fits9, fits10) %>% filter(str_detect(term, "x"))  %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  bind_cols(cmpd = rep(colnames(ints), 12))

pFDR <- all %>% filter(p.adj < 0.05) %>% select(p.value) %>% max

library(RColorBrewer)
# Different colour palettes
set2 <- rep(brewer.pal(8, "Set2"), each = 43, length.out = nrow(all))
accent <- rep(brewer.pal(8, "Accent") , each = 43, length.out = nrow(all))
rain <- rep(rainbow(12), each = 43, length.out = nrow(all))

# Plot (add col palette as necessary)
#plot(-log10(all$p.value), col = set2, pch = 19, cex = 0.6)

# x axis width
x = 1:nrow(all)

# draw empty plot
plot(NULL, xlim=c(0, nrow(all)), ylim=c(0, max(-log10(all$p.value))), xaxt='n',
     ylab='-log10(p-value)', xlab='')
points(x, -log10(all$p.value), pch=19, col=set2, cex = 0.6)

# axis labels
labs <- c("Meno", "Fast", "Alc", "Life alc.", "Smoke", "BMI", "Age", "BP", 
          "WC", "HRT", "Diabet.", "OC")
abline(h = -log(0.05), col = "grey")
abline(h = -log(pFDR), col = "red")

# Get label positions and add to axis
midpoint <- ncol(ints)/2
at.labs <- seq(midpoint, 2*midpoint*length(labs), 2*midpoint)
axis(1, at = at.labs, labels = labs, las=1, cex.axis = 0.7)

# Make variable for annotations and add to plot
all1 <- all %>% mutate(p.low = ifelse(-log10(p.value) > 12, p.value, NA))
text(1:nrow(all), -log10(all1$p.low), all$cmpd, pos = 4, cex = 0.7)


# Old ------------------------------

# Time to centrifugation vs fasting status
boxplot(CENTTIME ~ FASTING, data = meta, varwidth = T)
meta1 <- meta[meta$CENTTIME < 100, ]
library(gplots)
boxplot2(CENTTIME ~ FASTING, data = meta1, varwidth = T, col = "dodgerblue",
         xlab = "Fasting status", ylab = "TBC")

# Rawest feature data
raw <- read_delim("D:/J_ROTHWELL/X_AlignedCohorteE3NData_cpmg_ssCitPEG_0612.txt", delim = ";")

# Excel files with intermediate steps
# List of 54 compounds and IDs
#xl1 <- read_xlsx("C:/J_ROTHWELL/1505_E3N_Identification.xlsx")

# Previous version of intensities above without scaling
#dat2 <- read_tsv("C:/J_ROTHWELL/1507_XMetaboliteE3N_cpmg.txt")

# List of metabolites and corresponding clusters
#xl2 <- read_xlsx("C:/J_ROTHWELL/1507_ClusterAnnotation_E3N.xlsx")

# Metadata and intensities
#xl4 <- read_xlsx("C:/J_ROTHWELL/1603_MatriceY_CohorteE3N_Appar.xlsx", sheet = 6)

#ints <- read_tsv("C:/J_ROTHWELL/1507_XMetabolite_std_cpmg_E3N.txt")
#ints <- ints[, -1]
#ints1 <- read_tsv("C:/J_ROTHWELL/1507_XMetaboliteE3N_cpmg.txt", skip = 1)

# Food intake data
#path <- "Y:/RepPerso/Fabienne WILM/02_Demandes_ponctuelles/10_LIFEPATH/TABLES"
#list.files(path)

# Food intake and other data
#meta1 <- read_csv("D01_20171031_LIFEPATH.csv")
#meta2 <- read_csv("D01_20161018_LIFEPATH.csv")
#meta3 <- read_csv("D01_20150917_LIFEPATH.csv")



