# Barplots alcohol intake and ethanol
library(tidyverse)

# Read 1623 observations of 44 intensity variables (appears to be final scaled data) and metadata
ints <- read_tsv("1507_XMetabolite_std_cpmg_E3N.txt")

# Lifestyle data. Subset variables needed
meta <- read_csv("metadata.csv", na = "9999")
data <- left_join(meta, ints, by = "CODBMB")
mat <- data %>% group_by(CT, MENOPAUSE) %>% summarise(alcohol = mean(ALCOHOL)) %>% 
  spread(MENOPAUSE, alcohol) %>% as.matrix

library(gplots)
barplot2(mat[ , -1], beside = T, col = c("grey12", "grey82"), ylim = c(0, 20),
         names.arg = c("Pre-menopausal", "Post-menopausal"), legend = c("Controls", "Cases"))

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
fit <- kruskal.test(Ethanol ~ alcQ4, data = data)

# Plot boxplot of intakes for manuscript panel plot
plotB <- ggplot(data, aes(x = alcQ4, y = Ethanol)) + geom_boxplot(outlier.shape = NA, fill = "grey") +
  geom_jitter(width = 0.2, colour = "black", size = 0.8) + 
  annotate("text", label = "Kruskal-Wallis P = 0.70", x = 2.5, y = 1.5, size = 3) +
  ylim(-1.5, 1.5) + 
  theme_bw(base_size = 10) +
  xlab("Category reported alcohol intake (quartile)") +
  ylab("Plasma ethanol (scaled intensity)") +
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
