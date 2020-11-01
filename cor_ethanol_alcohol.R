# Correlation ethanol - alcohol intake
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
