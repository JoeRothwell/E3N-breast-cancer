# Correlation ethanol - alcohol intake
source("BC_prep_data.R")
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
