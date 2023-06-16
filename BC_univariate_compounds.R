# BC risk models for metabolites and also for lifestyle variables alone
# Also heatmap of differences
source("BC_prep_data.R")
library(broom)
library(survival)

# Univariate analyses & MS tables----

# CLR models to get odds ratios for metabolites: all, pre-menopausal only and post-menopausal only
# All subjects
fits1 <- apply(ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSBMB + 
         DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta))

# Pre-menopausal. Need to remove hormone treatment therapy variable
fits2 <- apply(ints[pre, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + #DURTHSBMB + 
                     CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[pre, ]))

# Pre-menopausal + adjusting for ethanol
meta.eth <- cbind(meta[pre, ], ETH = ints[pre, 14])
fits3 <- apply(ints[pre, -14], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + ETH +
                     CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta.eth))

# Post-menopausal
fits4 <- apply(ints[post, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSBMB + 
                     CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[post, ]))

# Fasting only for review Dec 2020
fits1a <- apply(ints[fast, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL +
                    DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[fast, ]))

fits1b <- apply(ints[fast0, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL +
                    CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[fast0, ]))

fits1c <- apply(ints[fast1, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL +
                    CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[fast1, ]))

# Sensitivity analysis for supplement table 5. Adjust Life_Alcohol_Pattern2, exclude cases <55 at diagnosis,
# exclude cases diagnosed first 2 years of follow up
fits5 <- apply(ints[pre, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
                     Life_Alcohol_Pattern_2 + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[pre, ]))

fits6 <- apply(ints[preS2, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
                     CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[preS2, ]))

fits7 <- apply(ints[pre0, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
                    CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[pre0, ]))

fits8 <- apply(ints[pos, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
                    CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[pos, ]))

# Follow-up time all subjects by quartile
fits8 <- apply(ints[fuQ1, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
              CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[fuQ1, ]))

fits9 <- apply(ints[fuQ2, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
              CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[fuQ2, ]))

fits10 <- apply(ints[fuQ3, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
              CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[fuQ3, ]))

fits11 <- apply(ints[fuQ4, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
              CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[fuQ4, ]))

# Duration of OC use
fits1oc <- apply(ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSBMB + 
              DURTHSBMB + CENTTIME + STOCKTIME + durOC + strata(MATCH) + x, data = meta))

fits2oc <- apply(ints[pre, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + #DURTHSBMB + 
                CENTTIME + STOCKTIME + durOC + strata(MATCH) + x, data = meta[pre, ]))

# Post-menopausal
fits4oc <- apply(ints[post, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSBMB + 
                CENTTIME + STOCKTIME + durOC + strata(MATCH) + x, data = meta[post, ]))

# Tables for manuscript
# Generate tidy output table from models
tidy.output <- function(mod) {
  
  library(broom) 
  # Function to exponentiate and round
  
  df <- map_df(mod, tidy) %>% filter(term == "x") %>% cbind(Compound = names(mod)) %>%
    left_join(cmpd.meta, by  = "Compound") %>%
    mutate(P.value = round(p.value, 4), FDR = round(p.adjust(p.value, method = "fdr"), 3))
  
  # Calculate FDR adjusted confidence intervals  
  fdr.hits <- sum(p.adjust(df$p.value, method = "fdr") < 0.05)
  alpha.raw <- 0.05
  alpha.fdr <- (fdr.hits * 0.05)/nrow(df)
  
  df$OR <- exp(df$estimate)
  df$ci.low <- exp( df$estimate - (df$std.error * qnorm(1 - (alpha.raw/2))) )
  df$ci.high <- exp( df$estimate + (df$std.error * qnorm(1 - (alpha.raw/2))) )
  
  if(fdr.hits > 0) {
    df$ci.low.fdr <- exp( df$estimate - (df$std.error * qnorm(1 - (alpha.fdr/2))) )
    df$ci.high.fdr <- exp( df$estimate + (df$std.error * qnorm(1 - (alpha.fdr/2))) )
  } else {
    df$ci.low.fdr <- df$ci.low
    df$ci.high.fdr <- df$ci.high
  }
  
  df <- df %>% mutate_at(c("OR", "ci.low", "ci.high", "ci.low.fdr", "ci.high.fdr"), ~round(., 2))
  
  # Paste CIs together
  df$ci95  <- paste("(", df$ci.low, "-", df$ci.high, ")", sep = "")
  df$ciFDR <- paste("(", df$ci.low.fdr, "-", df$ci.high.fdr, ")", sep = "")
  
  # Select columns and order
  output <- df %>% select(Compound = "display_name", "description", "OR",
                          "ci95", "ciFDR", "P.value", "FDR") %>% arrange(description)
  
}

all <- tidy.output(fits1)
pre <- tidy.output(fits2)
pre1 <- tidy.output(fits3)
post <- tidy.output(fits4)
preS1 <- tidy.output(fits7)
preS2 <- tidy.output(fits6)
pos <- tidy.output(fits8)
fast <- tidy.output(fits1a)
fast0 <- tidy.output(fits1b)
fast1 <- tidy.output(fits1c)

# Duration of OC use for reviewers' comments
all.oc <- tidy.output(fits1oc)
pre.oc <- tidy.output(fits2oc)
pos.oc <- tidy.output(fits4oc)

# Test for follow-up time
fu1 <- tidy.output(fits8)
fu2 <- tidy.output(fits9)
fu3 <- tidy.output(fits10)
fu4 <- tidy.output(fits1a)


# Retain only metabolite groups with at least one p-value < 0.05
tab <- bind_rows("All" = all, "Post" = post, "Pre" = pre, "Pre.eth" = pre1, .id = "Group") %>%
  arrange(Group, -OR) %>% 
  group_by(Compound) %>%
  filter(min(FDR) < 0.05) %>%
  select(Compound, description, everything()) 

# Copy and paste this into manuscript via Excel


# Revisions after reviewers' comments
# mean intensity by case-control status. Join unscaled metabolomics data and metadata (1572)
# Run data_pre up to line 81 only (unscaled)
alldat <- cbind(meta, ints0)
cmpds <- c("NAC", "Ethanol", "Histidine", "Glycerol", "Ornithine", "Leucine", "Albumin", "Glutamine",
           "Glutamate", "Pyruvate")

# Get overall means regardless of CC status (for scaling)
mn1 <- alldat %>% summarise_at(cmpds, mean, na.rm = T)
mn2 <- alldat %>% filter(MENOPAUSE == 1) %>% summarise_at(cmpds, mean, na.rm = T)
mn3 <- alldat %>% filter(MENOPAUSE == 0) %>% summarise_at(cmpds, mean, na.rm = T)

# Get means and SDs by case-control status
msd1 <- alldat %>% group_by(CT) %>% summarise_at(cmpds, mean, na.rm = T)
msd2 <- alldat %>% filter(MENOPAUSE == 1) %>% group_by(CT) %>% summarise_at(cmpds, mean, na.rm = T)
msd3 <- alldat %>% filter(MENOPAUSE == 0) %>% group_by(CT) %>% summarise_at(cmpds, mean, na.rm = T)

msd4 <- alldat %>% group_by(CT) %>% summarise_at(cmpds, sd, na.rm = T)
msd5 <- alldat %>% filter(MENOPAUSE == 1) %>% group_by(CT) %>% summarise_at(cmpds, sd, na.rm = T)
msd6 <- alldat %>% filter(MENOPAUSE == 0) %>% group_by(CT) %>% summarise_at(cmpds, sd, na.rm = T)

# Put data in correct order for means:
msd <- bind_rows("all" = msd1, "post" = msd2, "pre" = msd3, "pre1" = msd3, "all" = mn1, "post" = mn2,
                 "pre" = mn3, "pre1" = mn3, .id = "subgroup")

# Calculate scaled means
msdcalc <- msd %>%
  pivot_longer(-(CT:subgroup)) %>% pivot_wider(names_from = CT) %>% arrange(fct_inorder(name)) %>%
  mutate(ctrl.scale.grp = `0`/`NA`, case.scale.grp = `1`/`NA`)

# For SDs:
msd0 <- bind_rows("all" = msd4, "post" = msd5, "pre" = msd6, "pre1" = msd6, "all" = mn1, "post" = mn2,
                  "pre" = mn3, "pre1" = mn3, .id = "subgroup") %>% #replace_na(list(CT = 2)) %>%
  pivot_longer(-(CT:subgroup)) %>% pivot_wider(names_from = CT) %>% arrange(fct_inorder(name)) %>%
  mutate(sdT = `0`/`NA`, sdC = `1`/`NA`)


### For 2nd round of reviewers' comments: scale to all subjects, not group specific
msd.new <- bind_rows("all" = msd1, "post" = msd2, "pre" = msd3, "pre1" = msd3, "all" = mn1, "post" = mn1,
                     "pre" = mn1, "pre1" = mn1, .id = "subgroup") 

# Calculate scaled means
msdcalc1 <- msd.new %>%
  pivot_longer(-(CT:subgroup)) %>% pivot_wider(names_from = CT) %>% arrange(fct_inorder(name)) %>%
  mutate(ctrl.scale.all = `0`/`NA`, case.scale.all = `1`/`NA`)

# Scaled SDs
msd.new1 <- bind_rows("all" = msd4, "post" = msd5, "pre" = msd6, "pre1" = msd6, "all" = mn1, "post" = mn1,
                  "pre" = mn1, "pre1" = mn1, .id = "subgroup") 

msdcalc2 <- msd.new1 %>%
  pivot_longer(-(CT:subgroup)) %>% pivot_wider(names_from = CT) %>% arrange(fct_inorder(name)) %>%
  mutate(sdT = `0`/`NA`, sdC = `1`/`NA`)


# Old: forest plots (not used in manuscript)
t1 <- map_df(fits1, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% arrange(description)
t2 <- map_df(fits2, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% arrange(description)
t3 <- map_df(fits3, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta[-14, ]) %>% arrange(description)
t4 <- map_df(fits4, tidy) %>% filter(str_detect(term, "x")) %>% bind_cols(cmpd.meta) %>% arrange(description)

par(mfrow = c(1,2))
par(mar=c(5,4,1,2))
library(metafor)
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio per SD increase in conc", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t1$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

par(mar=c(5,4,1,2))
forest(t2$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio Q4 vs Q1", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t2$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

par(mfrow = c(1,2))
par(mar=c(5,4,1,2))
library(metafor)
forest(t3$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio per SD increase in conc", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t1$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

par(mar=c(5,4,1,2))
forest(t4$estimate, ci.lb = t2$conf.low, ci.ub = t2$conf.high, refline = 1,
       ylim = c(1, max(rowvec) + 3), xlab = "Odds ratio Q4 vs Q1", 
       transf = exp, rows = rowvec, efac = 0.5, pch = 18, cex = 0.8, psize = 1.5, 
       annosym = c("  (", " to ", ")"), slab = t2$display_name)

hh <- par("usr")
text(hh[1], max(rowvec) + 2, "Metabolite", pos = 4, cex = 0.8)
text(hh[2], max(rowvec) + 2, "OR [95% CI]", pos = 2, cex = 0.8)

# Funnel plots for metabolites
funnel(x = t2$estimate, sei = t2$std.error)
funnel(x = t2$estimate, sei = t2$std.error, yaxis = "vi")



# Output with stargazer (abandoned, easier to go through Excel)
library(stargazer)
stargazer(tab, summary = F, type = "html", out = "metabolite_table_selected_new.html")

# Manhattan plot
df <- bind_rows("All" = all, "Pre" = pre, "Post" = post, .id = "Group")

library(ggplot2)
ggplot(pre, aes(y = reorder(Compound, P.value), x = log10(P.value))) + 
  theme_minimal() + geom_point() + 
  #geom_vline(xintercept = -3, linetype = "dashed") +
  xlab("-log10(p-value)") +
  facet_grid(description ~ ., scales = "free_y", space = "free_y", switch= "x") +
  theme(axis.title.y = element_blank(), #axis.text.y = element_text(size=9),
        legend.position = c(0.25, 0.4),
        legend.box.background = element_rect(colour="grey")) #+
  #ggtitle("Metabolite associations with WCRF score (cal)")

# Plot data with Metafor
# Get vectors for row spacings using groups (may add compound classes later)
cmpds_ordered <- cmpd_meta %>% arrange(description) %>% mutate(row = 1:n() + (as.numeric(description)-1))
rowvec <- cmpds_ordered$row

# Investigation of Ethanol
boxplot(dat$Ethanol ~ dat$CT + dat$MENOPAUSE, varwidth = T, outline = F,
        names = c("Control, pre", "Case, pre", "Control, post", "Case, post"),
        col = "dodgerblue", ylab = "Plasma ethanol conc (scaled)")

boxplot(log(dat$ALCOHOL) ~ dat$CT + dat$MENOPAUSE, varwidth = T, outline = F,
        names = c("Control, pre", "Case, pre", "Control, post", "Case, post"),
        col = "hotpink", ylab = "Plasma ethanol conc (scaled)")

mod1 <- wilcox.test(ALCOHOL ~ CT, data = dat, subset = MENOPAUSE == 0)
mod2 <- wilcox.test(ALCOHOL ~ CT, data = dat, subset = MENOPAUSE == 1)


hist(log2(dat$Ethanol))
plot(log(dat$Ethanol), dat$ALCOHOL)
fit.e <- lm(Ethanol ~ ALCOHOL, data = dat)
summary(fit.e)

# For publication
library(ggsignif)
ggplot(dat, aes(x = as.factor(CT), y = Ethanol)) + geom_boxplot() + ylim(0, 2) +
  geom_signif(comparisons = c(0, 1), map_signif_level = T)

# Heatmap of differences---------------------

# Make a lookup df for menopausal status
match.men <- dat %>% select(MATCH, MENOPAUSE) %>% unique()

# Wide and long datasets for gplots and ggplot2
diff.wide <- dat %>% select(MATCH, CT,`3Hydroxybutyrate`:Succinate) %>% arrange(CT) %>%
  group_by(MATCH) %>% summarise_all(list(d = diff)) #%>% left_join(match.men, by = "MATCH")

diff.long <- gather(diff.wide, compound, difference, -MATCH)

# Vectors of pre and post-menopausal 
pre <- diff.wide$MENOPAUSE == 0
post <- diff.wide$MENOPAUSE == 1

hist(dat1$difference, breaks = 50)
which.max(dat1$difference)
which.min(dat1$difference)
dat1[10291, ]
dat1[34429, ]

# With heatmap2
mat1 <- as.matrix(diff.wide[pre, -1])
mat2 <- diff.wide[post, -1]

quantile.range <- quantile(as.matrix(diff.wide[, -(1:2)]), probs = seq(0, 1, 0.01))
palette.breaks <- seq(quantile.range["1%"], quantile.range["99%"], 0.1)
colpalette1 <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1)
colpalette2 <- colorRampPalette(c('#ef8a62','#f7f7f7','#67a9cf'))(length(palette.breaks) - 1)

library(gplots)
par(mar=c(10, 4.5, 0, 0.5))
heatmap.2(as.matrix(diff.wide[, -(1:2)]), trace = "none", col = colpalette2,
          breaks = palette.breaks, dendrogram = "col", margins = c(9, 3),
          offsetCol = 0, srtCol = 60, cexRow = 0.05)

# With ggplot
dat1 <- inner_join(diff.long, dat, by = "MATCH") %>% filter(CT == 0)

ggplot(dat1, aes(x= MATCH, y = compound, fill = difference)) + geom_tile() +
  facet_grid(. ~ MENOPAUSE, scales = "free", space = "free") +
  scale_fill_gradientn(colours = colpalette2)

# ---------------------------------------------------------------------------------------------------

# Conditional logistic regression to get odds ratios for lifestyle factors
# Adjusting for same co-variates as in original manuscript
library(survival)
# Test for associations between alcohol and BC
fit0 <- clogit(CT ~ scale(ALCOHOL) + BMI + SMK + DIABETE + RTH + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH), data = meta)
fit1 <- clogit(CT ~ scale(ALCOHOL) + BMI + SMK + DIABETE + RTH + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH), data = meta, subset = MENOPAUSE == 0) 
fit2 <- clogit(CT ~ scale(ALCOHOL) + BMI + SMK + DIABETE + RTH + DURTHSDIAG + CENTTIME + STOCKTIME + strata(MATCH), data = meta, subset = MENOPAUSE == 1) 

# fit0: 1.0784 [0.9704, 1.199]
# fit1: 0.9569 [0.77228, 1.186]
# fit2: 1.1035 [0.9739 , 1.250]

fit <- clogit(CT ~ scale(BMI) + SMK + DIABETE + #BP + 
                scale(RTH) + scale(ALCOHOL) + scale(DURTHSDIAG) + 
                scale(CENTTIME) + STOCKTIME + strata(MATCH), data = meta, subset = MENOPAUSE == 1) 
# output <- cbind(exp(coef(fit)), exp(confint(fit)))

library(broom)
t1 <- tidy(fit) %>% select(-(std.error:p.value))

library(metafor)
#dev.off()
par(mar=c(5,4,2,2))
forest(t1$estimate, ci.lb = t1$conf.low, ci.ub = t1$conf.high, refline = 1, 
       xlab = "Multivariable adjusted odds ratio",
       alim = c(0, 5),
       xlim = c(-5, 10),
       transf = exp, pch = 18, psize = 1.5, slab = t1$term,
       main = paste(fit$nevent, "case-control pairs"))
#xlim = c(-1, 3)
hh <- par("usr")
text(hh[1], nrow(t1) + 2, "Variable", pos = 4)
text(hh[2], nrow(t1) + 2, "OR [95% CI]", pos = 2)
# matching factors removed!

# Non-metabolite model lme4
library(lme4)
fit1 <- glmer(CT ~ BMI + SMK + DIABETE + (1|MATCH), data = meta1, family = "binomial")
summary(fit1)



# Manuscript figures----
# Formerly BC_manuscript_figs.R
library(survival)
library(broom)
library(ggrepel)
source("BC_prep_data.R")


# Figure 1A follow up times

ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
ints <- scale(ints0)
meta <- read.csv("metadata.csv", na = "9999")

# Stacked follow-up time histogram for paper
meta$menopause <- ifelse(meta$MENOPAUSE == 1, "Post-menopausal", "Pre-menopausal")

library(ggplot2)
plotA <- ggplot(meta, aes(x = DIAGSAMPLING, fill = menopause)) + 
  geom_histogram(position = "stack", colour = "black") + 
  theme_bw(base_size = 10) + 
  scale_fill_manual(values = c("grey60", "white"), guide = guide_legend(reverse = TRUE)) + 
  xlab("Time from blood collection to diagnosis (years)") + ylab("No. cases") +
  theme(legend.position = "top",
        #legend.position = c(0.8, 0.85), 
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.minor = element_blank()) #+
ggtitle("")


# Plot of alcohol intake versus ethanol concentration (also see Barplots_alcohol.R)
data <- bind_cols(meta, data.frame(ints))
data$alcQ4 <- cut_number(data$ALCOHOL, n=4, label = c("Q1", "Q2", "Q3", "Q4"))
#data$Subgroup <- as.factor(data$MENOPAUSE)
boxplot(Ethanol ~ alcQ4, data = data)
fit <- kruskal.test(Ethanol ~ alcQ4, data = data)

plotB <- ggplot(data, aes(x = alcQ4, y = Ethanol, fill = fct_rev(menopause))) + 
  geom_boxplot(outlier.shape = NA, position = position_dodge()) +
  scale_fill_manual(values = c("white", "grey")) +
  geom_point(position = position_jitterdodge(), size = 0.5) +
  annotate("text", label = "Kruskal-Wallis P = 0.70", x = 2.5, y = 1.5, size = 3) +
  ylim(-1.0, 1.1) + 
  theme_bw(base_size = 10) +
  xlab("Quartile of reported alcohol intake") +
  ylab("Plasma ethanol (scaled intensity)") +
  theme(panel.grid.minor = element_blank(),
        #legend.position = c(0.5, 0.1),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_blank(),
        legend.position = "top") #+
ggtitle("")


# Correlation heatmap for paper
library(corrplot)
colnames(ints)[1] <- "3-Hydroxybutyrate"
colnames(ints)[41] <- "N-acetyl glycoproteins"
cormat <- cor(ints, use = "pairwise.complete.obs")

library(ggcorrplot)
cordf <- as_tibble(cormat)
plotC <- ggcorrplot(cordf, hc.order = T, hc.method = "ward", legend.title = "Pearson\ncorrelation") + 
  theme_minimal() + scale_x_continuous(expand = c(0,0)) + ggtitle("") +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# For fasting participants only (reviewers revisions)
ints.fast <- ints[meta$FASTING == 1, ]
cormatF <- cor(ints.fast, use = "pairwise.complete.obs")
cordff <- as_tibble(cormatF)
plotD <-
  ggcorrplot(cordff, hc.order = T, hc.method = "ward", legend.title = "Pearson\ncorrelation") + 
  theme_minimal() + scale_x_continuous(expand = c(0,0)) + ggtitle("") +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Arrange plots: 2 steps
library(cowplot)
plotAB <- plot_grid(plotA, plotB, labels = c("A", "B"), ncol = 2, rel_widths = c(1,1))
plot_grid(plotAB, plotD, labels = c("A", "C"), nrow = 2, rel_heights = c(1,1.7))


# Corrplot and corrr: see BC_study_exploratory

# Figure 2: Univariate models and smile plots for manuscript
# New quicker way (not yet finished). Define models for subsets
library(survival)
#base <- CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + CENTTIME + STOCKTIME + strata(MATCH)
#mclr  <- function(x, dat) update(clogit(base, ~. + x, data = dat))

mclr0 <- function(x, dat) { 
  clogit(CT ~ x + BMI + SMK + DIABETE + RTH + ALCOHOL + DURTHSBMB + CENTTIME + STOCKTIME + 
           strata(MATCH), data = dat) 
}

mclr1 <- function(x, dat) { 
  clogit(CT ~ x + BMI + SMK + DIABETE + RTH + ALCOHOL + CENTTIME + STOCKTIME + strata(MATCH), 
         data = dat) 
}

mclr2 <- function(x, dat) { 
  clogit(CT ~ x + BMI + SMK + DIABETE + RTH + ALCOHOL + ETHANOL + CENTTIME + 
           STOCKTIME + strata(MATCH), data = dat) 
}

meta.eth <- cbind(meta[pre, ], ETHANOL = ints[pre, 14])

# Note: need to remove hormone treatment therapy variable DURTHSBMB
t1 <- apply(ints, 2, mclr0, dat = meta) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

t2 <- apply(ints[pre, ], 2, mclr1, dat = meta[pre, ]) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

t2a <- apply(ints[pre0, ], 2, mclr1, dat = meta[pre0, ]) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

t3 <- apply(ints[post, ], 2, mclr1, dat = meta[post, ]) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

t4 <- apply(ints[pre, -14], 2, mclr2, dat = meta.eth) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta[-14, ])

t6a <- apply(ints[agelo, ], 2, mclr0, dat = meta[agelo, ]) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta[-14, ])

t7a <- apply(ints[agehi, ], 2, mclr0, dat = meta[agehi, ]) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta[-14, ])

all <- bind_rows(all = t1, pre = t2, post = t3, eth = t4, .id = "analysis") %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

# All and post-menopausal for fasting subjects only
t8 <- apply(ints[fast, ], 2, mclr0, dat = meta[fast, ]) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

t9 <- apply(ints[fast1, ], 2, mclr1, dat = meta[fast1, ]) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

t10 <- apply(ints[fast0, ], 2, mclr1, dat = meta[fast0, ]) %>% map_df(tidy, exponentiate = T) %>%
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)



# Plot faceting by analysis
library(ggplot2)
ggplot(all, aes(x = estimate, y = -log10(p.value))) + geom_point(shape = 1) +
  facet_wrap(fct_inorder(analysis) ~ ., scales = "free_x") + theme_bw() +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  geom_hline(yintercept = -log10(0.05), size = 0.2, colour = "grey60")


# Create base plot to cut down code
base <- ggplot(t1, aes((estimate), log10(p.value))) + 
  #geom_point(shape = 21, fill = "lightgreen") + 
  geom_point(shape = 1) + 
  theme_bw(base_size = 10) +
  xlab("Odds ratio per SD increase concentration") + 
  ylab(expression(paste(italic(P), "-value"))) +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  geom_hline(yintercept = log10(0.05), size = 0.2, colour = "grey60") +
  theme(panel.grid.minor = element_blank() , panel.grid.major = element_blank(),
        axis.title.x = element_text(size = 9))

# All participants
p1 <- base %+% xlim(0.8, 1.2) +
  scale_y_reverse(limits = c(0, -2.5), breaks = c(0:-2), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t1[t1$p.value < 0.15, ]) + 
  #labs(title =  "All participants") + 
  labs(title = "Breast cancer-metabolite associations",
       subtitle =  "All participants") +
  annotate("text", x = 0.85, y=-1.37, size = 3, 
           label = expression(paste("Raw ", italic(P),"-threshold" )))

# Pre-menopausal
p2 <- base %+% t2 + xlim(0.3, 1.7) +
  scale_y_reverse(limits = c(0, -3.9), breaks = c(0:-3), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t2[t2$p.value < 0.04, ] ) +
  geom_hline(yintercept = log10(0.014), linetype = "dotted") +
  labs(title = "", subtitle =  "Pre-menopausal women")

# Post menopausal
p3 <-base %+% t3 + xlim(0.8, 1.2) + 
  scale_y_reverse(limits = c(0, -2.5), breaks = c(0:-2), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t3[t3$p.value < 0.12, ]) +
  #geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  labs(title = "", subtitle =  "Post-menopausal women") +
  annotate("text", x = 1.15, y=-1.37, size = 3, 
           label = expression(paste("Raw ", italic(P),"-threshold" )))

# Sensitivity analysis pre-menopausal adjusted ethanol
p4 <- base %+% t4 + xlim(c(0.3, 1.7)) +
  scale_y_reverse(limits = c(0, -3.9), breaks = c(0:-3), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t4[t4$p.value < 0.04, ] ) +
  geom_hline(yintercept = log10(0.001), linetype = "dotted") +
  labs(title = "", subtitle = "Pre-menopausal, adjusted for plasma ethanol") +
  annotate("text", x = 0.5, y=c(-1.4, -3.1), size = 3, #hjust = 0,
           label = c(expression(paste("Raw ", italic(P),"-threshold" )), 
                     expression(paste("FDR ", italic(P),"-threshold" ))))

library(cowplot)
plot_grid(p1, p3, p2, p4, labels = LETTERS[1:4], label_size = 12) 
#plot_grid(p1, p3, p2, p4) 
# Output at 8.18x8.41 inch


# Sensitivity analysis for reviewers fasting status Dec 2020

all <- bind_rows(`Pre-menopausal women` = t2, `Post-menopausal women` = t3, 
                 `Pre-menopausal women (fasting)` = t10, 
                 `Post-menopausal women (fasting)` = t9, .id = "analysis") %>% 
  filter(term == "x") %>% mutate(p.adj = p.adjust(p.value, method = "fdr"))

#nlist <- c("Pre-menopausal", "Post-menopausal", "Pre-meno. fasting", "Post-meno. fasting")

ggplot(all, aes(x = estimate, y = -log10(p.value))) + 
  geom_text_repel(aes(label = display_name), size = 3, segment.colour = "grey",
                  data = all[all$p.value < 0.06, ] ) +
  geom_point(shape = 1) +
  facet_wrap(fct_inorder(analysis)  ~ ., scales = "fixed") + 
  #facet_grid(. ~ fct_inorder(analysis), scales = "free_x") + 
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  geom_hline(yintercept = -log10(0.05), colour = "black", linetype = "dotted") +
  xlab("Odds ratio") + ylab(expression(paste(-log10, " (", italic(P), "-value)")))



