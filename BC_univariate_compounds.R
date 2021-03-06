# BC risk models for metabolites and also for lifestyle variables alone
# Also heatmap of differences
source("BC_prep_data.R")
library(broom)
library(survival)

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
