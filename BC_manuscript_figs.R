library(survival)
library(broom)
library(ggrepel)
source("BC_prep_data.R")


# Figure 1A follow up times

ints0 <- read_tsv("1510_XMetaboliteE3N_cpmg_unscaled.txt")
ints <- scale(ints0)
meta <- read.csv("metadata.csv", na = "9999")

# Follow-up time plot for paper
meta$menopause <- ifelse(meta$MENOPAUSE == 1, "Post-menopausal", "Pre-menopausal")

# Facetted (old)
library(ggplot2)
plot1 <- ggplot(meta, aes(x = DIAGSAMPLING)) + 
  geom_histogram(colour = "black", fill = "grey80")+ theme_minimal() + 
  scale_fill_manual(values = c("white", "grey80")) + 
  xlab("Time from blood collection to diagnosis (years)") + ylab("No. cases") +
  theme(legend.position = c(0.8, 0.8)) +
  facet_grid(.~menopause, scales = "free_y") + ggtitle("")

# Stacked (new)
plotA <- ggplot(meta, aes(x = DIAGSAMPLING, fill = menopause)) + 
  geom_histogram(position = "stack", colour = "black") + 
  theme_bw(base_size = 10) + 
  scale_fill_manual(values = c("grey60", "white")) + 
  xlab("Time from blood collection to diagnosis (years)") + ylab("No. cases") +
  theme(legend.position = c(0.8, 0.85),
        legend.title = element_blank()) +
  ggtitle("")


# Plot of alcohol intake versus ethanol concentration



# Correlation heatmap for paper
library(corrplot)
colnames(ints)[1] <- "3-Hydroxybutyrate"
colnames(ints)[41] <- "N-acetyl glycoproteins"
cormat <- cor(ints, use = "pairwise.complete.obs")

library(ggcorrplot)
cordf <- as_tibble(cormat)
plotC <- ggcorrplot(cordf, hc.order = T, hc.method = "ward", legend.title = "Scale") + theme_minimal() +
  scale_x_continuous(expand = c(0,0)) + ggtitle("") +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Arrange plots: 2 steps
library(cowplot)
plotAB <- plot_grid(plotA, plotB, labels = c("A", "B"), ncol = 2, rel_widths = c(1,1))
plot_grid(plotAB, plotC, labels = c("A", "C"), nrow = 2, rel_heights = c(1,1.7))



# Corrplot and corrr: see BC_study_exploratory





# Figure 2: Univariate models and smile plots for manuscript

fits0 <- apply(ints, 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta))

t1 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

# Create base plot to cut down code
base <- ggplot(t1, aes(exp(estimate), log10(p.value))) + geom_point(shape = 1) + 
  theme_bw(base_size = 10) +
  xlab("Odds ratio per SD increase concentration") + ylab(expression(italic(P)~ "-value")) +
  geom_vline(xintercept = 1, size = 0.2, colour = "grey60") + 
  geom_hline(yintercept = log10(0.05), size = 0.2, colour = "grey60") +
  theme(panel.grid.minor = element_blank() , 
        panel.grid.major = element_blank())

p1 <- 
  base %+% xlim(0.8, 1.2) +
  scale_y_reverse(limits = c(0, -2.5), breaks = c(0:-2), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t1[t1$p.value < 0.15, ]) + 
  labs(title = "A", subtitle = "All participants")

# Pre-menopausal
meta1 <- meta[pre, ]
# Note: need to remove hormone treatment therapy variable DURTHSBMB
fits1 <- apply(ints[pre, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta[pre, ]))

t2 <- map_df(fits1, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p2 <- 
  base %+% t2 + xlim(0.3, 1.7) +
  scale_y_reverse(limits = c(0, -3.9), breaks = c(0:-3), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t2[t2$p.value < 0.04, ] ) +
  geom_hline(yintercept = log10(0.014), linetype = "dotted") +
  labs(title = "B", subtitle = "Pre-menopausal")


# Post menopausal
meta2 <- meta[post, ]
fits2 <- apply(ints[post, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta2))

t3 <- map_df(fits2, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p3 <-
  base %+% t3 + xlim(0.8, 1.2) + 
  scale_y_reverse(limits = c(0, -2.5), breaks = c(0:-2), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t3[t3$p.value < 0.12, ]) +
  #geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  labs(title = "C", subtitle = "Post-menopausal")


# Sensitivity analysis adjusted ethanol
# All subjects
eth <- ints[pre, 14]
meta.eth <- cbind(meta1, ETHANOL = eth)

fits1e <- apply(ints[pre, -14], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL +
          ETHANOL + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta.eth))

t4 <- map_df(fits1e, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta[-14, ])


p4 <- 
  base %+% t4 + xlim(c(0.3, 1.7)) +
  scale_y_reverse(limits = c(0, -3.9), breaks = c(0:-3), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t4[t4$p.value < 0.04, ] ) +
  geom_hline(yintercept = log10(0.001), linetype = "dotted") +
  labs(title = "D", subtitle = "Pre-menopausal, adjusted for plasma ethanol")
  
library(cowplot)
plot_grid(p1, p3, p2, p4, #labels = c('A', 'B', 'C', 'D'), 
          label_size = 12)  


# Pre-menopausal, exclude first 2 years of follow up
meta3 <- meta[pre0, ]
# Note: need to remove hormone treatment therapy variable DURTHSBMB
fits1 <- apply(ints[pre0, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
    CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta3))

t5 <- map_df(fits1, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p5 <- 
  base %+% t5 + xlim(0.3, 1.7) + 
  scale_y_reverse(breaks = c(-3, -2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t5[t5$p.value < 0.04, ] ) +
  geom_hline(yintercept = c(log10(0.05), log10(0.014)), linetype = c("dashed", "dotted")) +
  ggtitle("Pre-menopausal, follow-up < 2 y excluded")

# Stratification by age
table(Age = meta$Age1, Menopause = meta$MENOPAUSE)

# Under 55
meta4 <- meta[agelo, ]
fits0 <- apply(ints[agelo, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
  DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta4))

t6 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p6 <- 
  base %+% t6 + #xlim(0.8, 1.2) +
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t6[t6$p.value < 0.3, ]) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") + ggtitle("Age < 55")

# Over 55
meta5 <- meta[agehi, ]
fits0 <- apply(ints[agehi, ], 2, function(x) clogit(CT ~ BMI + SMK + DIABETE + RTH + ALCOHOL + 
   DURTHSBMB + CENTTIME + STOCKTIME + strata(MATCH) + x, data = meta5))

t7 <- map_df(fits0, tidy) %>% filter(str_detect(term, "x")) %>% 
  mutate(p.valueFDR = p.adjust(p.value, method = "fdr")) %>% bind_cols(cmpd.meta)

p7 <- 
  base %+% t7 + #xlim(0.8, 1.2) +
  scale_y_reverse(breaks = c(-2, -1, 0), labels = function(x) 10^x) +
  geom_text_repel(aes(label = display_name), size = 3, data = t7[t7$p.value < 0.3, ]) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") + ggtitle("Age > 55")













