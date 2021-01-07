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
base <- ggplot(t1, aes((estimate), log10(p.value))) + geom_point(shape = 1) + 
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
  labs(title = "", subtitle =  "All participants") +
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
# Output at 8.18x8.41 inch

# Sensitivity analysis for reviewers fasting status Dec 2020

all <- bind_rows(Premenopausal = t2, Postmenopausal = t3, 
                 Premenopausal.fasting = t10, Postmenopausal.fasting = t9, .id = "analysis") %>% 
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
  geom_hline(yintercept = -log10(0.05), size = 0.2, colour = "grey60") +
  xlab("Odds ratio")











