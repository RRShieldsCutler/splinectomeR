#### Simulated dataest to demonstrate the splinectomeR distance test ####

library(ggplot2)
library(reshape2)
library(splinectomeR)

colpal <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

df_orig <- data.frame(matrix(ncol=3, nrow = 120))
colnames(df_orig) <- c('patient','time','response')
ids <- c(1,2,3,4,5,6,7,8,9,10)
timeseries <- c(0,2,4,6,8,10,12,14,16,18,20,22)
patient <- unlist(lapply(X = ids, FUN = function(xx) rep(xx, times=12)))
obs <- c(rnorm(n=10, mean = 100, sd = 2),  # Make a distribution at each timepoint
         rnorm(n=10, mean = 100, sd = 4),
         rnorm(n=10, mean = 100, sd = 6),
         rnorm(n=10, mean = 100, sd = 8),
         rnorm(n=10, mean = 100, sd = 6),
         rnorm(n=10, mean = 100, sd = 4),
         rnorm(n=10, mean = 100, sd = 6),
         rnorm(n=10, mean = 100, sd = 8),
         rnorm(n=10, mean = 100, sd = 6),
         rnorm(n=10, mean = 100, sd = 4),
         rnorm(n=10, mean = 100, sd = 3),
         rnorm(n=10, mean = 100, sd = 2))
df_orig$patient <- patient
df_orig$time <-rep(timeseries, times = 10)
df_orig$response <- obs
df <- df_orig

#### A bump in the middle ####
dfm2 <- df_orig
dfm3 <- df_orig
dfm4 <- df_orig
dfm5 <- df_orig
z <- c(0,0.02,0.12,0.35,0.45,.55,0.45,0.38,0.2,0.1,0.01,0)  # Values by which to perturb the data
shifter <- c(z,z,z,z,z,z,z,z,z,z)
dfm2$Version <- 'base'
# Make perturbed data at several magnitudes by scaling z
dfm3$response <- dfm2$response + (dfm2$response * (shifter/20))
dfm3$Version <- 'shift_1x'; dfm3$patient <- paste(dfm3$patient, '1x', sep = '_')
dfm4$response <- dfm2$response + (dfm2$response * (shifter/10))
dfm4$Version <- 'shift_2x'; dfm4$patient <- paste(dfm4$patient, '2x', sep = '_')
dfm5$response <- dfm2$response + (dfm2$response * (shifter/5))
dfm5$Version <- 'shift_4x'; dfm5$patient <- paste(dfm5$patient, '4x', sep = '_')

dfm_bound <- rbind(dfm2, dfm3, dfm4, dfm5)

dfm_loess <- ggplot(dfm_bound) + geom_smooth(aes(x=time, y=response, color=Version), se = F, method = 'loess') +
  geom_point(aes(x=time, y=response, color=Version), alpha=0.5, size=1.2) +
  theme_classic() + theme(axis.text = element_text(color='black')) +
  scale_color_manual(values=colpal)
ggsave(dfm_loess, filename = '~/Box Sync/knights_box/splinectomer/doc/manuscript/figures/FigXA_dfmloess.png', width = 4, height = 3, dpi=600)

dfm_lm <- ggplot(dfm_bound) + geom_smooth(aes(x=time, y=response, color=Version), se = F, method = 'lm') +
  geom_point(aes(x=time, y=response, color=Version), alpha=0.5, size=1.2) +
  theme_classic() + theme(axis.text = element_text(color='black')) +
  scale_color_manual(values=colpal)
ggsave(dfm_lm, filename = '~/Box Sync/knights_box/splinectomer/doc/manuscript/figures/FigXA_dfmlm.png', width = 4, height = 3, dpi=600)

permu_result1 <- permuspliner(data = dfm_bound, xvar='time', yvar='response',
                             cases = 'patient', category = 'Version',
                             perms = 999, groups = 'base,shift_1x')
permu_result2 <- permuspliner(data = dfm_bound, xvar='time', yvar='response',
                             cases = 'patient', category = 'Version',
                             perms = 999, groups = 'base,shift_2x')
permu_result4 <- permuspliner(data = dfm_bound, xvar='time', yvar='response',
                             cases = 'patient', category = 'Version',
                             perms = 999, groups = 'base,shift_4x')

permuspliner.plot.permdistance(permu_result4, xlabel = 'time')
permuspliner.plot.permsplines(permu_result4, xvar = 'time', yvar = 'response')

slide_result <- sliding_spliner(data = df_bound, xvar='time', yvar='response',
                                cases = 'patient', category = 'Version',
                                groups = 'base,shift_2x')
sliding_spliner.plot.pvals(slide_result, xvar = 'time')  # Show there are small regions of significant difference

# If you average the first and last two points, is there a difference across time?
dfm2_w <- dcast(dfm2[, -4], patient ~ time)
dfm2_w$change <- ((dfm2_w$`22` + dfm2_w$`20`) / 2) - ((dfm2_w$`0` + dfm2_w$`2`) / 2)

dfm5_w <- dcast(dfm5[, -4], patient ~ time)
dfm5_w$change <- ((dfm5_w$`22` + dfm5_w$`20`) / 2) - ((dfm5_w$`0` + dfm5_w$`2`) / 2)

wilcox.test(dfm5_w$change, dfm2_w$change, paired = F)
# Answer: No


#### An early and late divergence ####
df2 <- df_orig
df3 <- df_orig
df4 <- df_orig
df5 <- df_orig
z <- c(0,0.2,.4,0.4,0.2,0.05,-.1,-.2,-.38,-0.3,-.1,0)  # Scale to perturb the data
shifter2 <- c(z,z,z,z,z,z,z,z,z,z)
df2$Version <- 'base'
# Make perturbed data at several magnitudes by scaling z
df3$response <- df2$response + (df2$response * (shifter2/20))
df3$Version <- 'shift_1x'; df3$patient <- paste(df3$patient, '1x', sep = '_')
df4$response <- df2$response + (df2$response * (shifter2/10))
df4$Version <- 'shift_2x'; df4$patient <- paste(df4$patient, '2x', sep = '_')
df5$response <- df2$response + (df2$response * (shifter2/5))
df5$Version <- 'shift_4x'; df5$patient <- paste(df5$patient, '4x', sep = '_')

df_bound <- rbind(df2, df3, df4, df5)

df_loess <- ggplot(df_bound) + geom_smooth(aes(x=time, y=response, color=Version), se = F, method = 'loess') +
  geom_point(aes(x=time, y=response, color=Version), alpha=0.5, size=1.2) +
  theme_classic() + theme(axis.text = element_text(color='black')) +
  scale_color_manual(values=colpal)
ggsave(df_loess, filename = '~/Box Sync/knights_box/splinectomer/doc/manuscript/figures/FigXA_df_loess.png', width = 4, height = 3, dpi=600)

df_lm <- ggplot(df_bound) + geom_smooth(aes(x=time, y=response, color=Version), se = F, method = 'lm') +
  geom_point(aes(x=time, y=response, color=Version), alpha=0.5, size=1.2) +
  theme_classic() + theme(axis.text = element_text(color='black')) +
  scale_color_manual(values=colpal)
ggsave(df_lm, filename = '~/Box Sync/knights_box/splinectomer/doc/manuscript/figures/FigXA_df_lm.png', width = 4, height = 3, dpi=600)

permu_result1 <- permuspliner(data = df_bound, xvar='time', yvar='response',
                             cases = 'patient', category = 'Version',
                             perms = 999, groups = 'base,shift_1x')
permu_result2 <- permuspliner(data = df_bound, xvar='time', yvar='response',
                             cases = 'patient', category = 'Version',
                             perms = 999, groups = 'base,shift_2x')
permu_result4 <- permuspliner(data = df_bound, xvar='time', yvar='response',
                             cases = 'patient', category = 'Version',
                             perms = 999, groups = 'base,shift_4x')
# Run with 99 perms for plotting speed
permu_result_plot <- permuspliner(data = df_bound, xvar='time', yvar='response',
                              cases = 'patient', category = 'Version',
                              perms = 99, groups = 'base,shift_4x')
permuspliner.plot.permdistance(permu_result_plot, xlabel = 'time')
permuspliner.plot.permsplines(permu_result_plot, xvar = 'time', yvar = 'response')

# If you average the first and last two points, is there a difference across time?
df2_w <- dcast(df2[, -4], patient ~ time)
df2_w$change <- ((df2_w$`22` + df2_w$`20`) / 2) - ((df2_w$`0` + df2_w$`2`) / 2)

df5_w <- dcast(df5[, -4], patient ~ time)
df5_w$change <- ((df5_w$`22` + df5_w$`20`) / 2) - ((df5_w$`0` + df5_w$`2`) / 2)

wilcox.test(df5_w$change, df2_w$change, paired = F)
# Answer here, also no.
