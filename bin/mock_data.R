library(datasets)
library(ggplot2)
library(reshape2)
library(splinectomeR)

df_orig <- data.frame(matrix(ncol=3, nrow = 120))
colnames(df_orig) <- c('patient','time','response')
ids <- c(1,2,3,4,5,6,7,8,9,10)
timeseries <- c(0,2,4,6,8,10,12,14,16,18,20,22)
patient <- unlist(lapply(X = ids, FUN = function(xx) rep(xx, times=12)))
obs <- c(rnorm(n=10, mean = 100, sd = 2),
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

ggplot(df) + geom_smooth(aes(x=df$time, y=df$response), method = 'lm') +
  theme_classic() + labs(x='time',y='response')

z <- c(0,0,0.2,0.75,0.9,1,0.9,0.75,0.4,0.2,0,0)
shifter <- c(z,z,z,z,z,z,z,z,z,z)
df$resp_0.5z <- df$response + (df$response * (shifter/10))
df$resp_1z <- df$response + (df$response * (shifter/5))
df$resp_1.5z <- df$response + (df$response * (shifter/2))

ggplot(df) + geom_smooth(aes(x=df$time, y=df$response), se = F, method = 'loess') +
  geom_smooth(aes(x=df$time, y=df$resp_0.5z), se = F, method = 'loess') +
  geom_smooth(aes(x=df$time, y=df$resp_1z), se = F, method = 'loess') +
  geom_smooth(aes(x=df$time, y=df$resp_1.5z), se = F, method = 'loess') +
  theme_classic()

wilcox.test(df[df$Time <= 2, ]$weight, df[df$Time <= 2, ]$weight_1.5z, paired = F)


df2 <- df_orig
df3 <- df_orig
df4 <- df_orig
df5 <- df_orig
z <- c(0,0.2,.4,0.4,0.2,0.05,-.1,-.2,-.38,-0.3,-.1,0)
shifter2 <- c(z,z,z,z,z,z,z,z,z,z)
df2$set <- 'base'
df3$response <- df2$response + (df2$response * (shifter2/20))
df3$set <- 'base_z1'; df3$patient <- paste(df3$patient, 'z1', sep = '_')
df4$response <- df2$response + (df2$response * (shifter2/10))
df4$set <- 'base_z2'; df4$patient <- paste(df4$patient, 'z2', sep = '_')
df5$response <- df2$response + (df2$response * (shifter2/5))
df5$set <- 'base_z3'; df5$patient <- paste(df5$patient, 'z3', sep = '_')

df_bound <- rbind(df2, df3, df4, df5)

ggplot(df_bound) + geom_smooth(aes(x=time, y=response, color=set), se = T, method = 'loess') +
  theme_classic()

ggplot(df_bound) + geom_smooth(aes(x=time, y=response, color=set), se = T, method = 'lm') +
  theme_classic()

permu_result <- permuspliner(data = df_bound, xvar='time', yvar='response',
                             cases = 'patient', category = 'set',
                             perms = 99, groups = 'base,base_z1')
permuspliner.plot.permdistance(permu_result, xlabel = 'time')
permuspliner.plot.permsplines(permu_result, xvar = 'time', yvar = 'response')



