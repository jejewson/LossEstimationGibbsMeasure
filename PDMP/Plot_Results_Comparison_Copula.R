

library(ggplot2)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(GGally)
library(tidyr)


load(file = "ResultsZZ_copula_n1000_Gold.RData")
load(file = "ResultsPM_copula_block_n1000.RData")
load(file = "ResultsPM_copula_block_n1000_m1000.RData")
#load(file = "ResultsPM_copula_n1000.RData")
load(file = "ResultsZZ_copula_n1000.RData")


samples_pmcmc_block_m5 = 2/(1 + exp(-samples_pmcmc_block_m5)) - 1
samples_pmcmc_block_m10 = 2/(1 + exp(-samples_pmcmc_block_m10)) - 1
samples_pmcmc_block_m20 = 2/(1 + exp(-samples_pmcmc_block_m20)) - 1
samples_pmcmc_block_m50 = 2/(1 + exp(-samples_pmcmc_block_m50)) - 1
samples_pmcmc_block_m100 = 2/(1 + exp(-samples_pmcmc_block_m100)) - 1
samples_pmcmc_block_m1000 = 2/(1 + exp(-samples_pmcmc_block_m1000)) - 1

sample_copula_IntractLik_m2_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m2_rcpp)) - 1
sample_copula_IntractLik_m5_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m5_rcpp)) - 1
sample_copula_IntractLik_m10_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m10_rcpp)) - 1
sample_copula_IntractLik_m20_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m20_rcpp)) - 1
sample_copula_IntractLik_m50_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m50_rcpp)) - 1

sample_copula_IntractLik_gold_rcpp = 2/(1 + exp(-sample_copula_IntractLik_gold_rcpp)) - 1

gold_mean = colMeans(sample_copula_IntractLik_gold_rcpp)
gold_sd = colSds(sample_copula_IntractLik_gold_rcpp)



sims_pm_block_m5 = (1:length(samples_pmcmc_block_m5[,1]))*5
time_pm_block_m5 = seq(as.numeric(time.taken.block_m5)/length(samples_pmcmc_block_m5),as.numeric(time.taken.block_m5),by=as.numeric(time.taken.block_m5)/length(samples_pmcmc_block_m5))

running_mean_pm_block_m5 <- apply(samples_pmcmc_block_m5, 2, function(x) cumsum(x) / seq_along(x))


running_var_pm_block_m5 <- apply(samples_pmcmc_block_m5, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m5 <- sqrt(running_var_pm_block_m5)

sims_pm_block_m10 = (1:length(samples_pmcmc_block_m10[,1]))*10
time_pm_block_m10 = seq(as.numeric(time.taken.block_m10)/length(samples_pmcmc_block_m10),as.numeric(time.taken.block_m10),by=as.numeric(time.taken.block_m10)/length(samples_pmcmc_block_m10))

running_mean_pm_block_m10 <- apply(samples_pmcmc_block_m10, 2, function(x) cumsum(x) / seq_along(x))

running_var_pm_block_m10 <- apply(samples_pmcmc_block_m10, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m10 <- sqrt(running_var_pm_block_m10)

sims_pm_block_m20 = (1:length(samples_pmcmc_block_m20[,1]))*20
time_pm_block_m20 = seq(as.numeric(time.taken.block_m20)/length(samples_pmcmc_block_m20),as.numeric(time.taken.block_m20),by=as.numeric(time.taken.block_m20)/length(samples_pmcmc_block_m20))

running_mean_pm_block_m20 <- apply(samples_pmcmc_block_m20, 2, function(x) cumsum(x) / seq_along(x))

running_var_pm_block_m20 <- apply(samples_pmcmc_block_m20, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m20 <- sqrt(running_var_pm_block_m20)

sims_pm_block_m50 = (1:length(samples_pmcmc_block_m50[,1]))*50
time_pm_block_m50 = seq(as.numeric(time.taken.block_m50)/length(samples_pmcmc_block_m50),as.numeric(time.taken.block_m50),by=as.numeric(time.taken.block_m50)/length(samples_pmcmc_block_m50))

running_mean_pm_block_m50 <- apply(samples_pmcmc_block_m50, 2, function(x) cumsum(x) / seq_along(x))

running_var_pm_block_m50 <- apply(samples_pmcmc_block_m50, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m50 <- sqrt(running_var_pm_block_m50)

sims_pm_block_m100 = (1:length(samples_pmcmc_block_m100[,1]))*100
time_pm_block_m100 = seq(as.numeric(time.taken.block_m100)/length(samples_pmcmc_block_m100),as.numeric(time.taken.block_m100),by=as.numeric(time.taken.block_m100)/length(samples_pmcmc_block_m100))

running_mean_pm_block_m100 <- apply(samples_pmcmc_block_m100, 2, function(x) cumsum(x) / seq_along(x))


running_var_pm_block_m100 <- apply(samples_pmcmc_block_m100, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m100 <- sqrt(running_var_pm_block_m100)


sims_pm_block_m1000 = (1:length(samples_pmcmc_block_m1000[,1]))*1000
time_pm_block_m1000 = seq(60*as.numeric(time.taken.block_m1000)/length(samples_pmcmc_block_m1000),60*as.numeric(time.taken.block_m1000),by=60*as.numeric(time.taken.block_m1000)/length(samples_pmcmc_block_m1000))

running_mean_pm_block_m1000 <- apply(samples_pmcmc_block_m1000, 2, function(x) cumsum(x) / seq_along(x))


running_var_pm_block_m1000 <- apply(samples_pmcmc_block_m1000, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m1000 <- sqrt(running_var_pm_block_m1000)




num_sims_zz_m2 = length(skeleton_copula_IntractLik_m2_rcpp$skeleton_T)*1/skeleton_copula_IntractLik_m2_rcpp$mean_alpha*2
sims_zz_m2 = seq(1, num_sims_zz_m2, by = num_sims_zz_m2/length(sample_copula_IntractLik_m2_rcpp[,1]))
time_zz_m2 = seq(as.numeric(timer_MMD_copula_m2_rcpp)/length(sample_copula_IntractLik_m2_rcpp),as.numeric(timer_MMD_copula_m2_rcpp),by=as.numeric(timer_MMD_copula_m2_rcpp)/length(sample_copula_IntractLik_m2_rcpp))


running_mean_zz_m2 <- apply(sample_copula_IntractLik_m2_rcpp, 2, function(x) cumsum(x) / seq_along(x))

running_var_zz_m2 <- apply(sample_copula_IntractLik_m2_rcpp, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m2 <- sqrt(running_var_zz_m2)




num_sims_zz_m5 = length(skeleton_copula_IntractLik_m5_rcpp$skeleton_T)*1/skeleton_copula_IntractLik_m5_rcpp$mean_alpha*2
sims_zz_m5 = seq(1, num_sims_zz_m5, by = num_sims_zz_m5/length(sample_copula_IntractLik_m5_rcpp[,1]))
time_zz_m5 = seq(as.numeric(timer_MMD_copula_m5_rcpp)/length(sample_copula_IntractLik_m5_rcpp),as.numeric(timer_MMD_copula_m5_rcpp),by=as.numeric(timer_MMD_copula_m5_rcpp)/length(sample_copula_IntractLik_m5_rcpp))

running_mean_zz_m5 <- apply(sample_copula_IntractLik_m5_rcpp, 2, function(x) cumsum(x) / seq_along(x))

running_var_zz_m5 <- apply(sample_copula_IntractLik_m5_rcpp, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m5 <- sqrt(running_var_zz_m5)


num_sims_zz_m10 = length(skeleton_copula_IntractLik_m10_rcpp$skeleton_T)*1/skeleton_copula_IntractLik_m10_rcpp$mean_alpha*2
sims_zz_m10 = seq(1, num_sims_zz_m10, by = num_sims_zz_m10/length(sample_copula_IntractLik_m10_rcpp[,1]))
time_zz_m10 = seq(as.numeric(timer_MMD_copula_m10_rcpp)/length(sample_copula_IntractLik_m10_rcpp),as.numeric(timer_MMD_copula_m10_rcpp),by=as.numeric(timer_MMD_copula_m10_rcpp)/length(sample_copula_IntractLik_m10_rcpp))

running_mean_zz_m10 <- apply(sample_copula_IntractLik_m10_rcpp, 2, function(x) cumsum(x) / seq_along(x))

running_var_zz_m10 <- apply(sample_copula_IntractLik_m10_rcpp, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m10 <- sqrt(running_var_zz_m10)


num_sims_zz_m20 = length(skeleton_copula_IntractLik_m20_rcpp$skeleton_T)*1/skeleton_copula_IntractLik_m20_rcpp$mean_alpha*2
sims_zz_m20 = seq(1, num_sims_zz_m20, by = num_sims_zz_m20/length(sample_copula_IntractLik_m20_rcpp[,1]))
time_zz_m20 = seq(60*as.numeric(timer_MMD_copula_m20_rcpp)/length(sample_copula_IntractLik_m20_rcpp),60*as.numeric(timer_MMD_copula_m20_rcpp),by=60*as.numeric(timer_MMD_copula_m20_rcpp)/length(sample_copula_IntractLik_m20_rcpp))

running_mean_zz_m20 <- apply(sample_copula_IntractLik_m20_rcpp, 2, function(x) cumsum(x) / seq_along(x))

running_var_zz_m20 <- apply(sample_copula_IntractLik_m20_rcpp, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m20 <- sqrt(running_var_zz_m20)


num_sims_zz_m50 = length(skeleton_copula_IntractLik_m50_rcpp$skeleton_T)*1/skeleton_copula_IntractLik_m50_rcpp$mean_alpha*2
sims_zz_m50 = seq(1, num_sims_zz_m50, by = num_sims_zz_m50/length(sample_copula_IntractLik_m50_rcpp[,1]))
time_zz_m50 = seq(as.numeric(timer_MMD_copula_m50_rcpp)*60/length(sample_copula_IntractLik_m50_rcpp),as.numeric(timer_MMD_copula_m50_rcpp)*60,by=as.numeric(timer_MMD_copula_m50_rcpp)*60/length(sample_copula_IntractLik_m50_rcpp))

running_mean_zz_m50 <- apply(sample_copula_IntractLik_m50_rcpp, 2, function(x) cumsum(x) / seq_along(x))

running_var_zz_m50 <- apply(sample_copula_IntractLik_m50_rcpp, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m50 <- sqrt(running_var_zz_m50)







############### PM plots for the first moment vs number of model simulations

param = 1 # change parameter number as desired


# Create a dataframe for ggplot
data <- data.frame(
  x = c(sims_pm_block_m20, sims_pm_block_m50, sims_pm_block_m100, sims_pm_block_m1000),
  y = c(running_mean_pm_block_m20[,param] - gold_mean[param], 
        running_mean_pm_block_m50[,param] - gold_mean[param], 
        running_mean_pm_block_m100[,param] - gold_mean[param], 
        running_mean_pm_block_m1000[,param] - gold_mean[param]),
  group = rep(c("PM 20", "PM 50", "PM 100", "PM 1000"), 
              times = c(length(sims_pm_block_m20), length(sims_pm_block_m50), length(sims_pm_block_m100), 
                        length(sims_pm_block_m1000)))
)

# Ensure the factor levels are in the correct order for legend appearance
data$group <- factor(data$group, levels = c("PM 20", "PM 50", "PM 100", "PM 1000"))

# Define colors and line types
color_map <- c("PM 20" = "#FF6F61", "PM 50" = "#377EB8", "PM 100" = "#32CD32", "PM 1000" = "#6A3D9A")
line_map <- c("PM 20" = "solid", "PM 50" = "dashed", "PM 100" = "dotted", "PM 1000" = "dotdash")

# Create the plot
p = ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(y = "Estimate - Gold") +
  ylim(-0.10, 0.10) +  # Restrict y-axis limits
  xlim(0, 1e7) +  # Restrict x-axis limits
  theme_minimal(base_size = 20) +  # Clean theme
  theme(
    axis.title.x=element_blank(),
    legend.title = element_blank(),  # Remove legend title
    legend.position = c(1,1),  # Move legend to bottom
    legend.key.width = unit(2, "cm"),
    legend.justification = c(1, 1)
  )

x11()
p






############### PM plots for the second moment vs number of model simulations

param = 1 # change parameter number as desired


# Create a dataframe for ggplot
data <- data.frame(
  x = c(sims_pm_block_m20, sims_pm_block_m50, sims_pm_block_m100, sims_pm_block_m1000),
  y = c(running_sd_pm_block_m20[,param] - gold_sd[param], 
        running_sd_pm_block_m50[,param] - gold_sd[param], 
        running_sd_pm_block_m100[,param] - gold_sd[param], 
        running_sd_pm_block_m1000[,param] - gold_sd[param]),
  group = rep(c("PM 20", "PM 50", "PM 100", "PM 1000"), 
              times = c(length(sims_pm_block_m20), length(sims_pm_block_m50), 
                        length(sims_pm_block_m100), length(sims_pm_block_m1000)))
)

# Ensure the factor levels are in the correct order for legend appearance
data$group <- factor(data$group, levels = c("PM 20", "PM 50", "PM 100", "PM 1000"))

# Define colors and line types
color_map <- c("PM 20" = "#FF6F61", "PM 50" = "#377EB8", "PM 100" = "#32CD32", "PM 1000" = "#6A3D9A")
line_map <- c("PM 20" = "solid", "PM 50" = "dashed", "PM 100" = "dotted", "PM 1000" = "dotdash")

# Create the plot
p = ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Number of Model Simulations", y = "Estimate - Gold") +
  ylim(-0.05, 0.15) +  # Restrict y-axis limits
  xlim(0, 1e7) +  # Restrict x-axis limits
  theme_minimal(base_size = 20) +  # Clean theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = c(1,1),  # Move legend to bottom
    legend.key.width = unit(2, "cm"),
    legend.justification = c(1, 1)
  )

x11()
p





############### ZZ plots for the first moment vs number of model simulations


# Create a dataframe for ggplot
data <- data.frame(
  x = c(sims_zz_m2, sims_zz_m5, sims_zz_m20, sims_zz_m50),
  y = c(running_mean_zz_m2[,param] - gold_mean[param], 
        running_mean_zz_m5[,param] - gold_mean[param], 
        running_mean_zz_m20[,param] - gold_mean[param], 
        running_mean_zz_m50[,param] - gold_mean[param]),
  group = rep(c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ 50"), 
              times = c(length(sims_zz_m2), length(sims_zz_m5), 
                        length(sims_zz_m20), length(sims_zz_m50)))
)

# Explicitly order the factor levels
data$group <- factor(data$group, levels = c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ 50"))

# Define colors and line types
color_map <- c("ZZ 2" = "#FF6F61", "ZZ 5" = "#377EB8", "ZZ 20" = "#32CD32", "ZZ 50" = "#6A3D9A")
line_map <- c("ZZ 2" = "solid", "ZZ 5" = "dashed", "ZZ 20" = "dotted", "ZZ 50" = "dotdash")

# Create the plot
p <- ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  #labs(x = "Number of Model Simulations", y = "Estimate - Gold") +
  ylim(-0.1, 0.1) +  # Restrict y-axis limits
  xlim(0, 1e7) +  # Restrict x-axis limits
  theme_minimal(base_size = 20) +  # Clean theme
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.title = element_blank(),  # Remove legend title
    legend.position = c(1,1),  # Move legend to bottom
    legend.key.width = unit(2, "cm"),
    legend.justification = c(1, 1)
  )

x11()
p




############### ZZ plots for the second moment vs number of model simulations


# Create a dataframe for ggplot
data <- data.frame(
  x = c(sims_zz_m2, sims_zz_m5, sims_zz_m20, sims_zz_m50),
  y = c(running_sd_zz_m2[,param] - gold_sd[param], 
        running_sd_zz_m5[,param] - gold_sd[param], 
        running_sd_zz_m20[,param] - gold_sd[param], 
        running_sd_zz_m50[,param] - gold_sd[param]),
  group = rep(c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ 50"), 
              times = c(length(sims_zz_m2), length(sims_zz_m5), 
                        length(sims_zz_m20), length(sims_zz_m50)))
)

# Explicitly order the factor levels
data$group <- factor(data$group, levels = c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ 50"))

# Define colors and line types
color_map <- c("ZZ 2" = "#FF6F61", "ZZ 5" = "#377EB8", "ZZ 20" = "#32CD32", "ZZ 50" = "#6A3D9A")
line_map <- c("ZZ 2" = "solid", "ZZ 5" = "dashed", "ZZ 20" = "dotted", "ZZ 50" = "dotdash")

# Create the plot
p <- ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Number of Model Simulations") +
  ylim(-0.05, 0.15) +  # Restrict y-axis limits
  xlim(0, 1e7) +  # Restrict x-axis limits
  theme_minimal(base_size = 20) +  # Clean theme
  theme(
    axis.title.y=element_blank(),
    legend.title = element_blank(),  # Remove legend title
    legend.position = c(1,1),  # Move legend to bottom
    legend.key.width = unit(2, "cm"),
    legend.justification = c(1, 1)
  )

x11()
p






############### PM plots for the first moment vs time


# Create a dataframe for ggplot
data <- data.frame(
  x = c(time_pm_block_m20, time_pm_block_m50, time_pm_block_m100, time_pm_block_m1000),
  y = c(running_mean_pm_block_m20[,param] - gold_mean[param], 
        running_mean_pm_block_m50[,param] - gold_mean[param], 
        running_mean_pm_block_m100[,param] - gold_mean[param], 
        running_mean_pm_block_m1000[,param] - gold_mean[param]),
  group = rep(c("PM 20", "PM 50", "PM 100", "PM 1000"), 
              times = c(length(time_pm_block_m20), length(time_pm_block_m50), 
                        length(time_pm_block_m100), length(time_pm_block_m1000)))
)

# Ensure the factor levels are in the correct order for legend appearance
data$group <- factor(data$group, levels = c("PM 20", "PM 50", "PM 100", "PM 1000"))

# Define colors and line types
color_map <- c("PM 20" = "#FF6F61", "PM 50" = "#377EB8", "PM 100" = "#32CD32", "PM 1000" = "#6A3D9A")
line_map <- c("PM 20" = "solid", "PM 50" = "dashed", "PM 100" = "dotted", "PM 1000" = "dotdash")

# Create the plot
p = ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Time (minutes)", y = "Estimate - Gold") +
  ylim(-0.10, 0.10) +  # Restrict y-axis limits
  xlim(0, 20) +  # Restrict x-axis limits
  theme_minimal(base_size = 14) +  # Clean theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = "bottom",  # Move legend to bottom
    legend.key.width = unit(2, "cm")  # Increase legend key size
  )

x11()
p





############### PM plots for the second moment vs time


# Create a dataframe for ggplot
data <- data.frame(
  x = c(time_pm_block_m20, time_pm_block_m50, time_pm_block_m100, time_pm_block_m1000),
  y = c(running_sd_pm_block_m20[,param] - gold_sd[param], 
        running_sd_pm_block_m50[,param] - gold_sd[param], 
        running_sd_pm_block_m100[,param] - gold_sd[param], 
        running_sd_pm_block_m1000[,param] - gold_sd[param]),
  group = rep(c("PM 20", "PM 50", "PM 100", "PM 1000"), 
              times = c(length(time_pm_block_m20), length(time_pm_block_m50), 
                        length(time_pm_block_m100), length(time_pm_block_m1000)))
)

# Ensure the factor levels are in the correct order for legend appearance
data$group <- factor(data$group, levels = c("PM 20", "PM 50", "PM 100", "PM 1000"))

# Define colors and line types
color_map <- c("PM 20" = "#FF6F61", "PM 50" = "#377EB8", "PM 100" = "#32CD32", "PM 1000" = "#6A3D9A")
line_map <- c("PM 20" = "solid", "PM 50" = "dashed", "PM 100" = "dotted", "PM 1000" = "dotdash")

# Create the plot
p = ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Time (minutes)", y = "Estimate - Gold") +
  ylim(-0.05, 0.15) +  # Restrict y-axis limits
  xlim(0, 20) +  # Restrict x-axis limits
  theme_minimal(base_size = 14) +  # Clean theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = "bottom",  # Move legend to bottom
    legend.key.width = unit(2, "cm")  # Increase legend key size
  )

x11()
p




############### ZZ plots for the first moment vs time



# Create a dataframe for ggplot
data <- data.frame(
  x = c(time_zz_m2, time_zz_m5, time_zz_m20, time_zz_m50),
  y = c(running_mean_zz_m2[,param] - gold_mean[param], 
        running_mean_zz_m5[,param] - gold_mean[param], 
        running_mean_zz_m20[,param] - gold_mean[param], 
        running_mean_zz_m50[,param] - gold_mean[param]),
  group = rep(c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ 50"), 
              times = c(length(time_zz_m2), length(time_zz_m5), 
                        length(time_zz_m20), length(time_zz_m50)))
)

# Explicitly order the factor levels
data$group <- factor(data$group, levels = c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ 50"))

# Define colors and line types
color_map <- c("ZZ 2" = "#FF6F61", "ZZ 5" = "#377EB8", "ZZ 20" = "#32CD32", "ZZ 50" = "#6A3D9A")
line_map <- c("ZZ 2" = "solid", "ZZ 5" = "dashed", "ZZ 20" = "dotted", "ZZ 50" = "dotdash")

# Create the plot
p <- ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Time (minutes)", y = "Estimate - Gold") +
  ylim(-0.1, 0.1) +  # Restrict y-axis limits
  xlim(0, 20) +  # Restrict x-axis limits
  theme_minimal(base_size = 14) +  # Clean theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = "bottom",  # Move legend to bottom
    legend.key.width = unit(2, "cm")  # Increase legend key size
  )

x11()
p









############### ZZ plots for the second moment vs time



# Create a dataframe for ggplot
data <- data.frame(
  x = c(time_zz_m2, time_zz_m5, time_zz_m20, time_zz_m50),
  y = c(running_sd_zz_m2[,param] - gold_sd[param], 
        running_sd_zz_m5[,param] - gold_sd[param], 
        running_sd_zz_m20[,param] - gold_sd[param], 
        running_sd_zz_m50[,param] - gold_sd[param]),
  group = rep(c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ 50"), 
              times = c(length(time_zz_m2), length(time_zz_m5), 
                        length(time_zz_m20), length(time_zz_m50)))
)

# Explicitly order the factor levels
data$group <- factor(data$group, levels = c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ 50"))

# Define colors and line types
color_map <- c("ZZ 2" = "#FF6F61", "ZZ 5" = "#377EB8", "ZZ 20" = "#32CD32", "ZZ 50" = "#6A3D9A")
line_map <- c("ZZ 2" = "solid", "ZZ 5" = "dashed", "ZZ 20" = "dotted", "ZZ 50" = "dotdash")

# Create the plot
p <- ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Time (minutes)", y = "Estimate - Gold") +
  ylim(-0.05, 0.15) +  # Restrict y-axis limits
  xlim(0, 20) +  # Restrict x-axis limits
  theme_minimal(base_size = 14) +  # Clean theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = "bottom",  # Move legend to bottom
    legend.key.width = unit(2, "cm")  # Increase legend key size
  )

x11()
p















##### Posterior plot results - PM (n = 1000)



samples_zigzag = sample_copula_IntractLik_gold_rcpp

# True value
#rho = 0.5
#true <- -log(2/(rho+1) - 1)
true = 0.5

# Custom labels for the columns
custom_labels <- c(
  "rho"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(samples_zigzag, "ZZ Gold")
data2 <- convert_to_long(samples_pmcmc_block_m20, "PM 20")
data4 <- convert_to_long(samples_pmcmc_block_m100, "PM 100")
data5 <- convert_to_long(samples_pmcmc_block_m1000, "PM 1000")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)

# Reorder Source factor to control the order in the legend
combined_data$Source <- factor(combined_data$Source, levels = c("PM 20", "PM 100", "PM 1000", "ZZ Gold"))

# Plot densities using line types for each source and vertical dashed lines for true values
plot <- ggplot(combined_data, aes(x = Value)) +
  # Solid black line for Zig Zag (make sure it's in the legend)
  geom_density(data = subset(combined_data, Source == "ZZ Gold"), 
               aes(color = "ZZ Gold", linetype = "ZZ Gold"), 
               size = 1.8, key_glyph = draw_key_path) +  # Increased line width
  # Different line types for other sources
  geom_density(data = subset(combined_data, Source != "ZZ Gold"), 
               aes(color = Source, linetype = Source), size = 1.8, key_glyph = draw_key_path) +  # Increased line width
  # Vertical dashed lines for true values
  geom_vline(xintercept = true, color = "black", linetype = "dashed", size = 1.2) +  # Increased line width for vertical line
  theme_minimal(base_size = 20) +  # Increase base size for better readability
  labs(
    title = "n = 1000",
    x = expression(rho),  # Change x-axis label to "rho"
    y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(values = c("ZZ Gold" = "black", 
                                "PM 20" = "#FF6F61",  # Soft coral red
                                "PM 100" = "#377EB8", # Medium slate blue
                                "PM 1000" = "#32CD32")) +  # Lime green
  scale_linetype_manual(values = c("ZZ Gold" = "solid", 
                                   "PM 20" = "longdash", 
                                   "PM 100" = "twodash",
                                   "PM 1000" = "dashed")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.position = c(0,1),  # Move legend to bottom
    legend.key.width = unit(2, "cm"),
    legend.justification = c(0, 1),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Light gray grid lines
    panel.grid.minor = element_line(color = "gray95", size = 0.25)  # Faint minor grid lines
  ) + 
  xlim(-0.5, 1)+
  ylim(0, 3.5)

# Display the plot
x11()
print(plot)



#### plot results zigzag (n = 1000)


# True value
#rho = 0.5
#true <- -log(2/(rho+1) - 1)
true = 0.5

# Custom labels for the columns
custom_labels <- c(
  "theta"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(sample_copula_IntractLik_m2_rcpp, "ZZ 2")
data2 <- convert_to_long(sample_copula_IntractLik_m10_rcpp, "ZZ 10")
data4 <- convert_to_long(sample_copula_IntractLik_m50_rcpp, "ZZ 50")
data5 <- convert_to_long(sample_copula_IntractLik_gold_rcpp, "ZZ Gold")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)




# Reorder Source factor to control the order in the legend
combined_data$Source <- factor(combined_data$Source, levels = c("ZZ 2", "ZZ 10", "ZZ 50", "ZZ Gold"))

# Plot densities using line types for each source and vertical dashed lines for true values
plot <- ggplot(combined_data, aes(x = Value)) +
  # Solid black line for Zig Zag (make sure it's in the legend)
  geom_density(data = subset(combined_data, Source == "ZZ Gold"), 
               aes(color = "ZZ Gold", linetype = "ZZ Gold"), 
               size = 1.8, key_glyph = draw_key_path) +  # Increased line width
  # Different line types for other sources
  geom_density(data = subset(combined_data, Source != "ZZ Gold"), 
               aes(color = Source, linetype = Source), size = 1.8, key_glyph = draw_key_path) +  # Increased line width
  # Vertical dashed lines for true values
  geom_vline(xintercept = true, color = "black", linetype = "dashed", size = 1.2) +  # Increased line width for vertical line
  theme_minimal(base_size = 20) +  # Increase base size for better readability
  labs(
    title = "n = 1000",
    x = expression(rho),  # Change x-axis label to "rho"
    #y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(values = c("ZZ Gold" = "black", 
                                "ZZ 2" = "#FF6F61",  # Soft coral red
                                "ZZ 10" = "#377EB8", # Medium slate blue
                                "ZZ 50" = "#32CD32")) +  # Lime green
  scale_linetype_manual(values = c("ZZ Gold" = "solid", 
                                   "ZZ 2" = "longdash", 
                                   "ZZ 10" = "twodash",
                                   "ZZ 50" = "dashed")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y=element_blank(),
    legend.position = c(0,1),  # Move legend to bottom
    legend.key.width = unit(2, "cm"),
    legend.justification = c(0, 1),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Light gray grid lines
    panel.grid.minor = element_line(color = "gray95", size = 0.25)  # Faint minor grid lines
  ) + 
  xlim(-0.5, 1) +
  ylim(0, 3.5)

# Display the plot
x11()
print(plot)










##### density plots for n = 100

load(file = "ResultsZZ_copula_n100_Gold.RData")
load(file = "ResultsPM_copula_block_n100.RData")
#load(file = "ResultsPM_copula_n1000.RData")
load(file = "ResultsZZ_copula_n100.RData")



samples_pmcmc_block_m2 = 2/(1 + exp(-samples_pmcmc_block_m2)) - 1
samples_pmcmc_block_m5 = 2/(1 + exp(-samples_pmcmc_block_m5)) - 1
samples_pmcmc_block_m10 = 2/(1 + exp(-samples_pmcmc_block_m10)) - 1
samples_pmcmc_block_m50 = 2/(1 + exp(-samples_pmcmc_block_m50)) - 1
samples_pmcmc_block_m20 = 2/(1 + exp(-samples_pmcmc_block_m20)) - 1

sample_copula_IntractLik_m2_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m2_rcpp)) - 1
sample_copula_IntractLik_m5_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m5_rcpp)) - 1
sample_copula_IntractLik_m10_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m10_rcpp)) - 1
sample_copula_IntractLik_m20_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m20_rcpp)) - 1
sample_copula_IntractLik_m50_rcpp = 2/(1 + exp(-sample_copula_IntractLik_m50_rcpp)) - 1

sample_copula_IntractLik_gold_rcpp = 2/(1 + exp(-sample_copula_IntractLik_gold_rcpp)) - 1




##### Posterior plot results - PM (n = 100)



samples_zigzag = sample_copula_IntractLik_gold_rcpp

# True value
#rho = 0.5
#true <- -log(2/(rho+1) - 1)
true = 0.5

# Custom labels for the columns
custom_labels <- c(
  "rho"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(samples_zigzag, "ZZ Gold")
data2 <- convert_to_long(samples_pmcmc_block_m2, "PM 2")
data4 <- convert_to_long(samples_pmcmc_block_m10, "PM 10")
data5 <- convert_to_long(samples_pmcmc_block_m50, "PM 50")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)

# Reorder Source factor to control the order in the legend
combined_data$Source <- factor(combined_data$Source, levels = c("PM 2", "PM 10", "PM 50", "ZZ Gold"))

# Plot densities using line types for each source and vertical dashed lines for true values
plot <- ggplot(combined_data, aes(x = Value)) +
  # Solid black line for Zig Zag (make sure it's in the legend)
  geom_density(data = subset(combined_data, Source == "ZZ Gold"), 
               aes(color = "ZZ Gold", linetype = "ZZ Gold"), 
               size = 1.8, key_glyph = draw_key_path) +  # Increased line width
  # Different line types for other sources
  geom_density(data = subset(combined_data, Source != "ZZ Gold"), 
               aes(color = Source, linetype = Source), size = 1.8, key_glyph = draw_key_path) +  # Increased line width
  # Vertical dashed lines for true values
  geom_vline(xintercept = true, color = "black", linetype = "dashed", size = 1.2) +  # Increased line width for vertical line
  theme_minimal(base_size = 20) +  # Increase base size for better readability
  labs(
    title = "n = 100",
    x = expression(rho),  # Change x-axis label to "rho"
    y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(values = c("ZZ Gold" = "black", 
                                "PM 2" = "#FF6F61",  # Soft coral red
                                "PM 10" = "#377EB8", # Medium slate blue
                                "PM 50" = "#32CD32")) +  # Lime green
  scale_linetype_manual(values = c("ZZ Gold" = "solid", 
                                   "PM 2" = "longdash", 
                                   "PM 10" = "twodash",
                                   "PM 50" = "dashed")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.position = c(0,1),  # Move legend to bottom
    legend.key.width = unit(2, "cm"),
    legend.justification = c(0, 1),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Light gray grid lines
    panel.grid.minor = element_line(color = "gray95", size = 0.25)  # Faint minor grid lines
  ) + 
  xlim(-1.2, 1.2)+
  ylim(0, 1.75)

# Display the plot
x11()
print(plot)




#### plot results zigzag (n = 100)


# True value
#rho = 0.5
#true <- -log(2/(rho+1) - 1)
true = 0.5

# Custom labels for the columns
custom_labels <- c(
  "theta"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(sample_copula_IntractLik_m2_rcpp, "ZZ 2")
data2 <- convert_to_long(sample_copula_IntractLik_m10_rcpp, "ZZ 10")
data4 <- convert_to_long(sample_copula_IntractLik_m50_rcpp, "ZZ 50")
data5 <- convert_to_long(sample_copula_IntractLik_gold_rcpp, "ZZ Gold")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)




# Reorder Source factor to control the order in the legend
combined_data$Source <- factor(combined_data$Source, levels = c("ZZ 2", "ZZ 10", "ZZ 50", "ZZ Gold"))

# Plot densities using line types for each source and vertical dashed lines for true values
plot <- ggplot(combined_data, aes(x = Value)) +
  # Solid black line for Zig Zag (make sure it's in the legend)
  geom_density(data = subset(combined_data, Source == "ZZ Gold"), 
               aes(color = "ZZ Gold", linetype = "ZZ Gold"), 
               size = 1.8, key_glyph = draw_key_path) +  # Increased line width
  # Different line types for other sources
  geom_density(data = subset(combined_data, Source != "ZZ Gold"), 
               aes(color = Source, linetype = Source), size = 1.8, key_glyph = draw_key_path) +  # Increased line width
  # Vertical dashed lines for true values
  geom_vline(xintercept = true, color = "black", linetype = "dashed", size = 1.2) +  # Increased line width for vertical line
  theme_minimal(base_size = 20) +  # Increase base size for better readability
  labs(
    title = "n = 100",
    x = expression(rho),  # Change x-axis label to "rho"
    #y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(values = c("ZZ Gold" = "black", 
                                "ZZ 2" = "#FF6F61",  # Soft coral red
                                "ZZ 10" = "#377EB8", # Medium slate blue
                                "ZZ 50" = "#32CD32")) +  # Lime green
  scale_linetype_manual(values = c("ZZ Gold" = "solid", 
                                   "ZZ 2" = "longdash", 
                                   "ZZ 10" = "twodash",
                                   "ZZ 50" = "dashed")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title.y=element_blank(),
    legend.position = c(0,1),  # Move legend to bottom
    legend.key.width = unit(2, "cm"),
    legend.justification = c(0, 1),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Light gray grid lines
    panel.grid.minor = element_line(color = "gray95", size = 0.25)  # Faint minor grid lines
  ) + 
  xlim(-1.2, 1.2)+
  ylim(0, 1.75)

# Display the plot
x11()
print(plot)











