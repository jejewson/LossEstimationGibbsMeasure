

library(ggplot2)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(GGally)
library(tidyr)

load(file = "ResultsZZ_n100_Gold.RData")
load(file = "ResultsPMBlock.RData")
#load(file = "ResultsPM.RData")
load(file = "ResultsZZ.RData")


gold_mean = colMeans(sample_Gaussian_regression_IntractLik_gold)
gold_sd = colSds(sample_Gaussian_regression_IntractLik_gold)


thin = 10


sims_pm_block_m2 = (1:length(samples_pmcmc_block_m2[,1]))*2
length_block_m2 = length(sims_pm_block_m2)

sims_pm_block_m2 = sims_pm_block_m2[seq(thin, length_block_m2, by = thin)]
  
time_pm_block_m2 = c(0,seq(as.numeric(time.taken.m2)/length(samples_pmcmc_block_m2[,1]),as.numeric(time.taken.m2),by=as.numeric(time.taken.m2)/length(samples_pmcmc_block_m2[,1])))
time_pm_block_m2 = time_pm_block_m2[seq(thin, length_block_m2, by = thin)]

running_mean_pm_block_m2 <- apply(samples_pmcmc_block_m2, 2, function(x) cumsum(x) / seq_along(x))
running_mean_pm_block_m2 = running_mean_pm_block_m2[seq(thin, length_block_m2, by = thin),]

running_var_pm_block_m2 <- apply(samples_pmcmc_block_m2, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m2 <- sqrt(running_var_pm_block_m2)
running_sd_pm_block_m2 = running_sd_pm_block_m2[seq(thin, length_block_m2, by = thin),]


sims_pm_block_m5 = (1:length(samples_pmcmc_block_m5[,1]))*5
length_block_m5 = length(sims_pm_block_m5)

sims_pm_block_m5 = sims_pm_block_m2[seq(thin, length_block_m5, by = thin)]

time_pm_block_m5 = seq(as.numeric(time.taken.m5)/length(samples_pmcmc_block_m5[,1]),as.numeric(time.taken.m5),by=as.numeric(time.taken.m5)/length(samples_pmcmc_block_m5[,1]))
time_pm_block_m5 = time_pm_block_m5[seq(thin, length_block_m5, by = thin)]

running_mean_pm_block_m5 <- apply(samples_pmcmc_block_m5, 2, function(x) cumsum(x) / seq_along(x))
running_mean_pm_block_m5 = running_mean_pm_block_m5[seq(thin, length_block_m5, by = thin),]

running_var_pm_block_m5 <- apply(samples_pmcmc_block_m5, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m5 <- sqrt(running_var_pm_block_m5)
running_sd_pm_block_m5 = running_sd_pm_block_m5[seq(thin, length_block_m5, by = thin),]


sims_pm_block_m20 = (1:length(samples_pmcmc_block_m20[,1]))*20
length_block_m20 = length(sims_pm_block_m20)

sims_pm_block_m20 = sims_pm_block_m20[seq(thin, length_block_m20, by = thin)]

time_pm_block_m20 = seq(as.numeric(time.taken.m20)/length(samples_pmcmc_block_m20[,1]),as.numeric(time.taken.m20),by=as.numeric(time.taken.m20)/length(samples_pmcmc_block_m20[,1]))
time_pm_block_m20 = time_pm_block_m20[seq(thin, length_block_m20, by = thin)]

running_mean_pm_block_m20 <- apply(samples_pmcmc_block_m20, 2, function(x) cumsum(x) / seq_along(x))
running_mean_pm_block_m20 = running_mean_pm_block_m20[seq(thin, length_block_m20, by = thin),]

running_var_pm_block_m20 <- apply(samples_pmcmc_block_m20, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m20 <- sqrt(running_var_pm_block_m20)
running_sd_pm_block_m20 = running_sd_pm_block_m20[seq(thin, length_block_m20, by = thin),]


sims_pm_block_m100 = (1:length(samples_pmcmc_block_m100[,1]))*100
length_block_m100 = length(sims_pm_block_m100)

sims_pm_block_m100 = sims_pm_block_m100[seq(thin, length_block_m100, by = thin)]

time_pm_block_m100 = seq(as.numeric(time.taken.m100)*60/length(samples_pmcmc_block_m100[,1]),as.numeric(time.taken.m100)*60,by=as.numeric(time.taken.m100)*60/length(samples_pmcmc_block_m100[,1]))
time_pm_block_m100 = time_pm_block_m100[seq(thin, length_block_m100, by = thin)]

running_mean_pm_block_m100 <- apply(samples_pmcmc_block_m100, 2, function(x) cumsum(x) / seq_along(x))
running_mean_pm_block_m100 = running_mean_pm_block_m100[seq(thin, length_block_m100, by = thin),]

running_var_pm_block_m100 <- apply(samples_pmcmc_block_m100, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_pm_block_m100 <- sqrt(running_var_pm_block_m100)
running_sd_pm_block_m100 = running_sd_pm_block_m100[seq(thin, length_block_m100, by = thin),]



num_sims_zz_m2 = length(skeleton_Gaussian_regression_IntractLik_m2$skeleton_T)*1/skeleton_Gaussian_regression_IntractLik_m2$mean_alpha*2

sims_zz_m2 = seq(1, num_sims_zz_m2, by = num_sims_zz_m2/length(sample_Gaussian_regression_IntractLik_m2[,1]))
length_zz_m2 = length(sims_zz_m2)

sims_zz_m2 = sims_zz_m2[seq(thin, length_zz_m2, by = thin)]

time_zz_m2 = seq(as.numeric(timer_MMD_Gaussian_regression_m2)/length(sample_Gaussian_regression_IntractLik_m2[,1]),as.numeric(timer_MMD_Gaussian_regression_m2),by=as.numeric(timer_MMD_Gaussian_regression_m2)/length(sample_Gaussian_regression_IntractLik_m2[,1]))
time_zz_m2 = time_zz_m2[seq(thin, length_zz_m2, by = thin)]


running_mean_zz_m2 <- apply(sample_Gaussian_regression_IntractLik_m2, 2, function(x) cumsum(x) / seq_along(x))
running_mean_zz_m2 = running_mean_zz_m2[seq(thin, length_zz_m2, by = thin),]

running_var_zz_m2 <- apply(sample_Gaussian_regression_IntractLik_m2, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m2 <- sqrt(running_var_zz_m2)
running_sd_zz_m2 = running_sd_zz_m2[seq(thin, length_zz_m2, by = thin),]



num_sims_zz_m5 = length(skeleton_Gaussian_regression_IntractLik_m5$skeleton_T)*1/skeleton_Gaussian_regression_IntractLik_m5$mean_alpha*2
sims_zz_m5 = seq(1, num_sims_zz_m5, by = num_sims_zz_m5/length(sample_Gaussian_regression_IntractLik_m5[,1]))
length_zz_m5 = length(sims_zz_m5)

sims_zz_m5 = sims_zz_m5[seq(thin, length_zz_m5, by = thin)]

time_zz_m5 = seq(as.numeric(timer_MMD_Gaussian_regression_m5)/length(sample_Gaussian_regression_IntractLik_m5[,1]),as.numeric(timer_MMD_Gaussian_regression_m5),by=as.numeric(timer_MMD_Gaussian_regression_m5)/length(sample_Gaussian_regression_IntractLik_m5[,1]))
time_zz_m5 = time_zz_m5[seq(thin, length_zz_m5, by = thin)]


running_mean_zz_m5 <- apply(sample_Gaussian_regression_IntractLik_m5, 2, function(x) cumsum(x) / seq_along(x))
running_mean_zz_m5 = running_mean_zz_m5[seq(thin, length_zz_m5, by = thin),]

running_var_zz_m5 <- apply(sample_Gaussian_regression_IntractLik_m5, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m5 <- sqrt(running_var_zz_m5)
running_sd_zz_m5 = running_sd_zz_m5[seq(thin, length_zz_m5, by = thin),]


num_sims_zz_m10 = length(skeleton_Gaussian_regression_IntractLik_m10$skeleton_T)*1/skeleton_Gaussian_regression_IntractLik_m10$mean_alpha*2
sims_zz_m10 = seq(1, num_sims_zz_m10, by = num_sims_zz_m10/length(sample_Gaussian_regression_IntractLik_m10[,1]))
length_zz_m10 = length(sims_zz_m10)

sims_zz_m10 = sims_zz_m2[seq(thin, length_zz_m10, by = thin)]

time_zz_m10 = seq(as.numeric(timer_MMD_Gaussian_regression_m10)/length(sample_Gaussian_regression_IntractLik_m10[,1]),as.numeric(timer_MMD_Gaussian_regression_m10),by=as.numeric(timer_MMD_Gaussian_regression_m10)/length(sample_Gaussian_regression_IntractLik_m10[,1]))
time_zz_m10 = time_zz_m10[seq(thin, length_zz_m10, by = thin)]


running_mean_zz_m10 <- apply(sample_Gaussian_regression_IntractLik_m10, 2, function(x) cumsum(x) / seq_along(x))
running_mean_zz_m10 = running_mean_zz_m10[seq(thin, length_zz_m10, by = thin),]

running_var_zz_m10 <- apply(sample_Gaussian_regression_IntractLik_m10, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m10 <- sqrt(running_var_zz_m10)
running_sd_zz_m10 = running_sd_zz_m10[seq(thin, length_zz_m10, by = thin),]


num_sims_zz_m20 = length(skeleton_Gaussian_regression_IntractLik_m20$skeleton_T)*1/skeleton_Gaussian_regression_IntractLik_m20$mean_alpha*2
sims_zz_m20 = seq(1, num_sims_zz_m20, by = num_sims_zz_m20/length(sample_Gaussian_regression_IntractLik_m20[,1]))
length_zz_m20 = length(sims_zz_m20)

sims_zz_m20 = sims_zz_m20[seq(thin, length_zz_m20, by = thin)]

time_zz_m20 = seq(as.numeric(timer_MMD_Gaussian_regression_m20)/length(sample_Gaussian_regression_IntractLik_m20[,1]),as.numeric(timer_MMD_Gaussian_regression_m20),by=as.numeric(timer_MMD_Gaussian_regression_m20)/length(sample_Gaussian_regression_IntractLik_m20[,1]))
time_zz_m20 = time_zz_m20[seq(thin, length_zz_m20, by = thin)]

running_mean_zz_m20 <- apply(sample_Gaussian_regression_IntractLik_m20, 2, function(x) cumsum(x) / seq_along(x))
running_mean_zz_m20 = running_mean_zz_m20[seq(thin, length_zz_m20, by = thin),]

running_var_zz_m20 <- apply(sample_Gaussian_regression_IntractLik_m20, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m20 <- sqrt(running_var_zz_m20)
running_sd_zz_m20 = running_sd_zz_m20[seq(thin, length_zz_m20, by = thin),]



num_sims_zz_m50 = length(skeleton_Gaussian_regression_IntractLik_m50$skeleton_T)*1/skeleton_Gaussian_regression_IntractLik_m50$mean_alpha*2
sims_zz_m50 = seq(1, num_sims_zz_m50, by = num_sims_zz_m50/length(sample_Gaussian_regression_IntractLik_m50[,1]))
length_zz_m50 = length(sims_zz_m50)

sims_zz_m50 = sims_zz_m50[seq(thin, length_zz_m50, by = thin)]

time_zz_m50 = seq(as.numeric(timer_MMD_Gaussian_regression_m50)/length(sample_Gaussian_regression_IntractLik_m50[,1]),as.numeric(timer_MMD_Gaussian_regression_m50),by=as.numeric(timer_MMD_Gaussian_regression_m50)/length(sample_Gaussian_regression_IntractLik_m50[,1]))
time_zz_m50 = time_zz_m50[seq(thin, length_zz_m50, by = thin)]

running_mean_zz_m50 <- apply(sample_Gaussian_regression_IntractLik_m50, 2, function(x) cumsum(x) / seq_along(x))
running_mean_zz_m50 = running_mean_zz_m50[seq(thin, length_zz_m50, by = thin),]

running_var_zz_m50 <- apply(sample_Gaussian_regression_IntractLik_m50, 2, function(x) (cumsum(x^2) - (cumsum(x)^2) / seq_along(x)) / (seq_along(x) - 1))
running_sd_zz_m50 <- sqrt(running_var_zz_m50)
running_sd_zz_m50 = running_sd_zz_m50[seq(thin, length_zz_m50, by = thin),]








############### PM plots for the first moment vs number of model simulations

param = 1 # change parameter number as desired

# Create a dataframe for ggplot
data <- data.frame(
  x = c(sims_pm_block_m2, sims_pm_block_m5, sims_pm_block_m20, sims_pm_block_m100),
  y = c(running_mean_pm_block_m2[,param] - gold_mean[param], 
        running_mean_pm_block_m5[,param] - gold_mean[param], 
        running_mean_pm_block_m20[,param] - gold_mean[param], 
        running_mean_pm_block_m100[,param] - gold_mean[param]),
  group = rep(c("PM 2", "PM 5", "PM 20", "PM 100"), 
              times = c(length(sims_pm_block_m2), length(sims_pm_block_m5), 
                        length(sims_pm_block_m20), length(sims_pm_block_m100)))
)

# Ensure the factor levels are in the correct order for legend appearance
data$group <- factor(data$group, levels = c("PM 2", "PM 5", "PM 20", "PM 100"))

# Define colors and line types
color_map <- c("PM 2" = "#FF6F61", "PM 5" = "#377EB8", "PM 20" = "#32CD32", "PM 100" = "#6A3D9A")
line_map <- c("PM 2" = "solid", "PM 5" = "dashed", "PM 20" = "dotted", "PM 100" = "dotdash")

# Create the plot
p = ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(y = "Estimate - Gold") +
  ylim(-0.1, 0.1) +  # Restrict y-axis limits
  xlim(0, 1e7) +  # Restrict x-axis limits
  theme_minimal(base_size = 20) +  # Clean theme
  theme(
    axis.title.x=element_blank(),
    legend.title = element_blank(),  # Remove legend title
    legend.position = c(1,1), 
    legend.key.width = unit(2, "cm"),
    legend.justification = c(1, 1)
  )


x11()
p




############### PM plots for the second moment vs number of model simulations

param = 1 # change parameter number as desired



# Create a dataframe for ggplot
data <- data.frame(
  x = c(sims_pm_block_m2, sims_pm_block_m5, sims_pm_block_m20, sims_pm_block_m100),
  y = c(running_sd_pm_block_m2[,param] - gold_sd[param], 
        running_sd_pm_block_m5[,param] - gold_sd[param], 
        running_sd_pm_block_m20[,param] - gold_sd[param], 
        running_sd_pm_block_m100[,param] - gold_sd[param]),
  group = rep(c("PM 2", "PM 5", "PM 20", "PM 100"), 
              times = c(length(sims_pm_block_m2), length(sims_pm_block_m5), 
                        length(sims_pm_block_m20), length(sims_pm_block_m100)))
)

# Ensure the factor levels are in the correct order for legend appearance
data$group <- factor(data$group, levels = c("PM 2", "PM 5", "PM 20", "PM 100"))

# Define colors and line types
color_map <- c("PM 2" = "#FF6F61", "PM 5" = "#377EB8", "PM 20" = "#32CD32", "PM 100" = "#6A3D9A")
line_map <- c("PM 2" = "solid", "PM 5" = "dashed", "PM 20" = "dotted", "PM 100" = "dotdash")

# Create the plot
p = ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Number of model simulations", y = "Estimate - Gold") +
  ylim(-0.1, 0.1) +  # Restrict y-axis limits
  xlim(0, 1e7) +  # Restrict x-axis limits
  theme_minimal(base_size = 20) +  # Clean theme
  theme(
    legend.title = element_blank(),  # Remove legend title
    legend.position = c(1,1), 
    legend.key.width = unit(2, "cm"),
    legend.justification = c(1, 1)
  )


x11()
p




############### PM plots for the first moment vs TIME

param = 1 # change parameter number as desired



# Create a dataframe for ggplot
data <- data.frame(
  x = c(time_pm_block_m2, time_pm_block_m5, time_pm_block_m20, time_pm_block_m100),
  y = c(running_mean_pm_block_m2[,param] - gold_mean[param], 
        running_mean_pm_block_m5[,param] - gold_mean[param], 
        running_mean_pm_block_m20[,param] - gold_mean[param], 
        running_mean_pm_block_m100[,param] - gold_mean[param]),
  group = rep(c("PM 2", "PM 5", "PM 20", "PM 100"), 
              times = c(length(time_pm_block_m2), length(time_pm_block_m5), 
                        length(time_pm_block_m20), length(time_pm_block_m100)))
)

# Ensure the factor levels are in the correct order for legend appearance
data$group <- factor(data$group, levels = c("PM 2", "PM 5", "PM 20", "PM 100"))

# Define colors and line types
color_map <- c("PM 2" = "#FF6F61", "PM 5" = "#377EB8", "PM 20" = "#32CD32", "PM 100" = "#6A3D9A")
line_map <- c("PM 2" = "solid", "PM 5" = "dashed", "PM 20" = "dotted", "PM 100" = "dotdash")

# Create the plot
p = ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Time (mins)", y = "Estimate - Gold") +
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




############### PM plots for the second moment vs TIME

param = 1 # change parameter number as desired



# Create a dataframe for ggplot
data <- data.frame(
  x = c(time_pm_block_m2, time_pm_block_m5, time_pm_block_m20, time_pm_block_m100),
  y = c(running_sd_pm_block_m2[,param] - gold_sd[param], 
        running_sd_pm_block_m5[,param] - gold_sd[param], 
        running_sd_pm_block_m20[,param] - gold_sd[param], 
        running_sd_pm_block_m100[,param] - gold_sd[param]),
  group = rep(c("PM 2", "PM 5", "PM 20", "PM 100"), 
              times = c(length(time_pm_block_m2), length(time_pm_block_m5), 
                        length(time_pm_block_m20), length(time_pm_block_m100)))
)

# Ensure the factor levels are in the correct order for legend appearance
data$group <- factor(data$group, levels = c("PM 2", "PM 5", "PM 20", "PM 100"))

# Define colors and line types
color_map <- c("PM 2" = "#FF6F61", "PM 5" = "#377EB8", "PM 20" = "#32CD32", "PM 100" = "#6A3D9A")
line_map <- c("PM 2" = "solid", "PM 5" = "dashed", "PM 20" = "dotted", "PM 100" = "dotdash")

# Create the plot
p = ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Time (mins)", y = "Estimate - Gold") +
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
    legend.position = c(1,1), 
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
  ylim(-0.1, 0.1) +  # Restrict y-axis limits
  xlim(0, 1e7) +  # Restrict x-axis limits
  theme_minimal(base_size = 20) +  # Clean theme
  theme(
    axis.title.y=element_blank(),
    legend.title = element_blank(),  # Remove legend title
    legend.position = c(1,1), 
    legend.key.width = unit(2, "cm"),
    legend.justification = c(1, 1)
  )

x11()
p






############### ZZ plots for the first moment vs TIME


# Create a dataframe for ggplot
data <- data.frame(
  x = c(time_zz_m2, time_zz_m5, time_zz_m10, time_zz_m20),
  y = c(running_mean_zz_m2[,param] - gold_mean[param], 
        running_mean_zz_m5[,param] - gold_mean[param], 
        running_mean_zz_m10[,param] - gold_mean[param], 
        running_mean_zz_m20[,param] - gold_mean[param]),
  group = rep(c("ZZ 2", "ZZ 5", "ZZ 10", "ZZ 20"), 
              times = c(length(time_zz_m2), length(time_zz_m5), 
                        length(time_zz_m10), length(time_zz_m20)))
)

# Explicitly order the factor levels
data$group <- factor(data$group, levels = c("ZZ 2", "ZZ 5", "ZZ 10", "ZZ 20"))

# Define colors and line types
color_map <- c("ZZ 2" = "#FF6F61", "ZZ 5" = "#377EB8", "ZZ 10" = "#32CD32", "ZZ 20" = "#6A3D9A")
line_map <- c("ZZ 2" = "solid", "ZZ 5" = "dashed", "ZZ 10" = "dotted", "ZZ 20" = "dotdash")

# Create the plot
p <- ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Time (mins)", y = "Estimate - Gold") +
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




############### ZZ plots for the second moment vs TIME


# Create a dataframe for ggplot
data <- data.frame(
  x = c(time_zz_m2, time_zz_m5, time_zz_m10, time_zz_m20),
  y = c(running_sd_zz_m2[,param] - gold_sd[param], 
        running_sd_zz_m5[,param] - gold_sd[param], 
        running_sd_zz_m10[,param] - gold_sd[param], 
        running_sd_zz_m20[,param] - gold_sd[param]),
  group = rep(c("ZZ 2", "ZZ 5", "ZZ 10", "ZZ 20"), 
              times = c(length(time_zz_m2), length(time_zz_m5), 
                        length(time_zz_m10), length(time_zz_m20)))
)

# Explicitly order the factor levels
data$group <- factor(data$group, levels = c("ZZ 2", "ZZ 5", "ZZ 10", "ZZ 20"))

# Define colors and line types
color_map <- c("ZZ 2" = "#FF6F61", "ZZ 5" = "#377EB8", "ZZ 10" = "#32CD32", "ZZ 20" = "#6A3D9A")
line_map <- c("ZZ 2" = "solid", "ZZ 5" = "dashed", "ZZ 10" = "dotted", "ZZ 20" = "dotdash")

# Create the plot
p <- ggplot(data, aes(x = x, y = y, color = group, linetype = group)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = color_map) + 
  scale_linetype_manual(values = line_map) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) + # Gold reference line
  labs(x = "Time (mins)", y = "Estimate - Gold") +
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






######## density plots PM


# True values for each column
true <- c(4, 4, 3, 3, 2, 2, 1, 1, 0)

# Custom labels for the columns
custom_labels <- c(
  "beta_1", "beta_2", "beta_3", "beta_4", 
  "beta_5", "beta_6", "beta_7", "beta_8", "log sigma"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(sample_Gaussian_regression_IntractLik_gold, "ZZ Gold")
data2 <- convert_to_long(samples_pmcmc_block_m2, "PM 2")
data4 <- convert_to_long(samples_pmcmc_block_m10, "PM 10")
data5 <- convert_to_long(samples_pmcmc_block_m20, "PM 20")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)

# Create a data frame for true values
true_values <- data.frame(
  Column = paste0("V", 1:9),  # Match the column names after pivoting
  TrueValue = true
)

# Replace Column with custom labels
combined_data$CustomColumn <- factor(
  combined_data$Column,
  levels = paste0("V", 1:9),
  labels = custom_labels
)

true_values$CustomColumn <- factor(
  true_values$Column,
  levels = paste0("V", 1:9),
  labels = custom_labels
)

# Define the correct order for the groups (PM 2, PM 10, PM 20, ZZ Gold)
combined_data$Source <- factor(combined_data$Source, 
                               levels = c("PM 2", "PM 10", "PM 20", "ZZ Gold"))

# Plot densities using line types for each source and vertical dashed lines for true values
plot <- ggplot(combined_data, aes(x = Value)) +
  # Solid black line for Zig Zag (make sure it's in the legend)
  geom_density(data = subset(combined_data, Source == "ZZ Gold"), 
               aes(color = Source, linetype = Source), 
               size = 1.2) +  
  # Different line types for other sources
  geom_density(data = subset(combined_data, Source != "ZZ Gold"), 
               aes(color = Source, linetype = Source), size = 1.2) +
  # Vertical dashed lines for true values
  geom_vline(data = true_values, aes(xintercept = TrueValue), color = "black", linetype = "dashed", size = 0.8) +
  facet_wrap(~ CustomColumn, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(
    x = "Value",
    y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(
    values = c("PM 2" = "#FF6F61", 
               "PM 10" = "#377EB8", 
               "PM 20" = "#32CD32", 
               "ZZ Gold" = "black")
  ) +
  scale_linetype_manual(
    values = c("PM 2" = "longdash", 
               "PM 10" = "twodash", 
               "PM 20" = "dashed", 
               "ZZ Gold" = "solid")
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "top"
  ) 

# Display the plot
x11()
print(plot)






######## density plots ZZ


# True values for each column
true <- c(4, 4, 3, 3, 2, 2, 1, 1, 0)

# Custom labels for the columns
custom_labels <- c(
  "beta_1", "beta_2", "beta_3", "beta_4", 
  "beta_5", "beta_6", "beta_7", "beta_8", "log sigma"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(sample_Gaussian_regression_IntractLik_gold, "ZZ Gold")
data2 <- convert_to_long(sample_Gaussian_regression_IntractLik_m2, "ZZ 2")
data4 <- convert_to_long(sample_Gaussian_regression_IntractLik_m5, "ZZ 10")
data5 <- convert_to_long(sample_Gaussian_regression_IntractLik_m20, "ZZ 20")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)

# Create a data frame for true values
true_values <- data.frame(
  Column = paste0("V", 1:9),  # Match the column names after pivoting
  TrueValue = true
)

# Replace Column with custom labels
combined_data$CustomColumn <- factor(
  combined_data$Column,
  levels = paste0("V", 1:9),
  labels = custom_labels
)

true_values$CustomColumn <- factor(
  true_values$Column,
  levels = paste0("V", 1:9),
  labels = custom_labels
)

# Define the correct order for the groups (ZZ 2, ZZ 10, ZZ 20, ZZ Gold)
combined_data$Source <- factor(combined_data$Source, 
                               levels = c("ZZ 2", "ZZ 10", "ZZ 20", "ZZ Gold"))

# Plot densities using line types for each source and vertical dashed lines for true values
plot <- ggplot(combined_data, aes(x = Value)) +
  # Solid black line for Zig Zag (make sure it's in the legend)
  geom_density(data = subset(combined_data, Source == "ZZ Gold"), 
               aes(color = Source, linetype = Source), 
               size = 1.2) +  
  # Different line types for other sources
  geom_density(data = subset(combined_data, Source != "ZZ Gold"), 
               aes(color = Source, linetype = Source), size = 1.2) +
  # Vertical dashed lines for true values
  geom_vline(data = true_values, aes(xintercept = TrueValue), color = "black", linetype = "dashed", size = 0.8) +
  facet_wrap(~ CustomColumn, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(
    x = "Value",
    y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(
    values = c("ZZ 2" = "#FF6F61", 
               "ZZ 10" = "#377EB8", 
               "ZZ 20" = "#32CD32", 
               "ZZ Gold" = "black")
  ) +
  scale_linetype_manual(
    values = c("ZZ 2" = "longdash", 
               "ZZ 10" = "twodash", 
               "ZZ 20" = "dashed", 
               "ZZ Gold" = "solid")
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "top"
  ) 

# Display the plot
x11()
print(plot)








##### Posterior plot results - PM beta_1



samples_zigzag = sample_Gaussian_regression_IntractLik_gold[,1]

# True value
true = 4

# Custom labels for the columns
custom_labels <- c(
  "beta_1"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(samples_zigzag, "ZZ Gold")
data2 <- convert_to_long(samples_pmcmc_block_m2[,1], "PM 2")
data4 <- convert_to_long(samples_pmcmc_block_m5[,1], "PM 5")
data5 <- convert_to_long(samples_pmcmc_block_m20[,1], "PM 20")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)

# Reorder Source factor to control the order in the legend
combined_data$Source <- factor(combined_data$Source, levels = c("PM 2", "PM 5", "PM 20", "ZZ Gold"))

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
    x = expression(beta[1]),  # Change x-axis label to "rho"
    y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(values = c("ZZ Gold" = "black", 
                                "PM 2" = "#FF6F61",  # Soft coral red
                                "PM 5" = "#377EB8", # Medium slate blue
                                "PM 20" = "#32CD32")) +  # Lime green
  scale_linetype_manual(values = c("ZZ Gold" = "solid", 
                                   "PM 2" = "longdash", 
                                   "PM 5" = "twodash",
                                   "PM 20" = "dashed")) +
  theme(
    legend.position = c(1,1),  
    legend.key.width = unit(2, "cm"),
    legend.justification = c(1, 1),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Light gray grid lines
    panel.grid.minor = element_line(color = "gray95", size = 0.25)  # Faint minor grid lines
  ) + 
  xlim(3, 6)+
  ylim(0, 1.5)


# Display the plot
x11()
print(plot)



##### Posterior plot results - PM log sigma



samples_zigzag = sample_Gaussian_regression_IntractLik_gold[,9]

# True value
true = 0

# Custom labels for the columns
custom_labels <- c(
  "log sigma"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(samples_zigzag, "ZZ Gold")
data2 <- convert_to_long(samples_pmcmc_block_m2[,9], "PM 2")
data4 <- convert_to_long(samples_pmcmc_block_m5[,9], "PM 5")
data5 <- convert_to_long(samples_pmcmc_block_m20[,9], "PM 20")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)

# Reorder Source factor to control the order in the legend
combined_data$Source <- factor(combined_data$Source, levels = c("PM 2", "PM 5", "PM 20", "ZZ Gold"))

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
    x = expression(log(sigma)),
    y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(values = c("ZZ Gold" = "black", 
                                "PM 2" = "#FF6F61",  # Soft coral red
                                "PM 5" = "#377EB8", # Medium slate blue
                                "PM 20" = "#32CD32")) +  # Lime green
  scale_linetype_manual(values = c("ZZ Gold" = "solid", 
                                   "PM 2" = "longdash", 
                                   "PM 5" = "twodash",
                                   "PM 20" = "dashed")) +
  theme(
    legend.position = c(1,1),  
    legend.key.width = unit(2, "cm"),
    legend.justification = c(1, 1),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Light gray grid lines
    panel.grid.minor = element_line(color = "gray95", size = 0.25)  # Faint minor grid lines
  ) + 
xlim(-0.5, 1.5)+
ylim(0, 2.5)


# Display the plot
x11()
print(plot)









######## density plots ZZ


# True values for each column
true <- c(4, 4, 3, 3, 2, 2, 1, 1, 0)

# Custom labels for the columns
custom_labels <- c(
  "beta_1", "beta_2", "beta_3", "beta_4", 
  "beta_5", "beta_6", "beta_7", "beta_8", "log sigma"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(sample_Gaussian_regression_IntractLik_gold, "ZZ Gold")
data2 <- convert_to_long(sample_Gaussian_regression_IntractLik_m2, "ZZ 2")
data4 <- convert_to_long(sample_Gaussian_regression_IntractLik_m5, "ZZ 5")
data5 <- convert_to_long(sample_Gaussian_regression_IntractLik_m20, "ZZ 20")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)

# Create a data frame for true values
true_values <- data.frame(
  Column = paste0("V", 1:9),  # Match the column names after pivoting
  TrueValue = true
)

# Replace Column with custom labels
combined_data$CustomColumn <- factor(
  combined_data$Column,
  levels = paste0("V", 1:9),
  labels = custom_labels
)

true_values$CustomColumn <- factor(
  true_values$Column,
  levels = paste0("V", 1:9),
  labels = custom_labels
)

# Define the correct order for the groups (ZZ 2, ZZ 5, ZZ 20, ZZ Gold)
combined_data$Source <- factor(combined_data$Source, 
                               levels = c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ Gold"))

# Plot densities using line types for each source and vertical dashed lines for true values
plot <- ggplot(combined_data, aes(x = Value)) +
  # Solid black line for Zig Zag (make sure it's in the legend)
  geom_density(data = subset(combined_data, Source == "ZZ Gold"), 
               aes(color = Source, linetype = Source), 
               size = 1.2) +  
  # Different line types for other sources
  geom_density(data = subset(combined_data, Source != "ZZ Gold"), 
               aes(color = Source, linetype = Source), size = 1.2) +
  # Vertical dashed lines for true values
  geom_vline(data = true_values, aes(xintercept = TrueValue), color = "black", linetype = "dashed", size = 0.8) +
  facet_wrap(~ CustomColumn, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(
    x = "Value",
    y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(
    values = c("ZZ 2" = "#FF6F61", 
               "ZZ 5" = "#377EB8", 
               "ZZ 20" = "#32CD32", 
               "ZZ Gold" = "black")
  ) +
  scale_linetype_manual(
    values = c("ZZ 2" = "longdash", 
               "ZZ 5" = "twodash", 
               "ZZ 20" = "dashed", 
               "ZZ Gold" = "solid")
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "top"
  ) 

# Display the plot
x11()
print(plot)










##### Posterior plot results - ZZ beta_1



samples_zigzag = sample_Gaussian_regression_IntractLik_gold[,1]

# True value
true = 4

# Custom labels for the columns
custom_labels <- c(
  "beta_1"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(samples_zigzag, "ZZ Gold")
data2 <- convert_to_long(sample_Gaussian_regression_IntractLik_m2[,1], "ZZ 2")
data4 <- convert_to_long(sample_Gaussian_regression_IntractLik_m5[,1], "ZZ 5")
data5 <- convert_to_long(sample_Gaussian_regression_IntractLik_m20[,1], "ZZ 20")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)

# Reorder Source factor to control the order in the legend
combined_data$Source <- factor(combined_data$Source, levels = c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ Gold"))

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
    x = expression(beta[1]),  # Change x-axis label to "rho"
    #y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(values = c("ZZ Gold" = "black", 
                                "ZZ 2" = "#FF6F61",  # Soft coral red
                                "ZZ 5" = "#377EB8", # Medium slate blue
                                "ZZ 20" = "#32CD32")) +  # Lime green
  scale_linetype_manual(values = c("ZZ Gold" = "solid", 
                                   "ZZ 2" = "longdash", 
                                   "ZZ 5" = "twodash",
                                   "ZZ 20" = "dashed")) +
  theme(
    legend.position = c(1,1),  
    legend.key.width = unit(2, "cm"),
    axis.title.y=element_blank(),
    legend.justification = c(1, 1),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Light gray grid lines
    panel.grid.minor = element_line(color = "gray95", size = 0.25)  # Faint minor grid lines
  ) + 
xlim(3, 6)+
  ylim(0, 1.5)


# Display the plot
x11()
print(plot)





##### Posterior plot results - ZZ log sigma



samples_zigzag = sample_Gaussian_regression_IntractLik_gold[,9]

# True value
true = 0

# Custom labels for the columns
custom_labels <- c(
  "log(sigma)"
)

# Convert matrices into long-format data frames for ggplot
convert_to_long <- function(data, label) {
  as.data.frame(data) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Value") %>%
    mutate(Source = label)
}

data1 <- convert_to_long(samples_zigzag, "ZZ Gold")
data2 <- convert_to_long(sample_Gaussian_regression_IntractLik_m2[,9], "ZZ 2")
data4 <- convert_to_long(sample_Gaussian_regression_IntractLik_m5[,9], "ZZ 5")
data5 <- convert_to_long(sample_Gaussian_regression_IntractLik_m20[,9], "ZZ 20")

# Combine the data
combined_data <- bind_rows(data1, data2, data4, data5)

# Reorder Source factor to control the order in the legend
combined_data$Source <- factor(combined_data$Source, levels = c("ZZ 2", "ZZ 5", "ZZ 20", "ZZ Gold"))

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
    x = expression(log(sigma)),  # Change x-axis label to "rho"
    #y = "Density",
    color = NULL,  # Remove color legend title
    linetype = NULL  # Remove linetype legend title
  ) +
  scale_color_manual(values = c("ZZ Gold" = "black", 
                                "ZZ 2" = "#FF6F61",  # Soft coral red
                                "ZZ 5" = "#377EB8", # Medium slate blue
                                "ZZ 20" = "#32CD32")) +  # Lime green
  scale_linetype_manual(values = c("ZZ Gold" = "solid", 
                                   "ZZ 2" = "longdash", 
                                   "ZZ 5" = "twodash",
                                   "ZZ 20" = "dashed")) +
  theme(
    legend.position = c(1,1),  
    legend.key.width = unit(2, "cm"),
    axis.title.y=element_blank(),
    legend.justification = c(1, 1),
    panel.grid.major = element_line(color = "gray90", size = 0.5),  # Light gray grid lines
    panel.grid.minor = element_line(color = "gray95", size = 0.25)  # Faint minor grid lines
  ) + 
xlim(-0.5, 1.5)+
  ylim(0, 2.5)


# Display the plot
x11()
print(plot)












