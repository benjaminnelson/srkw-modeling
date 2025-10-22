## Same as previous version but with delay in the recovery start
# Load libraries
library("dplyr")
library("tidyr")
# library("brms")
library("ggplot2")
# library("loo")

## Clear everything ##
rm(list=ls()) 
gc() 

set.seed(1234)

setwd("C:/Users/nelso/OneDrive/BNelson/Consulting/Columbia Riverkeeper/srkw-modeling-columbia/salmon-modeling")

the_data<- read.csv("wcvi-2023.csv", header = TRUE)
catch_data<- read.csv("wcvi-2022-catch.csv", header = TRUE)
prior_data<- read.csv("wcvi-2023-prior.csv", header = TRUE)
source("helper_functions.r")

years<- the_data$Year
pops<- names(the_data %>% select(-Year, -AI.Total))

## Recovery scenarios
recovery_years<- c(10, 20, 30)
pct_increases<- c(0.128, 0.268, 0.409)

## Generate prior distribution for total abundance in 2023
log_devs<- log(prior_data$obs) - log(prior_data$pre)
devs<- prior_data$obs - prior_data$pre
sum_sq_devs<- sum((devs - mean(devs))^ 2)
prior_var<- (1/(nrow(prior_data) - 1)) * sum_sq_devs
prior_mean<- 1.02 # Pre-season WCVI (2023)
prior_sd<- sqrt(prior_var) # Std dev of pre/post residuals

## Plot WCVI time-series and prior for 2023
plot_cols<- viridis::viridis(length(pops))

par(mfrow=c(1,1), mar=c(4, 5, 2, 5) + 0.1)
plot(years, the_data$AI.Total, ann = FALSE, las = 1, type = "l", col = "white", ylim=c(0, 1.60), xlim = c(min(the_data$Year), max(the_data$Year) + 1))
for(i in 1:length(pops)){
  temp_stock<- names(the_data)[i]
  if(temp_stock == "URB") {
    lines(years, the_data[,(i+1)], col = alpha(plot_cols[i], 1.0), lwd = 2.5)
    text(2015, 1.55, font = 2, cex = 1.05, "Upriver Brights", col = alpha(plot_cols[i], 1.0))
  }
  if(temp_stock == "LYF") {
      lines(years, the_data[,(i+1)], col = alpha("grey5", 0.25), lwd = 2.5)
      # text(2015, 1.475, font = 2, cex = 1.05, "Lyons Ferry", col = alpha(plot_cols[i], 1.0))
    } 
  else {lines(years, the_data[,(i+1)], col = alpha("grey5", 0.25), lwd = 2.5)
    }
}
text(2015, 1.475, font = 2, cex = 1.05, "Total", col = alpha("black", 1.0))
# plot.window(xlim = range(the_data$Year), ylim = c(0, 2))
lines(the_data$Year, the_data$AI.Total, col="black", lwd = 2.5)
points(the_data$Year, the_data$AI.Total, col="black", pch = 16)
# points((max(years + 1)), prior_mean, pch = 9, col = "black")
# segments((max(years + 1)), (prior_mean + (1.96 * sqrt(prior_var))), (max(years + 1)), (prior_mean - (1.96 * sqrt(prior_var))), lwd = 1)
# abline(lm(the_data$AI.Total ~ the_data$Year), lwd = 2, lty = "dashed")
# axis(4, col.axis="black", las=1)
mtext(side = 2, line = 2.5, font = 2, "Index", outer = FALSE)

## Stan modeling block
model_data<- the_data %>% 
  select(-AI.Total) %>%
  pivot_longer(cols = -Year, names_to = "stock", values_to = "abundance") %>%
  mutate(
    abundance = if_else(abundance == 0, min(the_data[the_data != 0]), abundance),
    log_abundance = log(abundance),
    stock_id = dense_rank(stock)) %>%
  group_by(stock_id) %>%
  mutate(
    std_log_abundance = (log_abundance - mean(log_abundance)) / sd(log_abundance),
    log_mean_abundance = mean(log_abundance),
    log_sd_abundance = sd(log_abundance)) %>%
  ungroup() %>%
  arrange(stock_id, Year)
# Add sample number
model_data$sample_id<- seq(from = 1, to = nrow(model_data), by = 1)


# ## Populate model matrix with recovery scenarios
# model_array <- array(list(), dim = c(length(recovery_years), length(pct_increases)))
# for(i in 1:length(recovery_years)){
#   for(j in 1:length(pct_increases)){
#     print(c(i, j))
#     model_array[[i, j]] <- ar_model(mean_increase = pct_increases[j], increase_years = recovery_years[i], proj_years = 100);
#   }
# }
# 
# saveRDS(model_array, "C:/Users/nelso/Desktop/ar-models/ar-model-array.rds")

model_array<- readRDS("C:/Users/nelso/Desktop/ar-models/ar-model-array.rds")
# base_case<- ar_model(mean_increase = 0.0, increase_years = 10, proj_years = 100)
base_case<- readRDS("C:/Users/nelso/Desktop/ar-models/base-case.rds")

# Function to plot the URB and WCVI trends side-by-side
# par(mfrow = c(3, 2), oma = c(4, 3, 2, 1) + 0.1, mar = c(1, 1, 1, 1) + 0.1)
par(mfcol=c(2, 4), oma = c(4, 3, 2, 1) + 0.1, mar = c(1, 1, 0.1, 1) + 0.1)
# Plot mid low, med, high recovery scenarios for 20-year recovery
index_plot(ar_model = base_case)
mtext(side = 2, line = 2, font = 2, cex = 1.20, "WCVI abundance", outer = FALSE)
index_plot(ar_model = model_array[[2, 1]], write_csv = FALSE, pct_change = "13", yax = FALSE)
index_plot(ar_model = model_array[[2, 2]], write_csv = FALSE, pct_change = "27", yax = FALSE)
index_plot(ar_model = model_array[[2, 3]], write_csv = FALSE, pct_change = "41", yax = FALSE)
mtext(side = 1, line = 2, font = 2, cex = 1.20, "Year", outer = TRUE)
mtext(side = 3, line = 0.1, font = 2, cex = 1.20, "20-Year Recovery", outer = TRUE)
