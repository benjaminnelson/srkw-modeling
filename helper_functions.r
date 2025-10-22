z_back<- function(x, sd, mean) exp(x * stock_sd + stock_mean)

# Function for plotting indices
index_plot<- function(ar_model, n_proj = 100, spaghetti_sims = 10, wcvi_sims = 1000, write_csv = FALSE, pct_change = NULL, yax = TRUE){
  # Define model posteriors
  ppd_seq<- seq(from = 1, to = (length(pops) * n_proj), by = n_proj)
  posts<- as.data.frame(ar_model)
  etas<- matrixStats::colQuantiles(
    as.matrix(posts[ ,grepl("eta", colnames(posts))]), probs = c(0.50, 0.025, 0.975))[-(1:length(pops)), ]
  total_etas<- matrixStats::colQuantiles(
    as.matrix(posts[ ,grepl("total_abs_eta", colnames(posts))]), probs = c(0.50, 0.025, 0.975))
  ppds<- matrixStats::colQuantiles(
    as.matrix(posts[ ,grepl("eta_rep_proj", colnames(posts))]), probs = c(0.50, 0.025, 0.975, 0.05, 0.95, 0.25, 0.75))
  raw_ppds<- as.matrix(posts[ ,grepl("eta_rep_proj", colnames(posts))])
  means<- matrixStats::colQuantiles(
    as.matrix(posts[ ,grepl("mean_value", colnames(posts))]), probs = c(0.50, 0.025, 0.975))
  # print(means)
  raw_log_indices<- matrixStats::colQuantiles(
    as.matrix(posts[ ,grepl("raw_indices", colnames(posts))]), probs = c(0.50, 0.025, 0.975))
  total_log_indices<- matrixStats::colQuantiles(
    as.matrix(posts[ ,grepl("total_indices", colnames(posts))]), probs = c(0.50, 0.025, 0.975))
  raw_total_log_indicies<- posts[ , grepl("total_indices", colnames(posts))]
  if(write_csv != FALSE) {
    file_name<- paste0(as.character(pct_change), "_over_", n_proj, "_years.csv", sep = "")
    write.csv(raw_total_log_indicies, file_name, row.names = FALSE)
  }
  
  # Plot URB
  temp_stock<- model_data[model_data$stock == "URB", ]
  y_max<- 1.5 * max(ppds[ppd_seq[31]:(ppd_seq[31] + n_proj - 1), 3])
  # y_min<- 1.5 * min(temp_stock$std_log_abundance)
  # Observed
  plot(temp_stock$Year, temp_stock$std_log_abundance, ann = FALSE, pch = 16, las = 1, col = "black", ylim = c(-4, 8), yaxt = ifelse(yax == TRUE, "s", "n"), xaxt = "n", xlim = c(min(temp_stock$Year), (max(temp_stock$Year) + n_proj)))
  axis(1, labels = FALSE)
  axis(2, labels = FALSE)
  abline(h = 0, col = "grey50", lty = "dashed")
  lines(temp_stock$Year, temp_stock$std_log_abundance, col = "black")
  # Predicted
  # Spaghetti 
  ppd_samples<- sample(1:nrow(raw_ppds), spaghetti_sims, FALSE)
  for(j in 1:spaghetti_sims) {lines((max(years) + 1):(max(years) + n_proj), raw_ppds[ppd_samples[j], ppd_seq[31]:(ppd_seq[31] + n_proj - 1)], col = alpha("purple", 0.75), lwd = 0.50)}
  lines((max(temp_stock$Year + 1)):(max(temp_stock$Year) + n_proj), ppds[ppd_seq[31]:(ppd_seq[31] + n_proj - 1), 1], col = "purple", lwd = 2)
  if(yax == TRUE){mtext(side = 2, line = 2, font = 2, cex = 1.20, "URB abundance", outer = FALSE)}
  
  # Plot WCVI
  total_ppd_samples<- sample(1:nrow(raw_total_log_indicies), wcvi_sims, FALSE)
  plot(the_data$Year, the_data$AI.Total, ann = FALSE, las = 1, cex = 1.20, pch = 16, ylim = c(0, 10), yaxt = ifelse(yax == TRUE, "s", "n"), xlim = c(min(the_data$Year), (max(the_data$Year) + n_proj)))
  axis(2, labels = FALSE)
  lines(the_data$Year, the_data$AI.Total)
  for(j in 1:wcvi_sims) {lines((max(years) + 1):(max(years) + n_proj), raw_total_log_indicies[total_ppd_samples[j], ], col = alpha("blue", 0.02), lwd = 0.50)}
  lines((max(years) + 1):(max(years) + n_proj), total_log_indices[,1], lwd = 0.5, col = alpha("blue", 0.0005))
}


ar_model<- function(mean_increase = 0.0, increase_years = 10, proj_years = 50){
  ## Stan Modeling settings
  # Data/pars definitions/assignments
  S = length(unique(model_data$stock)) # Number of stocks
  T = max(model_data$Year - min(model_data$Year) + 1) # Maximum time-step in dataset
  x = matrix(model_data$std_log_abundance, nrow = T, ncol = S)
  n_proj = proj_years
  mean_increase = mean_increase
  increase_years = increase_years
  increase_delay = 5
  increase_stock_id = 31
  stock_mean<- unique(model_data$log_mean_abundance)
  stock_sd<- unique(model_data$log_sd_abundance)
  
  # Create a data list for Stan
  data_list <- list(
    x = x, 
    T = T, 
    S = S,
    stock_mean = stock_mean,
    stock_sd = stock_sd,
    n_proj = n_proj, 
    increase_delay = increase_delay,
    increase_stock_id = increase_stock_id,
    mean_increase = mean_increase, 
    increase_years = increase_years,
    prior_mean = prior_mean,
    prior_sd = prior_sd
  )
  
  # Model fitting
  chains = 3
  iter = 1000
  thin = 1
  
  ar <- rstan::stan("ar1-multi-prior-delay.stan",
                    data = data_list,
                    chains = chains,
                    warmup = iter / 2,
                    iter = iter,
                    thin = thin,
                    seed = 1234,
                    cores = chains)
  return(ar)
}

