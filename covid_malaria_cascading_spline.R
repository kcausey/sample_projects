#-------------------Header------------------------------------------------
# Author: Kate Causey
# Date: 3/4/2021
# Purpose: Cascade Spline Analysis for COVID-related disruptions to ACT and ITN Delivery
#          Adapted model for external partners at Malaria Atlas Project
#          Cascading random splines of health service disruption versus human mobility
#
#***************************************************************************

#------------------SET-UP--------------------------------------------------

# clear memory
rm(list=ls())

start_time <- Sys.time()

# load packages, install if missing

packages <- c("data.table","magrittr","ggplot2", "msm", 
              "dplyr", "gridExtra", "ggpubr", "wCorr", "dplyr", "matrixStats")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

# custom IHME package for MRBRT tool. See download instructions via email
library(mrtoolr)

# Directories and Inputs -------------------------------------------------------------

# path of the repository with code and input data
working_dir <- "{FILEPATH}" 

input_dir <- file.path(working_dir, "input_data")
output_dir <- file.path(working_dir, "output")

# read in location information and codes
hierarchy <- fread(file.path(input_dir, "location_metadata.csv"))
# subset to national-level locations
countries <- hierarchy[level == 3, location_id]

# read in predictor variable: mobility
mob <- fread(file.path(input_dir, "mobility.csv"))
mob[, date := as.Date(date)]

# calculate average mobility for PREMISE survey period
mob_premise <- mob[date>="2020-03-01" & date<"2020-07-01", .(mob_avg = mean(mob_avg)), by = "location_id"]

# Run MRBRT cascading spline model ---------------------------------------------------------------

pdf(file.path(output_dir, "spline_fits.pdf"), width = 17, height = 11)

# specify outcomes
ratios <- c("itn_ratio", "act_ratio")

for(ratio in ratios){
  
  model_dir <- file.path(output_dir, ratio)
  dir.create(model_dir, recursive = T)
    
  dt <- fread(paste0(input_dir, "/", ratio, ".csv"))
  
  # Merge on mobility
  dt <- merge(dt, mob_premise, all.x = T)
  
  # Global Model ------------------------------------------------------------
  
  dat_1 <- MRData()
  dat_1$load_df(
    data = dt,  col_obs = "log_ratio", col_obs_se = "log_ratio_se",
    col_covs = list("mob_avg"), col_study_id = "location_id"
  )
  
  
  model <- MRBRT(data = dat_1,
                 cov_models = list(LinearCovModel("mob_avg",
                                                  use_spline = TRUE,
                                                  use_re = FALSE, # no random effects
                                                  prior_spline_monotonicity = "decreasing", # disruption decreases as mobility disruption increases
                                                  spline_knots_type = "domain", # define evenly spaced knots
                                                  spline_knots = array(c(0, 0.2, 0.4, 1)), # 20%, 40% mobility disruption
                                                  spline_r_linear = TRUE, # right linear tail
                                                  spline_l_linear = FALSE, # free left tail
                                                  spline_degree = 2L, # cubic spline
                                                  prior_spline_maxder_gaussian = rbind(c(0, 0, 0), c(Inf, Inf, 0.1))))) # this prior suggests that the right linear tail should be flat (slope zero) unless the data suggests otherwise
  
  data_temp <- MRData(covs = list(mob_avg = array(seq(0,1,0.01)))) #set the range of the mobility variable to confirm knot placement
  
  model$attach_data(data_temp) 
  
  model$fit_model()
  
  #Predict Spline
  
  global_pred <- expand.grid(location_id = unique(dt$location_id),
                             mob_avg = seq(0,1,0.01)) %>% as.data.table()
  
  dat_pred <- MRData()
  
  dat_pred$load_df(
    data = global_pred,
    col_covs=list("mob_avg")
  )
  
  global_pred$pred <- model$predict(dat_pred) %>% exp()
  
  global_pred <- merge(global_pred, hierarchy[,.(location_id, location_name)])
  global_pred <- merge(global_pred, 
                       dt[, .(location_id, mob_avg_obs = mob_avg, obs = exp(log_ratio), lower, upper, weight)],
                       allow.cartesian = T)
  
  # plot fit against data, global model
  gg <- ggplot(global_pred, aes(x = mob_avg_obs, y = obs))+
    geom_point(aes(size = weight))+
    geom_line(aes(x= mob_avg, y = pred))+
    theme_bw()+
    labs(title = paste0("Global Model: ", ratio),
         x = "Mobility Disruption",
         y = "Health Services Disruption")
  
  print(gg)
  
  # Location specific models (stage 2) --------------------------------------
  
  # get stage 1 estimates as prior for stage 2
  samples_model <- mrtoolr::core$other_sampling$sample_simple_lme_beta(sample_size = 1000L, model = model)
  sds_model <- apply(samples_model, 2, sd)
  
  model_locs <- unique(dt$location_id)
  
  # set up datasets to predict out
  loc_pred <- data.table()
  mob_pred <- data.table()
  mob_pred_global <- mob[!(location_id %in% model_locs)][order(mob_avg)]
  
  # location counter
  l <- 1
  
  for(loc_id in model_locs){
    
    message(loc_id)
    
    # controls the strength of the global prior
    theta <- 3

    dt_loc <- dt[location_id == loc_id]
    
    dat_2 <- MRData()
    dat_2$load_df(
      data = dt_loc,  col_obs = "log_ratio", col_obs_se = "log_ratio_se",
      col_covs = list("mob_avg"), col_study_id = "location_id"
    )
    
    # use stage 1 as a prior for stage 2
    mod2_beta_prior <- rbind(model$beta_soln, sds_model * theta)
    
    # this prior suggests that the right linear tail should be flat (slope zero) unless the data suggests otherwise
    mod2_max_der_prior <- rbind(rep(0, ncol(mod2_beta_prior) - 1),
                                c(rep(Inf, ncol(mod2_beta_prior) -2), 0.1)) # sets 
    
    # location-specific spline
    model_2 <- MRBRT(data = dat_2,
                     cov_models = list(LinearCovModel("mob_avg",
                                                      use_spline = TRUE,
                                                      use_re = FALSE,
                                                      prior_spline_monotonicity = "decreasing",
                                                      spline_knots_type = "domain",
                                                      spline_knots = array(c(0, 0.2, 0.4, 1)),
                                                      spline_r_linear = TRUE,
                                                      spline_l_linear = FALSE, 
                                                      spline_degree = 2L,
                                                      prior_beta_gaussian = mod2_beta_prior),
                                                      prior_spline_maxder_gaussian = mod2_max_der_prior)))
    
    model_2$attach_data(data_temp)
    
    model_2$fit_model()
    
    #Predict spline for plotting
    
    dt_pred2 <- expand.grid(location_id = loc_id,
                            mob_avg = seq(0,1,0.01)) %>% as.data.table()
    
    dat_pred2 <- MRData()
    
    dat_pred2$load_df(
      data = dt_pred2,
      col_covs=list("mob_avg")
    )
    
    dt_pred2$pred <- model_2$predict(dat_pred2) %>% exp()
    
    dt_pred2 <- merge(dt_pred2, hierarchy[,.(location_id, location_name)])
    dt_pred2 <- merge(dt_pred2, dt_loc[, .(location_id, mob_avg_obs = mob_avg, obs = exp(log_ratio), 
                                           lower, upper)],
                      allow.cartesian = T)
    
    loc_pred <- rbind(loc_pred, dt_pred2)
    
    #Predict out based on observed mobility for location
    dt_pred3 <- mob[location_id == loc_id][order(mob_avg)]
    
    dat_pred3 <- MRData()
    
    dat_pred3$load_df(
      data = dt_pred3,
      col_covs=list("mob_avg")
    )
    
    dt_pred3$pred <- model_2$predict(dat_pred3) %>% exp()
    
    mob_pred <- rbind(mob_pred, dt_pred3)
    
    l <- l+1
    
  }
  
  # plot global and stage 2 for each location on 1 figure
  
  min_plot_val <- min(loc_pred$pred, global_pred$pred, loc_pred$obs)
  max_plot_val <- max(loc_pred$pred, global_pred$pred, loc_pred$obs)
  
  loc_pred[, plot_lower := lower]
  loc_pred[, plot_upper := upper]
  loc_pred[lower < min_plot_val, plot_lower := min_plot_val]
  loc_pred[upper > max_plot_val, plot_upper := max_plot_val]
  
  gg <- ggplot(loc_pred)+
    geom_hline(yintercept = 1, color = "red", size = 0.5)+
    geom_line(aes(x=mob_avg, y=pred))+
    geom_line(data = global_pred[location_id %in% c(model_locs)], aes(x=mob_avg, y=pred), color = "blue", alpha = 0.5)+
    geom_point(aes(x=mob_avg_obs, y=obs))+
    geom_errorbar(aes(x=mob_avg_obs, ymin = plot_lower, ymax = plot_upper), alpha = 0.5, size = 0.5)+
    facet_wrap(~location_name)+
    labs(title = paste0("Location Spline fit: ", ratio, " v mobility"),
         color = "Premise")+
    theme_bw()+
    scale_y_continuous(limits = c(min_plot_val, max_plot_val))
  
  print(gg)
  
  # Predict out annual 2020 from Global model---------------------
  
  dat_pred4 <- MRData()
  
  dat_pred4$load_df(
    data = mob_pred_global,
    col_covs=list("mob_avg")
  )
  
  mob_pred_global$pred <- model$predict(dat_pred4) %>% exp()
  
  # Daily results
  
  mob_pred <- rbind(mob_pred, mob_pred_global)
    
  mob_pred[mob_avg <= 0, "pred" := 1]
    
  mob_pred[, year_id := year(date)]
  
  write.csv(mob_pred, paste0(model_dir, "/mrbrt_", ratio, "_results_daily.csv"), row.names = F)
  
  # Annual results
    
  pred_annual <- mob_pred[, .(pred = mean(pred), mob_avg = mean(mob_avg)), by = c("location_id", "year_id")]
  
  pred_annual <- merge(pred_annual, hierarchy[, .(location_id, ihme_loc_id, location_name)])
  
  write.csv(pred_annual, paste0(model_dir, "/mrbrt_", ratio, "_results_annual.csv"),row.names = F)

  
}

dev.off()

# report amount of time for program
end_time <- Sys.time()
end_time - start_time
