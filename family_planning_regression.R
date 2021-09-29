
#-------------------Header------------------------------------------------
# Author: Kate Causey
# Date: 9/1/21
# Purpose: Stage 1 analysis for drivers of change of family planning 
#          Global, country-level analysis for Exemplars in Global Health
#          Outcomes: mCPR and demand satisfied (country-year)
#          Testing various dependent variables, effect size in linear 
#               regression, and relative importance via Shapley decomposition
#          
#***************************************************************************

#------------------SET-UP--------------------------------------------------

# clear memory
rm(list=ls())

# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "{FILEPATH}"
  h_root <- "{FILEPATH}"
  central_lib <- "{FILEPATH}"
  } else {
  j_root <- "{FILEPATH}"
  h_root <- "{FILEPATH}"
  central_lib <- "{FILEPATH}"
  }


# load packages, install if missing

lib.loc <- paste0(h_root,"R/",R.Version()$platform,"/",R.Version()$major,".",R.Version()$minor)
dir.create(lib.loc,recursive=T, showWarnings = F)
.libPaths(c(lib.loc,.libPaths()))

packages <- c("data.table","magrittr","ggplot2", "relaimpo", "openxlsx", "miceadds", "estimatr", "matrixStats", "ggrepel", "jtools")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

# Set-up -------------------------------------------------------------

# directories
parent_dir <- file.path(j_root, "{FILEPATH}")

out_dir <- file.path(parent_dir, "phase1", Sys.Date())
dir.create(out_dir, recursive = T)

dir.create(file.path(out_dir, "model_results"))

# specify years and ages to be included in analysis
years <- 2000:2019
ages_detailed <- 8:14
ages <- c(ages_detailed, 24)

# map of age ids to names
age_map <- data.table(age_id = ages, 
                      age_label = c("15 to 19",
                                    "20 to 14", 
                                    "25 to 29",
                                    "30 to 34",
                                    "35 to 39",
                                    "40 to 44",
                                    "45 to 49",
                                    "15 to 49"))

# central functions to pull internal IHME id codebook
source(file.path(central_lib, "current/r/get_location_metadata.R"))
locations <- get_location_metadata(35, gbd_round_id = 7, decomp_step = "iterative")

# limit to countries classified as LMIC in 2000
lmic_dt <- fread(file.path(j_root, "temp/kcausey/gates_ventures/data_sharing", "country_classification.csv"))
lmic_iso3 <- lmic_dt[income_grp_2000 %in% c("Low income", "Lower middle income"), iso]
countries <- locations[ihme_loc_id %in% lmic_iso3, location_id]

# pull data source counts for sensitivity analysis
data_sparse <- fread(file.path(parent_dir, "raw_data", "contra_source_counts_2000_2017.csv"))[ihme_loc_id %in% lmic_iso3]

data_rich_3 <- data_sparse[any_contra >=3, ihme_loc_id]
data_rich_4 <- data_sparse[any_contra >=4, ihme_loc_id]

# define function to run analysis
# provide: formula, dataset, name for plotting and saving, and whether covariates should be standardized
model_plot_save <- function(form, model_dt, model_title, std = F){
  
  start_time <- Sys.time()
  
  # standardize coefficents by standardizing the input covariates into z-scores
  if(std){
    model_title <- paste0(model_title, " standardized coefficients")
    
    for(this.var in setdiff(all.vars(form), "ihme_loc_id")){
      model_dt[, eval(this.var) := (get(this.var) - mean(get(this.var), na.rm = T))/sd(get(this.var), na.rm = T)]
    }
  }
  
  # account for clustered nature of dataset by country
  des <- svydesign(id=~location_id, data=model_dt, variables = model_dt, weights = NULL)
  
  # run linear regression
  mod_1 <- svyglm(form,
                    design = des, 
                    na.action = "na.omit")
    
    pred_dt <- na.omit(model_dt, cols = c("location_id", all.vars(form)))
    pred_dt[, pred := predict(mod_1, newdata = pred_dt) %>% as.numeric()]
    
    # highlight outliers    
    pred_dt[, resid := get(all.vars(form)[[1]]) - pred]    
    pred_dt[year_id %in% c(max(year_id)) & abs(resid) > 0.2, lab := ihme_loc_id]
    
    # plot to show observed versus predicted, highlighting outliers
    gg <- ggplot(pred_dt, aes(y = get(all.vars(form)[[1]]), x = pred, color = year_id))+
      geom_abline(aes(slope = 1, intercept = 0))+
      geom_point()+
      geom_label_repel(aes(label = lab))+
      labs(title =model_title,
           y = "Observed",
           x = "Predicted",
           color = "Year")+
      theme_bw()

    print(gg)
    
    # print positive and negative outliers of each model in absolute space
    print(head(pred_dt[year_id == max(year_id), .(location_name, resid)][order(-resid)], 8))
    print(tail(pred_dt[year_id == max(year_id), .(location_name, resid)][order(-resid)], 3))
    print(pred_dt[year_id == max(year_id), .N])
    
    # print positive and negative outliers of each model accounting for change over time
    chg_resid_dt <- merge(pred_dt[year_id == max(year_id), .(location_name, resid)], pred_dt[year_id == min(year_id), .(location_name, resid)],
                          by = "location_name",
                          suffixes = c("2", "1"))    
    chg_resid_dt[, chg_resid := resid2 - resid1]
    
    print(head(chg_resid_dt[, .(location_name, resid1, resid2, chg_resid)][order(-chg_resid)], 8))
    print(tail(chg_resid_dt[, .(location_name, resid1, resid2, chg_resid)][order(-chg_resid)], 3))
    print(chg_resid_dt[, .N])
    
    # use relaimpo package to calculate relative importance of each variable to the variance explained by the model
    bt <- boot.relimp(mod_1, b=100, rela = T)
    draws <- booteval.relimp(bt)$lmg.boot
    
      # calculate dataset of model coefficients, confidence intervals, and relative importance
      save_dt <- data.table(var = names(coef(mod_1)),
                            beta = coef(mod_1),
                            lower = confint(mod_1)[,1],
                            upper = confint(mod_1)[,2],
                            rel_importance_mean = c(NA, colMeans(draws)),
                            rel_importance_lower = c(NA, colQuantiles(draws, probs = 0.025)),
                            rel_importance_upper = c(NA, colQuantiles(draws, probs = 0.975)))
       
    write.csv(save_dt, file.path(out_dir, "model_results", paste0(model_title,".csv")), row.names = T)
    
    save_dt[, sig := "Not significant"]
    save_dt[upper < 0, sig := "Significantly negative"]
    save_dt[lower > 0, sig := "Significantly positive"]
    
    # generate plot of model coefficients and variance explained
    plot_dt <- melt.data.table(save_dt, id.vars = c("var", "sig"))
    plot_dt[, measure := "Relative Importance"]
    plot_dt[variable %in% c("beta", "lower", "upper"), measure := "Coefficient"]
    plot_dt[, metric := "mean"]
    plot_dt[grepl("lower", variable), metric := "lower"]
    plot_dt[grepl("upper", variable), metric := "upper"]
    
    plot_dt <- dcast.data.table(plot_dt, var + measure + sig ~ metric, value.var = "value")
    
    # set reasonable limits for plot
    cutoff <- 1
    plot_dt[lower > cutoff, lower := Inf]
    plot_dt[mean > cutoff, mean := Inf]
    plot_dt[upper > cutoff, upper := Inf]
    plot_dt[lower < (-1)*cutoff, lower := -Inf]
    plot_dt[mean  < (-1)*cutoff, mean := -Inf]
    plot_dt[upper  < (-1)*cutoff, upper := -Inf]
    
    rsquared <- round(attr(summ(mod_1), "rsq"), 2)*100
        
    gg <- ggplot(plot_dt[], aes(y = var, x = mean, fill = sig))+
      geom_col()+
      geom_segment(aes(yend = var, x = lower, xend = upper))+
      facet_wrap(~measure, scales = "free_x")+
      labs(title = model_title,
           subtitle = paste0("R-squared = ",rsquared, "%"),
           x = "",
           y = "",
           fill = "")+
      scale_fill_manual(values = c("Not significant" = "#adadad",
                                   "Significantly negative" = "#f2766d",
                                   "Significantly positive" = "#6da9e8"))+
      theme_bw()+
      theme(legend.position = "bottom")
    
    print(gg)
    
    print(paste0("Done with ", model_title, ", ", Sys.time()-start_time))

  
}


# Run model ---------------------------------------------------------------

# read in prepped dataset with dependent and independent variables
dt <- fread(file.path(out_dir, "square_dataset.csv"))

# shift conflict away from zero by the 1st percentile
conflict_cutpoint <- dt[conflict != 0, quantile(conflict, 0.01)]
conflict_lag_cutpoint <- dt[conflict_lag != 0, quantile(conflict_lag, 0.01)]
dah_cutpoint <- dt[fp_spend_dah != 0, quantile(fp_spend_dah, 0.01)]
dah_pc_cutpoint <- dt[fp_spend_dah_pc != 0, quantile(fp_spend_dah_pc, 0.01)]

dt[, conflict := conflict + conflict_cutpoint]
dt[, conflict_lag := conflict_lag + conflict_lag_cutpoint]
dt[, fp_spend_dah := fp_spend_dah + dah_cutpoint]
dt[, fp_spend_dah_pc := fp_spend_dah_pc + dah_pc_cutpoint]

# plot distributions of independent variables to evaluate normality assumption
pdf(file.path(out_dir, "predictor_distributions.pdf"), width = 8, height = 6)

for(var in c("demand_satisfied", "mcpr", "prop_married_15_24", "urban", "ldi", "edu_y_pc_15_49_f", 
             "fp_spend_pub_pc", "abortion_legality", "muslim_prop", "mean_repro_age", "conflict",
             "fp_spend_dah_pc", "fp_spend_gov_pc", "fp_spend_oop_pc", "fp_spend_all_pc")){
  
  hist(dt[age_group_id==24, get(var)], main = var)
  hist(log(dt[age_group_id==24, get(var)]), main = paste0("log(", var, ")"))
}


dev.off()

# take log of variables with approximate log-normal distribution after visual examination
dt[, log_ldi := log(ldi)]
dt[, log_u5mr := log(u5mr)]
dt[, log_hrh := log(hrh)]
dt[, log_fp_spend_pub_pc := log(fp_spend_pub_pc)]
dt[, log_conflict := log(conflict)]
dt[, log_conflict_lag := log(conflict_lag)]
dt[, log_fp_spend_dah_pc := log(fp_spend_dah_pc)]
dt[, log_fp_spend_gov_pc := log(fp_spend_gov_pc)]
dt[, log_fp_spend_oop_pc := log(fp_spend_oop_pc)]
dt[, log_fp_spend_all_pc := log(fp_spend_all_pc)]


# Model iterations --------------------------------------------------------

pdf(file.path(out_dir, "model_results_mcpr.pdf"), width = 8, height = 6)

model_plot_save(form = mcpr ~  prop_married_15_24 + urban + log_ldi + edu_y_pc_15_49_f + 
                       muslim_prop + mean_repro_age + log_fp_spend_pub_pc + abortion_legality,
                model_dt = dt[age_group_id == 24],
                model_title = "All age MCPR base model")

model_plot_save(form = mcpr ~  prop_married_15_24 + urban + log_ldi + edu_y_pc_15_49_f + 
                       muslim_prop + mean_repro_age + log_fp_spend_pub_pc + abortion_legality,
                model_dt = dt[age_group_id == 24],
                model_title = "All age MCPR base model",
                std = T)

model_plot_save(form = mcpr ~  prop_married_15_24 + urban + log_ldi + edu_y_pc_15_49_f + 
                       muslim_prop + mean_repro_age + log_fp_spend_dah_pc + log_fp_spend_gov_pc + 
                       log_fp_spend_oop_pc + abortion_legality,
                model_dt = dt[age_group_id == 24],
                model_title = "All age MCPR spend disaggregated")


model_plot_save(form = mcpr ~  prop_married_15_24 + urban + log_ldi + edu_y_pc_15_49_f + 
                       muslim_prop + mean_repro_age + log_fp_spend_pub_pc + abortion_legality +
                       log_conflict_lag,
                model_dt = dt[age_group_id == 24],
                model_title = "All age MCPR, add conflict 10-year average")

model_plot_save(form = mcpr ~  prop_married_15_24 + urban + log_ldi + edu_y_pc_15_49_f + 
                       muslim_prop + mean_repro_age + log_fp_spend_pub_pc + abortion_legality +
                       log_conflict_lag,
                model_dt = dt[age_group_id == 24 & ihme_loc_id %in% data_rich_3],
                model_title = "All age MCPR, data rich (3+ country-years)")

model_plot_save(form = mcpr ~  prop_married_15_24 + urban + log_ldi + edu_y_pc_15_49_f + 
                       muslim_prop + mean_repro_age + log_fp_spend_pub_pc + abortion_legality +
                       log_conflict_lag,
                model_dt = dt[age_group_id == 24 & ihme_loc_id %in% data_rich_4],
                model_title = "All age MCPR, data rich (4+ country-years)")

dev.off()


pdf(file.path(out_dir, "model_results_demand.pdf"), width = 8, height = 6)

model_plot_save(form = demand_satisfied ~  prop_married_15_24 + urban + log_ldi + edu_y_pc_15_49_f + 
                       muslim_prop + mean_repro_age + log_fp_spend_pub_pc + abortion_legality +
                       log_conflict_lag,
                model_dt = dt[age_group_id == 24],
                model_title = "All age demand_satisfied, add conflict 10-year average")

dev.off()
