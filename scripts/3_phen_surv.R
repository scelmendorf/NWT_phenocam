#code to fit discrete time logistic survival model to nwt phenocam data
# cross-validate; generate predictions
#sce 8/12/2021

#todo
#sce check that the plots show all retained terms in the BEST mod
# make comparisons of soil vs air in terms of RMSE
# best predictions

#make peak and eos models start at start, peak

#consider doing lasso regression instead


####setup####
library (tidyverse)
library(GGally)
library (cowplot)
library (gridExtra)
library (grid)


####Read in data and munge a bit to proper format####
# met data comes from script 1a_infill_met.R
met_data <- readRDS("data/all_met_infilled_seq_6_no_freeze_fill.rds") %>%
  # mutate(soil30_moist_is_infilled = ifelse(is.na(soilmoisture_b_30cm_avg_ave_infilled), FALSE, TRUE)) %>%
  # mutate(soilmoisture_b_30cm_avg_ave = ifelse(is.na(soilmoisture_b_30cm_avg_ave), soilmoisture_b_30cm_avg_ave_infilled,
  #                                             soilmoisture_b_30cm_avg_ave
  # )) %>%
  # mutate(soil5_temp_is_infilled = ifelse(is.na(soiltemp_5cm_avg_ave_infilled), FALSE, TRUE)) %>%
  # mutate(soiltemp_5cm_avg_ave = ifelse(is.na(soiltemp_5cm_avg_ave), soiltemp_5cm_avg_ave_infilled,
  #                                      soiltemp_5cm_avg_ave
  # )) %>%
  # mutate(aT_is_infilled = ifelse(is.na(aT_infilled), FALSE, TRUE)) %>%
  # mutate(airtemp_avg_ave = ifelse(is.na(airtemp_avg_ave), aT_infilled,
  #                                 airtemp_avg_ave)) %>%
  # #one more that did not get filled
  # mutate(airtemp_avg_ave = ifelse(sensornode == 'sn_14' & ymd =='2020-05-15',
  #                                 -0.512236, airtemp_avg_ave))%>%
  mutate(
    year = lubridate::year(ymd),
    yday = lubridate::yday(ymd)
  ) 

# all_phen comes from script 0_prep_data.R
# need this for snowmelt dates
all_phen <- read.csv("data/phen_clim_all.csv")
all_phen <- all_phen %>%
  select(year, sensornode, snowmelt_doy_infilled) %>%
  distinct() %>%
  mutate(sensornode = paste0("sn_", sensornode))

#read the phenometrics, made from 1_analyze_saddle_phenocams
sumry = read.csv("data_deriv/phenometrics.csv")

# add sf y or n, add phenophase y or n
met_data <- met_data %>%
  #add in snow
  left_join(., all_phen) %>%
  rowwise() %>%
  mutate(snow = dplyr::if_else(yday > snowmelt_doy_infilled, 1, 0)) %>%
  ungroup() %>%
  #sce add this in to see if we can get units on same scale
  mutate(soilmoisture_b_5cm_avg_ave = 100*soilmoisture_b_5cm_avg_ave)

#sce needs to figure out what happened w sn1,17,21 in 2018?
#dtr used for snowmelt where not observed, resulting in a bit smaller 
# sample for 2018.

((ggplot(met_data %>%
           filter(year>2017) %>%
           mutate(year = factor(lubridate::year(ymd)),
                            doy = lubridate::yday(ymd)), aes(x = doy))+
    geom_line(aes(y = soilmoisture_b_5cm_avg_ave, group = year, color = year),
               size = 0.5, alpha = 1))+
    ggtitle ("soil moisture 5cm for each node in each year \n vertical lines show snowmelt date")+
    geom_vline(
               aes(xintercept=snowmelt_doy_infilled, color = year))+
    facet_wrap(~sensornode)) %>%
  ggsave(file='plots/raw_data_plots/soil_moist5_seasonality_infilled.jpg', ., width = 8, height = 5)

((ggplot(met_data %>%
           filter(year>2017) %>%
           mutate(year = factor(lubridate::year(ymd)),
                  doy = lubridate::yday(ymd)), aes(x = doy))+
    geom_line(aes(y = soiltemp_5cm_avg_ave, group = year, color = year),
              size = 0.5, alpha = 1))+
    ggtitle ("soil moisture 5cm for each node in each year \n vertical lines show snowmelt date")+
    geom_vline(
      aes(xintercept=snowmelt_doy_infilled, color = year))+
    facet_wrap(~sensornode)) %>%
  ggsave(file='plots/raw_data_plots/soil_temp5_seasonality_infilled.jpg', ., width = 8, height = 5)




## define functions to process the sumry data into
# a consistent format per phenophase
make_phen_data <- function(sumry_data, phenophase) {
  phen_data <- sumry_data %>%
    select(sensornode, !!phenophase, year) %>%
    rename(doy = !!phenophase) %>%
    mutate(sensornode = paste0("sn_", sensornode)) %>%
    mutate(year = as.numeric(year)) %>%
    filter(year != 2017) %>%
    # this one has no snowmelt day in 2018
    filter(!(year == 2018 & sensornode == "sn_17"))
  return(phen_data)
}

phen_data_sos = make_phen_data(sumry, "sos") %>%
  rename(doy_event = doy)

phen_data_pop = make_phen_data(sumry, "pop") %>%
  rename(doy_event = doy)

phen_data_eos = make_phen_data(sumry, "eos") %>%
  rename(doy_event = doy)

#add drivers; truncate data after
# event which cannot be a driver (future cannot determine past)
phen_data_sos = phen_data_sos %>%
  left_join(., (met_data %>%
                  select(airtemp_avg_ave,
                         soiltemp_5cm_avg_ave, 
                         soilmoisture_b_5cm_avg_ave,
                         sensornode, ymd, yday,year, snow,
                         snowmelt_doy_infilled))) %>%
  filter(yday<=doy_event)

phen_data_pop = phen_data_pop %>%
  left_join(., (met_data %>%
                  select(airtemp_avg_ave,
                         soiltemp_5cm_avg_ave, 
                         soilmoisture_b_5cm_avg_ave,
                         sensornode, ymd, yday,year, snow,
                         snowmelt_doy_infilled))) %>%
  filter(yday<=doy_event)


phen_data_eos = phen_data_eos %>%
  left_join(., (met_data %>%
                  select(airtemp_avg_ave,
                         soiltemp_5cm_avg_ave, 
                         soilmoisture_b_5cm_avg_ave,
                         sensornode, ymd, yday,year, snow,
                         snowmelt_doy_infilled))) %>%
  filter(yday<=doy_event)

#check for complete predictor sets####
#make sure you have the right number of days ie aren't missing
# any driver days
check_complete_sos = phen_data_sos %>%
  group_by(sensornode, year) %>%
  summarize(nday = dplyr::n())

check_complete_sos  = phen_data_sos %>%
  select(sensornode, year, doy_event) %>%
  distinct() %>%
  left_join(., check_complete_sos) %>%
  filter(nday !=doy_event)

if (nrow(check_complete_sos) >0){
  phen_data = phen_data_sos %>%
    anti_join(., check_complete_sos %>%
                select(year, sensornode)%>%
                distinct())
}

check_complete_pop = phen_data_pop %>%
  group_by(sensornode, year) %>%
  summarize(nday = dplyr::n())

check_complete_pop  = phen_data_pop %>%
  select(sensornode, year, doy_event) %>%
  distinct() %>%
  left_join(., check_complete_pop) %>%
  filter(nday !=doy_event)

if (nrow(check_complete_pop) >0){
  phen_data = phen_data_pop %>%
    anti_join(., check_complete_pop %>%
                select(year, sensornode)%>%
                distinct())
}

check_complete_eos = phen_data_eos %>%
  group_by(sensornode, year) %>%
  summarize(nday = dplyr::n())

check_complete_eos  = phen_data_eos %>%
  select(sensornode, year, doy_event) %>%
  distinct() %>%
  left_join(., check_complete_eos) %>%
  filter(nday !=doy_event)

if (nrow(check_complete_eos) >0){
  phen_data = phen_data_eos %>%
    anti_join(., check_complete_eos %>%
                select(year, sensornode)%>%
                distinct())
}


#format long aka pseudo-observation format####
#Make into long format for time-to-event analysis
# adding in pseudo observations from before the event
#https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf

format_phen <- function(df){
  #df = data frame from make_phen_data
  # always days prior to snowmelt are remove; assumption is 0 prob of
  # event while after snow
  df = df %>%
    unite(subject, sensornode, year, remove = FALSE) %>%
    #note that this implies they respond to the climate forcing on the doy of
    # the observed event, maybe technically we should move it back 1d but
    # we also likely don't have that accurate phenodates anyhow
    mutate(status = ifelse(yday == doy_event, 1, 0),
           last_day = yday-1,
           days_since_snow = yday - snowmelt_doy_infilled) %>%
    mutate(days_since_snow = ifelse(days_since_snow<0, 0, days_since_snow),
           sensorid = as.numeric(as.factor(sensornode))) %>%
    arrange(sensornode, year, yday)%>%
    group_by(sensornode, year) %>%
    ##DDs based on 5 cm soil temp, threshold of 0
    arrange(sensornode, year, yday)%>%
    group_by(sensornode, year) %>%
    mutate(DD = ifelse(soiltemp_5cm_avg_ave>0, soiltemp_5cm_avg_ave, 0),
           GDD = cumsum(DD),
           DD_air = ifelse(airtemp_avg_ave>0, airtemp_avg_ave, 0),
           GDD_air = cumsum(DD_air),
           moist_7d = slider::slide_dbl(soilmoisture_b_5cm_avg_ave, mean, .before = 7, .after = 0),
           temp_7d_air = slider::slide_dbl(airtemp_avg_ave, mean, .before = 7, .after = 0),
           temp_7d = slider::slide_dbl(soiltemp_5cm_avg_ave, mean, .before = 7, .after = 0)) %>%
    ungroup() 
}

phen_data_sos = format_phen(phen_data_sos)
phen_data_pop = format_phen(phen_data_pop)
phen_data_eos = format_phen(phen_data_eos)

#include the met data after the event to viz preds and calc
# predicted transition dates which could be predicted to occur
# after the real one
phen_data_sos_pred = make_phen_data(sumry, "sos") %>%
  rename(doy_event = doy)%>%
  left_join(., (met_data %>%
                  select(airtemp_avg_ave,soiltemp_5cm_avg_ave, 
                         soilmoisture_b_5cm_avg_ave,
                         sensornode, ymd, yday,year, snow,
                         snowmelt_doy_infilled)))

phen_data_pop_pred = make_phen_data(sumry, "pop") %>%
  rename(doy_event = doy)%>%
  left_join(., (met_data %>%
                  select(airtemp_avg_ave,soiltemp_5cm_avg_ave, 
                         soilmoisture_b_5cm_avg_ave,
                         sensornode, ymd, yday,year, snow,
                         snowmelt_doy_infilled)))

phen_data_eos_pred = make_phen_data(sumry, "eos") %>%
  rename(doy_event = doy)%>%
  left_join(., (met_data %>%
                  select(airtemp_avg_ave,soiltemp_5cm_avg_ave, 
                         soilmoisture_b_5cm_avg_ave,
                         sensornode, ymd, yday,year, snow,
                         snowmelt_doy_infilled)))

phen_data_sos_pred = format_phen(phen_data_sos_pred)
phen_data_pop_pred = format_phen(phen_data_pop_pred)
phen_data_eos_pred = format_phen(phen_data_eos_pred)


#https://stats.stackexchange.com/questions/245038/variable-selection-combining-stepwise-regresion-and-cross-validation
# https://medium.com/analytics-vidhya/regularization-and-cross-validation-how-to-choose-the-penalty-value-lambda-1217fa4351e5

#define full models for prediction
# note definitely some things are collinear - but retained bc they explain
# excess variation

#full model
sos_all = glm(formula = status ~ GDD + days_since_snow+
                moist_7d + yday + temp_7d,
              family = binomial(link = "logit"),
              data = phen_data_sos %>%
                ungroup %>%
                filter(days_since_snow>0))

# #lasso regression will shrink params to zero
# library (glmnet)
# data = phen_data_sos %>%
#   ungroup %>%
#   filter(days_since_snow>0)
# x_sce=data %>% select(GDD, days_since_snow, moist_7d, yday, temp_7d) %>%
#   as.matrix()
# x_sce_scale = scale(x_sce)
# #https://glmnet.stanford.edu/articles/glmnetFamily.html
# cv <- cv.glmnet(x_sce, data$status, alpha = 1)
# newfit <- glmnet(x_sce, data$status, family = binomial(), lambda = cv$lambda.min)
# coef(newfit)
# 
# cv_scale <- cv.glmnet(x_sce_scale, data$status, alpha = 1)
# newfit_scale <- glmnet(x_sce_scale, data$status, family = binomial(), lambda = cv_scale$lambda.min)
# coef(newfit_scale)
# 
# 
# set.seed(1)
# x <- matrix(rnorm(500), ncol = 5)
# y <- rowSums(x[, 1:2]) + rnorm(100)
# biny <- ifelse(y > 0, 1, 0)  # binary data
# newfit <- glmnet(x, y, family = gaussian())
# newfit <- glmnet(x, biny, family = binomial())

#air
sos_air_all = glm(formula = status ~ GDD_air + yday + days_since_snow +temp_7d_air,
                  family = binomial(link = "logit"),
                  data = phen_data_sos %>%
                    ungroup %>%
                    filter(days_since_snow>0))

pop_all = glm(formula = status ~ GDD + days_since_snow + moist_7d + yday + temp_7d,
              family = binomial(link = "logit"),
              data = phen_data_pop %>%
                ungroup %>%
                filter(days_since_snow>0))

pop_air_all = glm(formula = status ~ GDD_air+ days_since_snow + yday + temp_7d_air,
                  family = binomial(link = "logit"),
                  data = phen_data_pop %>%
                    ungroup %>%
                    filter(days_since_snow>0))

eos_all = glm(formula = status ~ GDD + days_since_snow + moist_7d + yday + temp_7d,
              family = binomial(link = "logit"),
              data = phen_data_eos %>%
                ungroup %>%
                filter(days_since_snow>0))

eos_air_all = glm(formula = status ~ GDD_air + yday + days_since_snow + temp_7d_air,
                  family = binomial(link = "logit"),
                  data = phen_data_eos %>%
                    ungroup %>%
                    filter(days_since_snow>0))


#backwards stepwise
model.aic.backward_sos <- step(sos_all, direction = "backward", trace = 1)
summary (model.aic.backward_sos)

model.aic.backward_sos_air <- step(sos_air_all, direction = "backward", trace = 1)
summary (model.aic.backward_sos_air)

model.aic.backward_pop <- step(pop_all, direction = "backward", trace = 1)
summary (model.aic.backward_pop) # warmer springs make it earlier, but cooler summers make it later?

# model.aic.both_pop <- step(pop_all, direction = "both", trace = 1)
# summary (model.aic.both_pop)

# library (MuMIn)
# pop_all = glm(formula = status ~ GDD + days_since_snow + moist_7d + yday + temp_7d,
#               family = binomial(link = "logit"),
#               data = phen_data_pop %>%
#                 ungroup %>%
#                 filter(days_since_snow>0))
# dd <- dredge(pop_all)


model.aic.backward_pop_air <- step(pop_air_all, direction = "backward", trace = 1)
summary (model.aic.backward_pop_air) # warmer springs make it earlier, but cooler summers make it later?


model.aic.backward_eos <- step(eos_all, direction = "backward", trace = 1)
summary (model.aic.backward_eos)

model.aic.both_eos <- step(eos_all, direction = "both", trace = 1)
summary (model.aic.both_eos)

model.aic.backward_eos_air <- step(eos_air_all, direction = "backward", trace = 1)
summary (model.aic.backward_eos_air)


#predictive accuracy as a 2nd test against overfitting
#for that we need to calc survival
#product(1 - predicted hazards) from predicted on holdout set
#define function to cv spatial vs temporal w the surv models

#https://data.princeton.edu/wws509/notes/c7.pdf

#define function to cv and make oos predictions from the surv models
model_cv_surv = function(
  phen_data_fit, phen_data_pred, model){
  out = list()
  for (i in unique (phen_data_fit$sensornode)){
    phen_data_fit_i = phen_data_fit %>%
      filter(sensornode != i)
    phen_data_pred_i = phen_data_pred %>%
      filter(sensornode == i)
    mod = update (model, data = phen_data_fit_i %>%
                    filter(days_since_snow>0))
    surv  = 
      phen_data_pred_i %>%
      ungroup %>%
      filter(days_since_snow>0) %>%
      mutate(hazard = c(predict(mod, newdata =.,
                                type = 'response'))) %>%
      mutate(haz_compl = 1-hazard) %>%
      group_by(subject) %>%
      arrange(yday) %>%
      mutate(surv = cumprod(haz_compl)) %>%
      group_by(subject) %>%
      arrange(yday) %>%
      filter(surv<0.5) %>%
      slice(1) %>%
      select(subject, yday, year) %>%
      rename(doy_event_pred = yday) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, doy_event)%>%
                  distinct(),
                by = c("subject", "year"))
    out[[i]] = surv
  }
  out = out %>%bind_rows()
  spatial = Metrics::rmse(out$doy_event, out$doy_event_pred)
  
  #repeat drop by year
  out = list()
  for (i in unique (phen_data_fit$year)){
    phen_data_fit_i = phen_data_fit %>%
      filter(year != i)
    phen_data_pred_i = phen_data_pred %>%
      filter(year == i)
    mod = update (model, data = phen_data_fit_i %>%
                    filter(days_since_snow>0))
    surv  = 
      phen_data_pred_i %>%
      ungroup %>%
      filter(days_since_snow>0) %>%
      mutate(hazard = c(predict(mod, newdata =.,
                                type = 'response'))) %>%
      mutate(haz_compl = 1-hazard) %>%
      group_by(subject) %>%
      arrange(yday) %>%
      mutate(surv = cumprod(haz_compl)) %>%
      group_by(subject) %>%
      arrange(yday) %>%
      filter(surv<0.5) %>%
      slice(1) %>%
      select(subject, yday, year) %>%
      rename(doy_event_pred = yday) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, doy_event)%>%
                  distinct(),
                by = c("subject", "year"))
    out[[i]] = surv
  }
  out = out %>%bind_rows()
  temporal = Metrics::rmse(out$doy_event, out$doy_event_pred)
  return(list (spatial = spatial, temporal = temporal))
}

#define function to cv and make oos predictions from threshold models
model_cv_threshold = function(
  phen_data_fit, phen_data_pred, thresh){
  out = list()
  #drop spatial
  for (i in unique (phen_data_fit$sensornode)){
    phen_data_fit_i = phen_data_fit %>%
      filter(sensornode != i) %>%
      ungroup %>%
      filter(yday == doy_event)
    phen_data_pred_i = phen_data_pred %>%
      filter(sensornode == i)
    mod = lm (formula =as.formula(paste(thresh, paste(1), sep = '~')),
              phen_data_fit_i)
    
    surv = phen_data_pred_i %>%
      select(subject, year) %>%
      distinct() %>%
      mutate(mod_pred = c(predict(mod, newdata =.,
                                  type = 'response'))) %>%
      select(subject, mod_pred, year) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, !!thresh, doy_event, yday)%>%
                  distinct(),
                by = c("subject", "year")) %>%
      #should be able to do this with !! but it seems unhappy
      rename (varcol = !!thresh) %>%
      filter(varcol> mod_pred) %>%
      group_by(subject) %>%
      arrange(yday) %>%
      slice(1) %>%
      select(subject, yday, doy_event) %>%
      rename(doy_event_pred = yday)
    out[[i]] = surv
  }
  out = out %>%bind_rows()
  spatial = Metrics::rmse(out$doy_event, out$doy_event_pred)
  
  #repeat drop by year
  out = list()
  for (i in unique (phen_data_fit$year)){
    phen_data_fit_i = phen_data_fit %>%
      filter(year != i) %>%
      ungroup %>%
      filter(yday == doy_event)
    phen_data_pred_i = phen_data_pred %>%
      filter(year == i) 
    
    mod = lm (formula =as.formula(paste(thresh, paste(1), sep = '~')),
              phen_data_fit_i)
    surv = phen_data_pred_i %>%
      select(subject, year) %>%
      distinct() %>%
      mutate(mod_pred = c(predict(mod, newdata =.,
                                  type = 'response'))) %>%
      select(subject, mod_pred, year) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, !!thresh, doy_event, yday)%>%
                  distinct(),
                by = c("subject", "year")) %>%
      #should be able to do this with !! but it seems unhappy
      rename (varcol = !!thresh) %>%
      filter(varcol> mod_pred) %>%
      group_by(subject) %>%
      arrange(yday) %>%
      slice(1) %>%
      select(subject, yday, doy_event) %>%
      rename(doy_event_pred = yday)
    
    out[[i]] = surv
  }
  out = out %>%bind_rows()
  temporal = Metrics::rmse(out$doy_event, out$doy_event_pred)
  return(list (spatial = spatial, temporal = temporal))
}

#compare model fits EOS####

#drop each term in model
#full model
eos_fitted = glm(formula = formula(model.aic.backward_eos),
                 family = binomial(link = "logit"),
                 data = phen_data_eos %>%
                   ungroup %>%
                   filter(days_since_snow>0))

#both spatial and temporal predictions improved by adding yday to the model
eos_full_rmse = model_cv_surv(phen_data_eos, phen_data_eos_pred,
                              eos_fitted)


eos_fitted_air = glm(formula = formula(model.aic.backward_eos_air),
                     family = binomial(link = "logit"),
                     data = phen_data_eos %>%
                       ungroup %>%
                       filter(days_since_snow>0))
#both spatial and temporal predictions improved by adding yday to the model
eos_full_rmse_air = model_cv_surv(phen_data_eos, phen_data_eos_pred,
                                  eos_fitted_air)


#surv wo days_since_snow
eos_fitted_drop_days_since_snow = glm(formula = update(formula(model.aic.backward_eos), ~. -days_since_snow),
                                      family = binomial(link = "logit"),
                                      data = phen_data_eos %>%
                                        ungroup %>%
                                        filter(days_since_snow>0))

eos_drop_days_since_snow_rmse = model_cv_surv(phen_data_eos, phen_data_eos_pred,
                                              eos_fitted_drop_days_since_snow)
#surv wo moist
eos_fitted_drop_moist_7d = glm(formula = update(formula(model.aic.backward_eos), ~. -moist_7d),
                               family = binomial(link = "logit"),
                               data = phen_data_eos %>%
                                 ungroup %>%
                                 filter(days_since_snow>0))

eos_drop_moist_7d_rmse = model_cv_surv(phen_data_eos, phen_data_eos_pred,
                                       eos_fitted_drop_moist_7d)
#surv wo yday
eos_fitted_drop_yday= glm(formula = update(formula(model.aic.backward_eos), ~. -yday),
                          family = binomial(link = "logit"),
                          data = phen_data_eos %>%
                            ungroup %>%
                            filter(days_since_snow>0))

eos_drop_yday_rmse = model_cv_surv(phen_data_eos, phen_data_eos_pred,
                                   eos_fitted_drop_yday)
#surv wo temp
eos_fitted_drop_temp_7d= glm(formula = update(formula(model.aic.backward_eos), ~. -temp_7d),
                             family = binomial(link = "logit"),
                             data = phen_data_eos %>%
                               ungroup %>%
                               filter(days_since_snow>0))

eos_drop_temp_7d_rmse = model_cv_surv(phen_data_eos, phen_data_eos_pred,
                                      eos_fitted_drop_temp_7d)

#repeat for air
#surv wo yday
eos_fitted_drop_yday_air= glm(formula = update(formula(model.aic.backward_eos_air), ~. -yday),
                              family = binomial(link = "logit"),
                              data = phen_data_eos %>%
                                ungroup %>%
                                filter(days_since_snow>0))

eos_drop_yday_rmse_air = model_cv_surv(phen_data_eos, phen_data_eos_pred,
                                       eos_fitted_drop_yday_air)

#surv wo GDD
eos_fitted_drop_GDD_air= glm(formula = update(formula(model.aic.backward_eos_air), ~. -GDD),
                             family = binomial(link = "logit"),
                             data = phen_data_eos %>%
                               ungroup %>%
                               filter(days_since_snow>0))

eos_drop_GDD_rmse_air = model_cv_surv(phen_data_eos, phen_data_eos_pred,
                                      eos_fitted_drop_GDD_air)



#3 threshold/null models
eos_rmse_days_since_snow = model_cv_threshold(phen_data_eos, phen_data_eos_pred, 'days_since_snow')
eos_rmse_GDD =model_cv_threshold(phen_data_eos, phen_data_eos_pred, 'GDD')
eos_rmse_GDD_air =model_cv_threshold(phen_data_eos, phen_data_eos_pred, 'GDD_air')
eos_rmse_null =model_cv_threshold(phen_data_eos%>%
                                    mutate(yday1 = yday), phen_data_eos_pred%>%
                                    mutate(yday1 = yday), 'yday1')



#plots
all_rmse_eos = list()
for (mod in c(
  'eos_full_rmse', 
  'eos_full_rmse_air',
  'eos_drop_yday_rmse', 
  'eos_drop_yday_rmse_air', 
  'eos_drop_GDD_rmse_air', 
  'eos_drop_temp_7d_rmse', 
  'eos_drop_moist_7d_rmse', 
  'eos_drop_days_since_snow_rmse', 
  'eos_rmse_null', 
  'eos_rmse_GDD',
  'eos_rmse_GDD_air',
  'eos_rmse_days_since_snow')){
  all_rmse_eos[[mod]] = as.data.frame(get(mod))
}

all_rmse_eos = all_rmse_eos%>%
  data.table::rbindlist(., idcol = 'mod', use.names = TRUE)%>%
  pivot_longer(cols = c(spatial, temporal),
               names_to = 'type',
               values_to = 'rmse') %>%
  mutate(test = ifelse(grepl('drop', mod), 'drop', 'threshold'))%>%
  #relabel for nice plots
  mutate(
    mod = case_when(
      mod == 'eos_full_rmse' ~ 'time-to-event \n (soil temperature)',
      mod == 'eos_full_rmse_air' ~ 'time-to-event \n (air temperature)',
      mod == 'eos_rmse_days_since_snow' ~ 'Threshold \n (days past snowmelt)',
      mod == 'eos_rmse_null' ~ 'Threshold \n (day of year)',
      mod == 'eos_rmse_GDD' ~ 'Threshold \n (GDD soil)',
      mod == 'eos_rmse_GDD_air' ~ 'Threshold \n (GDD air)',
      mod == 'eos_drop_yday_rmse' ~ '- day of year',
      mod == 'eos_drop_temp_7d_rmse' ~ '- 7d running mean soil temperature',
      mod == 'eos_drop_moist_7d_rmse' ~ '- 7d running mean soil moisture',
      mod == 'eos_drop_days_since_snow_rmse' ~ '- days since snowmelt',
      mod == 'eos_drop_yday_rmse_air' ~ '- day of year (air temperature model)',
      mod == 'eos_drop_GDD_rmse_air' ~ '- GDD (air temperature model)',
      TRUE ~ mod
    )
  )

eos_vs_thresh <- ggplot(all_rmse_eos %>%
                          mutate(
                            mod = fct_reorder(mod, rmse, .fun = mean, .desc = FALSE)
                          )%>%
                          filter(test == 'threshold'), aes(fill=type, y=rmse, x=mod)) + 
  geom_bar(position="dodge", stat="identity")+
  #ggtitle("Comparison of time-to-event against \n threshold models for end of season")+
  ggtitle("senescence")+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1))+
  scale_fill_grey()+ guides(fill=guide_legend(title="test set"))+
  ylab('RMSE')

#define function to plot effects
plot_conditional_effect = function(dataset,
                                   phenometric = 'greenup', model,
                                   all_rmse, soil_or_air = 'soil'
){
  #dataset = phen_data_sos
  #min_event = min_sos
  #max_event = max_sos
  #model = sos_fitted
  #term = 'GDD'
  #all_rmse = all_rmse_sos
  
  #convert to deltas
  full = all_rmse %>%
    filter(grepl('full', model))
  
  if (soil_or_air == 'air'){
    full = full %>%
      filter(grepl('air', model))
  }else{
    full = full %>%
      filter(!grepl('air', model))
  }
  all_rmse_delta = all_rmse %>%
    mutate(full_spatial = full$spatial,
           full_temporal = full$temporal)%>%
    mutate(spatial = spatial - full_spatial,
           temporal = temporal - full_temporal) %>%
    select(model, spatial, temporal) %>%
    pivot_longer(cols = c(spatial, temporal),
                 names_to = 'type',
                 values_to = 'delta_rmse')
  
  if (soil_or_air == 'air'){
    all_rmse_delta = all_rmse_delta %>%
      filter(grepl('air', model))
  }else{
    all_rmse_delta = all_rmse_delta %>%
      filter(!grepl('air', model))
  }
  
  min_event = min(dataset$doy_event)
  max_event = max(dataset$doy_event)
  avg_values = dataset%>%
    filter(yday>=min_event&yday<max_event)
  plot_list = list()
  
  for (term in all.vars(formula(model)[-2])){
    nice_titles = if(term == 'yday'){
      'day of year'
    }else if (term == 'temp_7d'){
      'weekly soil temperature'
    }else if (term == 'days_since_snow'){
      'days past snowmelt'
    }else if (term == 'moist_7d'){
      'weekly soil moisture'
    }else if (term == 'GDD'){
      'GDD (soil)'
    }else if (term == 'GDD_air'){
      'GDD (air)'
    }else if (term == 'temp_7d_air'){
      'weekly air temperature'
    }
    vars = all.vars(formula(model)[-2])
    vars = vars[!vars==term]
    out = list()
    for (v in vars){
      out[[v]] = mean (avg_values[[v]])
    }
    preds = data.frame(dataset%>%
                         filter(yday>=min_event&yday<max_event)%>%
                         select(!!term))%>%
      distinct()
    for (v in vars){
      preds[[v]] = out[[v]]
    }
    preds$hazard = predict(model, preds, type = 'response')
    
    p1 <- ggplot (preds, aes_string(x=term, y='hazard'))+
      geom_line()+
      ylab(NULL)+
      theme_classic()+ylim(0,1)+
      xlab(nice_titles)
    
    inset.plot<-ggplot(all_rmse_delta  %>%
                         filter(grepl('drop', model)&grepl(term, model)), aes(fill=type, y=delta_rmse, x=model))+
      
      geom_bar(position="dodge", stat="identity", color = 'black')+
      theme(axis.text.x = element_text(angle = 90)) +
      theme_classic()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      scale_fill_grey()+ 
      ylim(-1,16)+
      theme(legend.position = 'none')+
      ylab(expression(Delta*"RMSE"))
    
    plot.with.inset <-
      ggdraw() +
      draw_plot(p1) +
      draw_plot(inset.plot, x = 0.5, y = .7, width = .5, height = .3)
    
    plot_list[[term]] <- plot.with.inset
  }
  return (plot_list)
}


all_rmse_eos = list()
for (mod in c(
  'eos_full_rmse', # best is 9,10,
  'eos_full_rmse_air',
  'eos_drop_yday_rmse', # drop yday it's 14, 25
  'eos_drop_yday_rmse_air', # drop yday it's 14, 25
  'eos_drop_GDD_rmse_air', # drop yday it's 14, 25
  'eos_drop_temp_7d_rmse', # drop temp its 12, 13
  'eos_drop_moist_7d_rmse', # drop moist it's 13, 16
  'eos_drop_days_since_snow_rmse', #days since snow is still 9,10
  'eos_rmse_null', # null is 13, 15
  'eos_rmse_GDD',
  'eos_rmse_GDD_air',
  'eos_rmse_days_since_snow')){
  all_rmse_eos[[mod]] = as.data.frame(get(mod))
}

all_rmse_eos = all_rmse_eos%>%
  data.table::rbindlist(, idcol = 'model')

#all_rmse_sos = all_rmse_sos %>% data.table::rbindlist(, idcol = 'model')
eos_effects <- plot_conditional_effect(dataset =phen_data_eos,
                                       model =eos_fitted,
                                       all_rmse = all_rmse_eos,
                                       soil_or_air = 'soil')

#labels
x.grob <- textGrob("senescence", 
                   gp=gpar(col="black", fontsize=20))
y.grob <- textGrob("h(x)", 
                   gp=gpar(col="black", fontsize=15), rot=90)

eos_all<-grid.arrange(arrangeGrob(cowplot::plot_grid(plotlist =eos_effects), left = y.grob, top = x.grob))
jpeg('plots/eos_all.jpg', width =600, height =600)
ggdraw(eos_all)
dev.off()

# repeat for peak#####
#drop each in turn
#
pop_fitted = glm(formula = formula(model.aic.backward_pop),
                 family = binomial(link = "logit"),
                 data = phen_data_pop %>%
                   ungroup %>%
                   filter(days_since_snow>0))

#both spatial and temporal predictions improved by adding yday to the model
pop_full_rmse = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                              pop_fitted)

#
pop_fitted_air = glm(formula = formula(model.aic.backward_pop_air),
                     family = binomial(link = "logit"),
                     data = phen_data_pop %>%
                       ungroup %>%
                       filter(days_since_snow>0))

#both spatial and temporal predictions improved by adding yday to the model
pop_full_rmse_air = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                  pop_fitted_air)



pop_fitted_drop_GDD = glm(formula = update(formula(model.aic.backward_pop), ~. -GDD),
                          family = binomial(link = "logit"),
                          data = phen_data_pop %>%
                            ungroup %>%
                            filter(days_since_snow>0))

pop_drop_GDD_rmse = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                  pop_fitted_drop_GDD)


pop_fitted_drop_days_since_snow = glm(formula = update(formula(model.aic.backward_pop), ~. -days_since_snow),
                                      family = binomial(link = "logit"),
                                      data = phen_data_pop %>%
                                        ungroup %>%
                                        filter(days_since_snow>0))

pop_drop_days_since_snow_rmse = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                              pop_fitted_drop_days_since_snow)




pop_fitted_drop_yday= glm(formula = update(formula(model.aic.backward_pop), ~. -yday),
                          family = binomial(link = "logit"),
                          data = phen_data_pop %>%
                            ungroup %>%
                            filter(days_since_snow>0))

pop_drop_yday_rmse = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                   pop_fitted_drop_yday)

pop_fitted_drop_temp_7d= glm(formula = update(formula(model.aic.backward_pop), ~. -temp_7d),
                             family = binomial(link = "logit"),
                             data = phen_data_pop %>%
                               ungroup %>%
                               filter(days_since_snow>0))

pop_drop_temp_7d_rmse = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                      pop_fitted_drop_temp_7d)


pop_fitted_drop_moist_7d= glm(formula = update(formula(model.aic.backward_pop), ~. -moist_7d),
                              family = binomial(link = "logit"),
                              data = phen_data_pop %>%
                                ungroup %>%
                                filter(days_since_snow>0))
pop_drop_moist_7d_rmse = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                       pop_fitted_drop_moist_7d)
#air drop 1 term
pop_fitted_drop_temp_7d_air= glm(formula = update(formula(model.aic.backward_pop_air), ~. -temp_7d_air),
                                 family = binomial(link = "logit"),
                                 data = phen_data_pop %>%
                                   ungroup %>%
                                   filter(days_since_snow>0))
pop_drop_temp_7d_rmse_air = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                          pop_fitted_drop_temp_7d_air)

pop_fitted_drop_GDD_air= glm(formula = update(formula(model.aic.backward_pop_air), ~. -GDD_air),
                             family = binomial(link = "logit"),
                             data = phen_data_pop %>%
                               ungroup %>%
                               filter(days_since_snow>0))
pop_drop_GDD_air_rmse_air = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                          pop_fitted_drop_GDD_air)

pop_fitted_drop_days_since_snow_air= glm(formula = update(formula(model.aic.backward_pop_air), ~.
                                                          -days_since_snow),
                                         family = binomial(link = "logit"),
                                         data = phen_data_pop %>%
                                           ungroup %>%
                                           filter(days_since_snow>0))
pop_drop_days_since_snow_rmse_air = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                                  pop_fitted_drop_days_since_snow_air)

pop_fitted_drop_yday_air= glm(formula = update(formula(model.aic.backward_pop_air), ~.
                                               -yday),
                              family = binomial(link = "logit"),
                              data = phen_data_pop %>%
                                ungroup %>%
                                filter(days_since_snow>0))
pop_drop_yday_rmse_air = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                       pop_fitted_drop_yday_air)

pop_fitted_drop_temp_7d_air= glm(formula = update(formula(model.aic.backward_pop_air), ~.
                                                  -temp_7d_air),
                                 family = binomial(link = "logit"),
                                 data = phen_data_pop %>%
                                   ungroup %>%
                                   filter(days_since_snow>0))
pop_drop_temp_7d_rmse_air = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                          pop_fitted_drop_temp_7d_air)



#null/threshold models
pop_rmse_days_since_snow =model_cv_threshold(phen_data_pop, phen_data_pop_pred, 'days_since_snow')
pop_rmse_GDD =model_cv_threshold(phen_data_pop, phen_data_pop_pred, 'GDD')
pop_rmse_GDD_air =model_cv_threshold(phen_data_pop, phen_data_pop_pred, 'GDD_air')
pop_rmse_null =model_cv_threshold(phen_data_pop%>%
                                    mutate(yday1 = yday), phen_data_pop_pred%>%
                                    mutate(yday1 = yday), 'yday1')

#plots
all_rmse_pop = list()
for (mod in c(
  'pop_full_rmse', 
  'pop_full_rmse_air', 
  'pop_drop_yday_rmse', 
  'pop_drop_yday_rmse_air', 
  'pop_drop_temp_7d_rmse', 
  'pop_drop_temp_7d_rmse_air', 
  'pop_drop_days_since_snow_rmse', 
  'pop_drop_GDD_rmse', 
  'pop_drop_GDD_air_rmse_air', 
  'pop_drop_moist_7d_rmse', 
  'pop_rmse_null', 
  'pop_rmse_GDD',
  'pop_rmse_GDD_air',
  'pop_rmse_days_since_snow')){
  all_rmse_pop[[mod]] = as.data.frame(get(mod))
}

all_rmse_pop = all_rmse_pop %>%
  data.table::rbindlist(., idcol = 'mod', use.names = TRUE)%>%
  pivot_longer(cols = c(spatial, temporal),
               names_to = 'type',
               values_to = 'rmse') %>%
  mutate(test = ifelse(grepl('drop', mod), 'drop', 'threshold'))%>%
  #relabel for nice plots
  mutate(
    mod = case_when(
      mod == 'pop_full_rmse' ~ 'time-to-event \n (soil temperature)',
      mod == 'pop_full_rmse_air' ~ 'time-to-event \n (air temperature)',
      mod == 'pop_rmse_days_since_snow' ~ 'Threshold \n (days past snowmelt)',
      mod == 'pop_rmse_null' ~ 'Threshold \n (day of year)',
      mod == 'pop_rmse_GDD' ~ 'Threshold \n (GDD soil)',
      mod == 'pop_rmse_GDD_air' ~ 'Threshold \n (GDD air)',
      mod == 'pop_drop_yday_rmse' ~ '- day of year',
      mod == 'pop_drop_GDD_rmse' ~ '- GDD (soil)',
      mod == 'pop_drop_GDD_air_rmse_air' ~ '- GDD (air)',
      mod == 'pop_drop_temp_7d_rmse' ~ '- 7d running mean soil temperature',
      mod == 'pop_drop_moist_7d_rmse' ~ '- 7d running mean soil moisture',
      mod == 'pop_drop_days_since_snow_rmse' ~ '- days since snowmelt',
      mod == 'pop_drop_yday_rmse_air' ~ '- day of year (air temperature model)',
      mod == 'pop_drop_temp_7d_rmse_air' ~ '- 7d running mean air temperature',
      mod == 'pop_drop_days_since_snow_rmse_air' ~ '- days since snowmelt (air temperature model)',
      TRUE ~ mod
    )
  )

pop_vs_thresh <- ggplot(all_rmse_pop %>%
                          mutate(
                            mod = fct_reorder(mod, rmse, .fun = mean, .desc = FALSE)
                          )%>%
                          filter(test == 'threshold'), aes(fill=type, y=rmse, x=mod)) + 
  geom_bar(position="dodge", stat="identity")+
  ggtitle("peak")+
  theme_classic()+
  ylim(0, 25)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1))+
  scale_fill_grey()+ guides(fill=guide_legend(title="test set"))+
  ylab('RMSE')


#plots
all_rmse_pop = list()
for (mod in c(
  'pop_full_rmse', 
  'pop_full_rmse_air', 
  'pop_drop_yday_rmse', 
  'pop_drop_yday_rmse_air', 
  'pop_drop_temp_7d_rmse', 
  'pop_drop_temp_7d_rmse_air', 
  'pop_drop_days_since_snow_rmse', 
  'pop_drop_GDD_rmse', 
  'pop_drop_GDD_air_rmse_air', 
  'pop_drop_moist_7d_rmse', 
  'pop_rmse_null', 
  'pop_rmse_GDD',
  'pop_rmse_GDD_air',
  'pop_rmse_days_since_snow')){
  all_rmse_pop[[mod]] = as.data.frame(get(mod))
}

all_rmse_pop = all_rmse_pop %>%
  data.table::rbindlist(, idcol = 'model')

#all_rmse_sos = all_rmse_sos %>% data.table::rbindlist(, idcol = 'model')
pop_effects <- plot_conditional_effect(dataset =phen_data_pop,
                                       model =pop_fitted,
                                       all_rmse = all_rmse_pop,
                                       soil_or_air = 'soil')

#labels
x.grob <- textGrob("peak", 
                   gp=gpar(col="black", fontsize=20))
y.grob <- textGrob("h(x)", 
                   gp=gpar(col="black", fontsize=15), rot=90)

pop_all<-grid.arrange(arrangeGrob(cowplot::plot_grid(plotlist =pop_effects), left = y.grob, top = x.grob))
jpeg('plots/pop_all.jpg', width =600, height =600)
ggdraw(pop_all)
dev.off()


###end repeat for peak

## run for sos

sos_fitted = glm(formula = formula(model.aic.backward_sos),
                 family = binomial(link = "logit"),
                 data = phen_data_sos %>%
                   ungroup %>%
                   filter(days_since_snow>0))

ggplot (broom::tidy(sos_fitted, conf.int = TRUE) %>%
          #rename(phenophase = name) %>%
          filter(term !='(Intercept)'),
        aes(xmin=conf.low, xmax = conf.high, y=term))+#, color = spec))+
  geom_vline(aes(xintercept =0))+
  geom_point(aes(x =estimate))+
  geom_linerange(position=position_dodge(width=0.5))


#both spatial and temporal predictions improved by adding yday to the model
sos_full_rmse = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                              sos_fitted)

sos_fitted_air = glm(formula = formula(model.aic.backward_sos_air),
                     family = binomial(link = "logit"),
                     data = phen_data_sos %>%
                       ungroup %>%
                       filter(days_since_snow>0))
ggplot (broom::tidy(sos_fitted_air, conf.int = TRUE) %>%
          #rename(phenophase = name) %>%
          filter(term !='(Intercept)'),
        aes(xmin=conf.low, xmax = conf.high, y=term))+#, color = spec))+
  geom_vline(aes(xintercept =0))+
  geom_point(aes(x =estimate))+
  geom_linerange(position=position_dodge(width=0.5))



#both spatial and temporal predictions improved by adding yday to the model
sos_full_rmse_air = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                  sos_fitted_air)


sos_fitted_drop_GDD = glm(formula = update(formula(model.aic.backward_sos), ~. -GDD),
                          family = binomial(link = "logit"),
                          data = phen_data_sos %>%
                            ungroup %>%
                            filter(days_since_snow>0))

sos_drop_GDD_rmse = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                  sos_fitted_drop_GDD)
sos_fitted_drop_days_since_snow = glm(formula = update(formula(model.aic.backward_sos), ~. -days_since_snow),
                                      family = binomial(link = "logit"),
                                      data = phen_data_sos %>%
                                        ungroup %>%
                                        filter(days_since_snow>0))

sos_drop_days_since_snow_rmse = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                              sos_fitted_drop_days_since_snow)

sos_fitted_drop_days_since_snow_air = glm(formula = update(formula(model.aic.backward_sos_air ), ~. -days_since_snow),
                                          family = binomial(link = "logit"),
                                          data = phen_data_sos %>%
                                            ungroup %>%
                                            filter(days_since_snow>0))

sos_drop_days_since_snow_rmse_air = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                                  sos_fitted_drop_days_since_snow_air)



sos_fitted_drop_temp_7d = glm(formula = update(formula(model.aic.backward_sos), ~. -temp_7d),
                              family = binomial(link = "logit"),
                              data = phen_data_sos %>%
                                ungroup %>%
                                filter(days_since_snow>0))

sos_drop_temp_7d_rmse = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                      sos_fitted_drop_temp_7d)

sos_fitted_drop_temp_7d_air = glm(formula = update(formula(model.aic.backward_sos_air), ~. -temp_7d_air),
                                  family = binomial(link = "logit"),
                                  data = phen_data_sos %>%
                                    ungroup %>%
                                    filter(days_since_snow>0))

sos_drop_temp_7d_rmse_air = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                          sos_fitted_drop_temp_7d_air)


sos_fitted_drop_yday= glm(formula = update(formula(model.aic.backward_sos), ~. -yday),
                          family = binomial(link = "logit"),
                          data = phen_data_sos %>%
                            ungroup %>%
                            filter(days_since_snow>0))

sos_drop_yday_rmse = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                   sos_fitted_drop_yday)

sos_fitted_drop_yday_air= glm(formula = update(formula(model.aic.backward_sos_air), ~. -yday),
                              family = binomial(link = "logit"),
                              data = phen_data_sos %>%
                                ungroup %>%
                                filter(days_since_snow>0))

sos_drop_yday_rmse_air = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                       sos_fitted_drop_yday_air)




#both improve a lot over the null model of yday
#null/threshold models
sos_rmse_days_since_snow =model_cv_threshold(phen_data_sos, phen_data_sos_pred, 'days_since_snow')
sos_rmse_GDD =model_cv_threshold(phen_data_sos, phen_data_sos_pred, 'GDD')
sos_rmse_GDD_air =model_cv_threshold(phen_data_sos, phen_data_sos_pred, 'GDD_air')
sos_rmse_null =model_cv_threshold(phen_data_sos%>%
                                    mutate(yday1 = yday), phen_data_sos_pred%>%
                                    mutate(yday1 = yday), 'yday1')

#plots
all_rmse_sos = list()
for (mod in c(
  'sos_full_rmse', 
  'sos_drop_yday_rmse', 
  'sos_drop_GDD_rmse', 
  'sos_drop_days_since_snow_rmse', 
  'sos_drop_temp_7d_rmse', 
  'sos_rmse_null', 
  'sos_rmse_GDD',
  'sos_rmse_GDD_air',
  'sos_rmse_days_since_snow',
  'sos_full_rmse_air', 
  'sos_drop_yday_rmse_air', 
  'sos_drop_temp_7d_rmse_air', 
  'sos_drop_days_since_snow_rmse_air')){ 
  all_rmse_sos[[mod]] = as.data.frame(get(mod))
}


all_rmse_sos = all_rmse_sos%>%
  data.table::rbindlist(., idcol = 'mod', use.names = TRUE) %>%
  pivot_longer(cols = c(spatial, temporal),
               names_to = 'type',
               values_to = 'rmse') %>%
  mutate(test = ifelse(grepl('drop', mod), 'drop', 'threshold'))%>%
  #relabel for nice plots
  mutate(
    mod = case_when(
      mod == 'sos_full_rmse' ~ 'time-to-event \n (soil temperature)',
      mod == 'sos_full_rmse_air' ~ 'time-to-event \n (air temperature)',
      mod == 'sos_rmse_days_since_snow' ~ 'Threshold \n (days past snowmelt)',
      mod == 'sos_rmse_null' ~ 'Threshold\n (day of year)',
      mod == 'sos_rmse_GDD' ~ 'Threshold \n (GDD soil)',
      mod == 'sos_rmse_GDD_air' ~ 'Threshold \n (GDD air)',
      mod == 'sos_drop_yday_rmse' ~ '- day of year',
      mod == 'sos_drop_GDD_rmse' ~ '- GDD \n (soil temperature model)',
      mod == 'sos_drop_temp_7d_rmse' ~ '- 7d running mean soil temperature',
      mod == 'sos_drop_temp_7d_rmse_air' ~ '- 7d running mean air temperature',
      mod == 'sos_drop_moist_7d_rmse' ~ '- 7d running mean soil moisture',
      mod == 'sos_drop_days_since_snow_rmse' ~ '- days since snowmelt',
      mod == 'sos_drop_yday_rmse_air' ~ '- day of year (air temperature model)',
      mod == 'sos_drop_days_since_snow_rmse_air' ~ '- days since snowmelt (air temperature model)',
      mod == 'sos_drop_GDD_rmse_air' ~ '- GDD (air temperature model)',
      TRUE ~ mod
    )
  )

sos_vs_thresh <- ggplot(all_rmse_sos %>%
                          mutate(
                            mod = fct_reorder(mod, rmse, .fun = mean, .desc = FALSE)
                          )%>%
                          filter(test == 'threshold'), aes(fill=type, y=rmse, x=mod)) + 
  geom_bar(position="dodge", stat="identity")+
  ggtitle("greenup")+
  theme_classic()+
  ylim(0, 25)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1))+
  scale_fill_grey()+ guides(fill=guide_legend(title="test set"))+
  ylab('RMSE')

# conditional effects plots

min_sos = min(phen_data_sos$doy_event)
max_sos = max(phen_data_sos$doy_event)

#regrab the rmses to make deltas
all_rmse_sos = list()
for (mod in c(
  'sos_full_rmse', # best is 9,10
  'sos_drop_yday_rmse', # drop yday it's 14, 25
  'sos_drop_GDD_rmse', # drop moist it's 13, 16
  'sos_drop_days_since_snow_rmse', # drop moist it's 13, 16
  'sos_drop_temp_7d_rmse', # drop moist it's 13, 16
  'sos_rmse_null', # null is 13, 15
  'sos_rmse_GDD',
  'sos_rmse_GDD_air',
  'sos_rmse_days_since_snow',
  'sos_full_rmse_air', # best is 9,10
  'sos_drop_yday_rmse_air', # drop yday it's 14, 25
  'sos_drop_temp_7d_rmse_air', # drop moist it's 13, 16
  'sos_drop_days_since_snow_rmse_air')){ # drop moist it's 13, 16
  all_rmse_sos[[mod]] = as.data.frame(get(mod))
}

all_rmse_sos = all_rmse_sos %>% data.table::rbindlist(, idcol = 'model')

#all_rmse_sos = all_rmse_sos %>% data.table::rbindlist(, idcol = 'model')
sos_effects <- plot_conditional_effect(dataset =phen_data_sos,
                                       model =sos_fitted,
                                       all_rmse = all_rmse_sos,
                                       soil_or_air = 'soil') 


x.grob <- textGrob("greenup", 
                   gp=gpar(col="black", fontsize=20))
y.grob <- textGrob("h(x)", 
                   gp=gpar(col="black", fontsize=15), rot=90)

sos_all<-grid.arrange(arrangeGrob(cowplot::plot_grid(plotlist =sos_effects), left = y.grob, top = x.grob))
jpeg('plots/sos_all.jpg', width =600, height =600)
ggdraw(sos_all)
dev.off()

#####

#all together
# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  eos_vs_thresh + theme(legend.box.margin = margin(0, 0, 0, 12))
)

#tst = sos_vs_thresh+scale_x_discrete(labels = scales::wrap_format(16))

#ggtitle("Comparison of time-to-event against \n threshold models for start of season")+
prow <- cowplot::plot_grid(sos_vs_thresh+theme(legend.position = 'none',axis.title.y = element_blank()),
                           pop_vs_thresh+theme(legend.position = 'none', axis.title.y = element_blank()),
                           eos_vs_thresh+theme(legend.position = 'none', axis.title.y = element_blank()), nrow = 1)

y.grob <- textGrob("RMSE", 
                   gp=gpar(col="black", fontsize=10), rot=90)
p1 = plot_grid(prow, legend, rel_widths = c(3,.4))

rmse_all<-grid.arrange(arrangeGrob(p1, left = y.grob))
jpeg('plots/rmse_all.jpg', width =800, height =200)
ggdraw(rmse_all)
dev.off()



