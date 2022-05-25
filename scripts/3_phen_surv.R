#code to fit discrete time logistic survival model to nwt phenocam data
# cross-validate; generate predictions
#sce 8/12/2021

#todo

# setup -------------------------------------------------------------------

library (tidyverse)
library(GGally)
library (cowplot)
library (grid)
library (gridExtra)
library(car)  # result from Anova() from package car for summary stats
library (patchwork)

# read data ---------------------------------------------------------------


# munge a bit to proper format####
# met data comes from script 1a_infill_met.R
met_data <- readr::read_csv("data_deriv/met_data_infilled.csv") %>%
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

#read the phenometrics, made from 2_analyze_saddle_phenocams
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

#just to get the DDs etc calculated for figs
phen_data_met_all = make_phen_data(sumry %>%
                                     mutate(met = 365) , "met") %>%
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

phen_data_met_all = phen_data_met_all %>%
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


# analyze surv ------------------------------------------------------------


#define functions for survival
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
    # mutate(GDD_scale = scale(GDD),
    #        GDD_air_scale = scale(GDD_air),
    #        moist_7d_scale = scale(moist_7d), 
    #        temp_7d_air_scale = scale(temp_7d_air), 
    #        temp_7d_scale = scale(temp_7d),
    #        yday_scale = scale(yday),
    #        days_since_snow_scale = scale(days_since_snow))%>%
    ungroup() 
}

phen_data_sos = format_phen(phen_data_sos)
phen_data_pop = format_phen(phen_data_pop)
phen_data_eos = format_phen(phen_data_eos)

phen_data_met_all = format_phen(phen_data_met_all)


# add met plots -----------------------------------------------------------

temp_plot <- ggplot (phen_data_met_all,
                     aes(x=yday, y=temp_7d))+
  geom_vline(data = sumry %>%filter(year>2017), aes(xintercept = sos), color = 'green', alpha = 0.2)+
  geom_vline(data = sumry %>%filter(year>2017), aes(xintercept = pop), color = 'blue', alpha = 0.2)+
  geom_vline(data = sumry %>%filter(year>2017), aes(xintercept = eos), color = 'brown', alpha = 0.2)+
  geom_line(aes(group = subject))+
  facet_wrap(~year, nrow =1)+ theme_classic()+
  xlab(NULL) + ylab ('soil temperature °C')



temp_plot_air <-ggplot (phen_data_met_all,
                        aes(x=yday, y=temp_7d_air))+
  geom_vline(data = sumry %>%filter(year>2017), aes(xintercept = sos), color = 'green', alpha = 0.2)+
  geom_vline(data = sumry %>%filter(year>2017), aes(xintercept = pop), color = 'blue', alpha = 0.2)+
  geom_vline(data = sumry %>%filter(year>2017), aes(xintercept = eos), color = 'brown', alpha = 0.2)+
  geom_line(aes(group = subject))+
  facet_wrap(~year, nrow =1)+ theme_classic() +
  xlab(NULL) + ylab ('air temperature °C')



moist_plot <-ggplot (phen_data_met_all,
                     aes(x=yday, y=moist_7d))+
  geom_vline(data = sumry %>%filter(year>2017), aes(xintercept = sos), color = 'green', alpha = 0.2)+
  geom_vline(data = sumry %>%filter(year>2017), aes(xintercept = pop), color = 'blue', alpha = 0.2)+
  geom_vline(data = sumry %>%filter(year>2017), aes(xintercept = eos), color = 'brown', alpha = 0.2)+
  geom_line(aes(group = subject), alpha = 0.9)+
  facet_wrap(~year, nrow =1)+ theme_classic()+
  xlab ('day of year') + ylab ('soil moisture %')

met_plot <- cowplot::plot_grid(temp_plot_air, temp_plot, moist_plot, nrow =3)

ggsave(met_plot, file = 'plots/met_with_phen.jpg', height = 8, width =7)


# summary anovas ----------------------------------------------------------

#some summary stats

car::Anova(lm (sos~ sensornode + year,
               data = sumry %>% filter (year>2017) %>%
                 mutate (year = factor(year), sensornode = factor (sensornode))),
           type = 'III')

car::Anova(lm (pop ~ sensornode + year,
               data = sumry %>% filter (year>2017) %>%
                 mutate (year = factor(year), sensornode = factor (sensornode))),
           type = 'III')


car::Anova(lm (eos ~ sensornode + year,
               data = sumry %>% filter (year>2017) %>%
                 mutate (year = factor(year), sensornode = factor (sensornode))),
           type = 'III')



# continue surv formatting ------------------------------------------------


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


# full models ---------------------------------------------------------------


#define full models for prediction
# note definitely some things are collinear - but retained bc they explain
# excess variation

#full model
# here filtered on days_since_snow >0 because limited growth prior to snow

#solves the same as this
sos_all = glm(formula = status ~ #GDD +
                snowmelt_doy_infilled+
                  moist_7d + yday + temp_7d,
                family = binomial(link = "logit"),
                data = phen_data_sos %>%
                  ungroup ()%>%
                  filter(days_since_snow>0))


#air
sos_air_all = glm(formula = status ~ GDD_air + yday + snowmelt_doy_infilled +temp_7d_air,
                  family = binomial(link = "logit"),
                  data = phen_data_sos %>%
                    ungroup ()%>%
                    filter(days_since_snow>0))

pop_all = glm(formula = status ~ GDD + snowmelt_doy_infilled + moist_7d + yday + temp_7d,
              family = binomial(link = "logit"),
              data = phen_data_pop %>%
                ungroup ()%>%
                filter(days_since_snow>0))


#can't just put snow in there or get fitted
# probs that are essentially 1 bc it never occurs before snowmelt
# so constrait it by just putting things after snowmelt in
# pop_all_2 = glm(formula = status ~ GDD + snowmelt_doy_infilled + moist_7d + yday + temp_7d + snow,
#               family = binomial(link = "logit"),
#               data = phen_data_pop %>%
#                 ungroup ())


pop_air_all = glm(formula = status ~ GDD_air+ snowmelt_doy_infilled + yday + temp_7d_air,
                  family = binomial(link = "logit"),
                  data = phen_data_pop %>%
                    ungroup %>%
                    filter(days_since_snow>0))

eos_all = glm(formula = status ~ GDD + snowmelt_doy_infilled + moist_7d + yday + temp_7d,
              family = binomial(link = "logit"),
              data = phen_data_eos %>%
                ungroup %>%
                filter(days_since_snow>0))

eos_air_all = glm(formula = status ~ GDD_air + yday + days_since_snow + temp_7d_air,
                  family = binomial(link = "logit"),
                  data = phen_data_eos %>%
                    ungroup %>%
                    filter(days_since_snow>0))


# select variables, stepwise ----------------------------------------------


#backwards stepwise
model.aic.backward_sos <- step(sos_all, direction = "backward", trace = 1)
summary (model.aic.backward_sos)

model.aic.backward_sos_air <- step(sos_air_all, direction = "backward", trace = 1)
summary (model.aic.backward_sos_air)

model.aic.backward_pop <- step(pop_all, direction = "backward", trace = 1)
summary (model.aic.backward_pop) # warmer springs make it earlier, but cooler summers make it later?

model.aic.backward_pop_air <- step(pop_air_all, direction = "backward", trace = 1)
summary (model.aic.backward_pop_air) # warmer springs make it earlier, but cooler summers make it later?

model.aic.backward_eos <- step(eos_all, direction = "backward", trace = 1)
summary (model.aic.backward_eos)

model.aic.backward_eos_air <- step(eos_air_all, direction = "backward", trace = 1)
summary (model.aic.backward_eos_air)


# define functions for validation -----------------------------------------


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


# compare model fits EOS --------------------------------------------------


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

#days_since_snow + moist_7d + yday + temp_7d

#surv wo days_since_snow
eos_fitted_drop_snowmelt_doy_infilled = glm(formula = update(formula(model.aic.backward_eos), ~. -snowmelt_doy_infilled),
                                      family = binomial(link = "logit"),
                                      data = phen_data_eos %>%
                                        ungroup %>%
                                        filter(days_since_snow>0))

eos_drop_snowmelt_doy_infilled_rmse = model_cv_surv(phen_data_eos, phen_data_eos_pred,
                                                    eos_fitted_drop_snowmelt_doy_infilled )
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

#repeat for air #GDD_air + yday
#surv wo yday
#GDD_air + yday
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

eos_full_rmse # 
eos_drop_yday_rmse # 
eos_drop_temp_7d_rmse # 
eos_drop_moist_7d_rmse # 
eos_drop_snowmelt_rmse #
eos_rmse_null # 
eos_rmse_GDD
eos_rmse_GDD_air
eos_rmse_days_since_snow

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
  'eos_drop_snowmelt_rmse', 
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
  #relable for nice plots
  mutate(
    mod = case_when(
      mod == 'eos_full_rmse' ~ 'time-to-event (soil temperature)',
      mod == 'eos_full_rmse_air' ~ 'time-to-event (air temperature)',
      mod == 'eos_rmse_days_since_snow' ~ 'Threshold (days past snowmelt)',
      mod == 'eos_rmse_null' ~ 'Threshold (day of year)',
      mod == 'eos_rmse_GDD' ~ 'Threshold (GDD soil)',
      mod == 'eos_rmse_GDD_air' ~ 'Threshold (GDD air)',
      mod == 'eos_drop_yday_rmse' ~ '- day of year',
      mod == 'eos_drop_temp_7d_rmse' ~ '- 7d running mean soil temperature',
      mod == 'eos_drop_moist_7d_rmse' ~ '- 7d running mean soil moisture',
      mod == 'eos_drop_snowmelt_rmse' ~ '- snowmelt',
      mod == 'eos_drop_yday_rmse_air' ~ '- day of year (air temperature model)',
      mod == 'eos_drop_GDD_rmse_air' ~ '- GDD (air temperature model)',
      TRUE ~ mod
    )
  )

eos_vs_thresh <- ggplot(all_rmse_eos %>%
                          mutate (mod = gsub('\\(', '\n(', mod))%>%
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
  scale_fill_grey()+ guides(fill=guide_legend(title="holdout set"))+
  ylab('RMSE')



# plotting functions ------------------------------------------------------


#define function to plot effects
plot_conditional_effect = function(dataset,
                                   phenometric = 'greenup', model,
                                   all_rmse, soil_or_air = 'soil', linecolor = 'green'
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
    }else if (term == 'snowmelt_doy_infilled'){
      'snowmelt doy'
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
      out[[v]] = mean (avg_values[[v]], na.rm = TRUE)
    }
    preds = data.frame(dataset%>%
                         filter(yday>=min_event&yday<max_event)%>%
                         select(!!term))%>%
      dplyr::filter(!is.na(!!sym(term)))%>%
      distinct()
    for (v in vars){
      preds[[v]] = out[[v]]
    }
    preds$hazard = predict(model, preds, type = 'response')
    
    p1 <- ggplot (preds, aes_string(x=term, y='hazard'))+
      geom_line(color = linecolor, size =2)+
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
      #ylim(-1,20)+
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
  'eos_full_rmse', 
  'eos_full_rmse_air',
  'eos_drop_yday_rmse', 
  'eos_drop_yday_rmse_air', 
  'eos_drop_GDD_rmse_air', 
  'eos_drop_temp_7d_rmse', 
  'eos_drop_moist_7d_rmse', 
  'eos_drop_snowmelt_doy_infilled_rmse', 
  'eos_rmse_null', 
  'eos_rmse_GDD',
  'eos_rmse_GDD_air',
  'eos_rmse_days_since_snow')){
  all_rmse_eos[[mod]] = as.data.frame(get(mod))
}

all_rmse_eos = all_rmse_eos%>%
  data.table::rbindlist(., idcol = 'model')

#all_rmse_sos = all_rmse_sos %>% data.table::rbindlist(, idcol = 'model')
eos_effects <- plot_conditional_effect(dataset =phen_data_eos,
                                       model =eos_fitted,
                                       all_rmse = all_rmse_eos,
                                       soil_or_air = 'soil', linecolor = 'brown')

p1 = cowplot::plot_grid(plotlist =eos_effects)

x.grob <- textGrob("senescence", 
                   gp=gpar(col="black", fontsize=20))
y.grob <- textGrob("h(x)", 
                   gp=gpar(col="black", fontsize=15), rot=90)

eos_all<-grid.arrange(arrangeGrob(p1, left = y.grob, top = x.grob))
jpeg('plots/eos_all.jpg', width =600, height =600)
ggdraw(eos_all)
dev.off()


# compare model fits peak (pop) -------------------------------------------

# repeat for peak#####
#drop each in turn
#GDD + snowmelt_doy_infilled + moist_7d + yday + temp_7d
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

#GDD + days_since_snow + moist_7d + yday + temp_7d

pop_fitted_drop_GDD = glm(formula = update(formula(model.aic.backward_pop), ~. -GDD),
                          family = binomial(link = "logit"),
                          data = phen_data_pop %>%
                            ungroup %>%
                            filter(days_since_snow>0))

pop_drop_GDD_rmse = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                  pop_fitted_drop_GDD)


pop_fitted_drop_snowmelt_doy_infilled = glm(formula = update(formula(model.aic.backward_pop), ~. -snowmelt_doy_infilled),
                                            family = binomial(link = "logit"),
                                            data = phen_data_pop %>%
                                              ungroup %>%
                                              filter(days_since_snow>0))

pop_drop_snowmelt_doy_infilled_rmse = model_cv_surv(phen_data_pop, phen_data_pop_pred,
                                                    pop_fitted_drop_snowmelt_doy_infilled)


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

pop_fitted_drop_snowmelt_doy_infilled_air= glm(formula = update(formula(model.aic.backward_pop_air), ~.
                                                          -snowmelt_doy_infilled),
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
  'pop_drop_snowmelt_doy_infilled_rmse', 
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
  #relable for nice plots
  mutate(
    mod = case_when(
      mod == 'pop_full_rmse' ~ 'time-to-event (soil temperature)',
      mod == 'pop_full_rmse_air' ~ 'time-to-event  (air temperature)',
      mod == 'pop_rmse_days_since_snow' ~ 'Threshold (days past snowmelt)',
      mod == 'pop_rmse_null' ~ 'Threshold (day of year)',
      mod == 'pop_rmse_GDD' ~ 'Threshold (GDD soil)',
      mod == 'pop_rmse_GDD_air' ~ 'Threshold (GDD air)',
      mod == 'pop_drop_yday_rmse' ~ '- day of year',
      mod == 'pop_drop_GDD_rmse' ~ '- GDD (soil)',
      mod == 'pop_drop_GDD_air_rmse_air' ~ '- GDD (air)',
      mod == 'pop_drop_temp_7d_rmse' ~ '- 7d running mean soil temperature',
      mod == 'pop_drop_moist_7d_rmse' ~ '- 7d running mean soil moisture',
      mod == 'pop_drop_snowmelt_doy_infilled_rmse' ~ '- snowmelt',
      mod == 'pop_drop_yday_rmse_air' ~ '- day of year (air temperature model)',
      mod == 'pop_drop_temp_7d_rmse_air' ~ '- 7d running mean air temperature',
      mod == 'pop_drop_snowmelt_doy_infilled_rmse_air' ~ '- snowmelt (air temperature model)',
      TRUE ~ mod
    )
  )

pop_vs_thresh <- ggplot(all_rmse_pop %>%
                          mutate (mod = gsub('\\(', '\n(', mod))%>%
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
  scale_fill_grey()+ guides(fill=guide_legend(title="holdout set"))+
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
  'pop_drop_snowmelt_doy_infilled_rmse', 
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
  data.table::rbindlist(., idcol = 'model')
.
#all_rmse_sos = all_rmse_sos %>% data.table::rbindlist(, idcol = 'model')
pop_effects <- plot_conditional_effect(dataset =phen_data_pop,
                                       model =pop_fitted,
                                       all_rmse = all_rmse_pop,
                                       soil_or_air = 'soil',
                                       linecolor = 'blue')

p1 = cowplot::plot_grid(plotlist =pop_effects, nrow =1)


x.grob <- textGrob("peak", 
                   gp=gpar(col="black", fontsize=20))
y.grob <- textGrob("h(x)", 
                   gp=gpar(col="black", fontsize=15), rot=90)

pop_all<-grid.arrange(arrangeGrob(p1, left = y.grob, top = x.grob))
jpeg('plots/pop_all.jpg', width =600, height =600)
ggdraw(pop_all)
dev.off()

###end repeat for peak

# repeat for SOS ----------------------------------------------------------

## run for sos

sos_fitted = glm(formula = formula(model.aic.backward_sos),
                 family = binomial(link = "logit"),
                 data = phen_data_sos %>%
                   ungroup %>%
                   filter(days_since_snow>0))


#both spatial and temporal predictions improved by adding yday to the model
sos_full_rmse = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                              sos_fitted)

sos_fitted_air = glm(formula = formula(model.aic.backward_sos_air),
                     family = binomial(link = "logit"),
                     data = phen_data_sos %>%
                       ungroup %>%
                       filter(days_since_snow>0))

#both spatial and temporal predictions improved by adding yday to the model
sos_full_rmse_air = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                  sos_fitted_air)

#GDD + days_since_snow + yday + temp_7d
sos_fitted_drop_GDD = glm(formula = update(formula(model.aic.backward_sos), ~. -GDD),
                          family = binomial(link = "logit"),
                          data = phen_data_sos %>%
                            ungroup %>%
                            filter(days_since_snow>0))

sos_drop_GDD_rmse = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                  sos_fitted_drop_GDD)
sos_fitted_snowmelt_doy_infilled = glm(formula = update(formula(model.aic.backward_sos), ~. -snowmelt_doy_infilled),
                                      family = binomial(link = "logit"),
                                      data = phen_data_sos %>%
                                        ungroup %>%
                                        filter(days_since_snow>0))

sos_drop_snowmelt_doy_infilled_rmse = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                              sos_fitted_drop_days_since_snow)

sos_fitted_drop_snowmelt_doy_infilled_air = glm(formula = update(formula(model.aic.backward_sos_air ), ~. -snowmelt_doy_infilled),
                                          family = binomial(link = "logit"),
                                          data = phen_data_sos %>%
                                            ungroup %>%
                                            filter(days_since_snow>0))

sos_drop_snowmelt_doy_infilled_rmse_air = model_cv_surv(phen_data_sos, phen_data_sos_pred,
                                                  sos_fitted_drop_snowmelt_doy_infilled_air)


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
  'sos_drop_snowmelt_doy_infilled_rmse',
  'sos_drop_temp_7d_rmse', 
  'sos_rmse_null', 
  'sos_rmse_GDD',
  'sos_rmse_GDD_air',
  'sos_rmse_days_since_snow',
  'sos_full_rmse_air', 
  'sos_drop_yday_rmse_air', 
  'sos_drop_temp_7d_rmse_air', 
  'sos_drop_snowmelt_doy_infilled_rmse_air')){ 
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
      mod == 'sos_full_rmse' ~ 'time-to-event (soil temperature)',
      mod == 'sos_full_rmse_air' ~ 'time-to-event (air temperature)',
      mod == 'sos_rmse_days_since_snow' ~ 'Threshold (days past snowmelt)',
      mod == 'sos_rmse_null' ~ 'Threshold (day of year)',
      mod == 'sos_rmse_GDD' ~ 'Threshold (GDD soil)',
      mod == 'sos_rmse_GDD_air' ~ 'Threshold (GDD air)',
      mod == 'sos_drop_yday_rmse' ~ '- day of year',
      mod == 'sos_drop_GDD_rmse' ~ '- GDD (soil temperature model)',
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
                          mutate (mod = gsub('\\(', '\n(', mod))%>%
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
  scale_fill_grey()+ guides(fill=guide_legend(title="holdout set"))+
  ylab('RMSE')

# conditional effects plots

min_sos = min(phen_data_sos$doy_event)
max_sos = max(phen_data_sos$doy_event)

#regrab the rmses to make deltas
all_rmse_sos = list()
for (mod in c(
  'sos_full_rmse', 
  'sos_drop_yday_rmse', 
  'sos_drop_GDD_rmse', 
  'sos_drop_snowmelt_doy_infilled_rmse', 
  'sos_drop_temp_7d_rmse', 
  'sos_rmse_null', 
  'sos_rmse_GDD',
  'sos_rmse_GDD_air',
  'sos_rmse_days_since_snow',
  'sos_full_rmse_air', 
  'sos_drop_yday_rmse_air', 
  'sos_drop_temp_7d_rmse_air', 
  'sos_drop_snowmelt_doy_infilled_rmse_air')){ 
  all_rmse_sos[[mod]] = as.data.frame(get(mod))
}

all_rmse_sos = all_rmse_sos %>% data.table::rbindlist(., idcol = 'model')

#all_rmse_sos = all_rmse_sos %>% data.table::rbindlist(, idcol = 'model')
sos_effects <- plot_conditional_effect(dataset =phen_data_sos,
                                       model =sos_fitted,
                                       all_rmse = all_rmse_sos,
                                       soil_or_air = 'soil',
                                       linecolor = 'green')

p1 = cowplot::plot_grid(plotlist =sos_effects)

x.grob <- textGrob("greenup", 
                   gp=gpar(col="black", fontsize=20))
y.grob <- textGrob("h(x)", 
                   gp=gpar(col="black", fontsize=15), rot=90)

sos_all<-grid.arrange(arrangeGrob(p1, left = y.grob, top = x.grob))
jpeg('plots/sos_all.jpg', width =600, height =600)
ggdraw(sos_all)
dev.off()

#####

# multiplots --------------------------------------------------------------

#v1
#all together
# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  eos_vs_thresh + theme(legend.box.margin = margin(0, 0, 0, 12),
                        legend.text = element_text(size=14),
                        legend.title= element_text(size=20))
)


prow <- cowplot::plot_grid(sos_vs_thresh+theme(legend.position = 'none',axis.title.y = element_blank(),
                                               axis.text.x = element_text(size = 14),
                                               plot.title = element_text(size = 25)),
                           pop_vs_thresh+theme(legend.position = 'none', axis.title.y = element_blank(),
                                               axis.text.x = element_text(size = 14),
                                               plot.title = element_text(size = 25)),
                           eos_vs_thresh+theme(legend.position = 'none', axis.title.y = element_blank(),
                                               axis.text.x = element_text(size = 14),
                                               plot.title = element_text(size = 25)), nrow = 1)

y.grob <- textGrob("\nRMSE", 
                   gp=gpar(col="black", fontsize=20), rot=90)
p1 = plot_grid(prow, legend, rel_widths = c(3,.4))

rmse_all<-grid.arrange(arrangeGrob(p1, left = y.grob))

jpeg('plots/rmse_all.jpg', width =1000, height =600)
ggdraw(rmse_all)
dev.off()

#rearrange
all_plots <-
  sos_effects[["snowmelt_doy_infilled"]] + 
  sos_effects[["temp_7d"]] +
  #sos_effects[["GDD"]]+theme(axis.title.x = element_blank()) +
  plot_spacer() + # moisture
  sos_effects[["yday"]] + 
  plot_spacer() + #GDD
  pop_effects[["snowmelt_doy_infilled"]] + 
  pop_effects[["temp_7d"]] +
  pop_effects[["GDD"]]+
  pop_effects[["yday"]] + 
  pop_effects[["moist_7d"]] + 
  eos_effects[["snowmelt_doy_infilled"]] + 
  eos_effects[["temp_7d"]] +
  plot_spacer() + # GDD
  eos_effects[["yday"]] + 
  eos_effects[["moist_7d"]] + plot_layout(ncol=5)
  
  


# 
# 
# #orig
# all_plots <-
#   #sos_effects[["GDD"]]+theme(axis.title.x = element_blank()) +
#   plot_spacer() +
#   sos_effects[["snowmelt_doy_infilled"]] + 
#   plot_spacer() +
#   sos_effects[["yday"]] + 
#   sos_effects[["temp_7d"]] +
#   pop_effects[["GDD"]] + pop_effects[["snowmelt_doy_infilled"]] + 
#   pop_effects[["moist_7d"]]  + 
#   pop_effects[["yday"]] + 
#   pop_effects[["temp_7d"]] + 
#   #plot_spacer() +
#   sos_effects[["GDD"]]+geom_rect(aes(xmin =0.05, xmax = 650,
#                                      ymin =0.1, ymax =1), 
#                                  alpha = 1, fill = 'white')+
#   eos_effects[["snowmelt_doy_infilled"]]+
#   eos_effects[["moist_7d"]]  + 
#   eos_effects[["yday"]] + 
#   eos_effects[["temp_7d"]] + plot_layout(ncol=5)


y.grob <- textGrob("h(x)", 
                   gp=gpar(col="black", fontsize=15), rot=90)

gt <- patchwork::patchworkGrob(all_plots )

all_plots<-grid.arrange(arrangeGrob(gt , left = y.grob))
jpeg('plots/all_effects.jpg', width =700, height =600, quality =100)
ggdraw(all_plots)
dev.off()


#sce test this thing
###
all_mods = list(
  sos=
    data.frame(summary (model.aic.backward_sos)$coefficients)%>%
    mutate(param = row.names(.)),
  pop = data.frame(summary (model.aic.backward_pop)$coefficients)%>%
    mutate(param = row.names(.)),
  eos = data.frame(summary (model.aic.backward_eos)$coefficients)%>%
    mutate(param = row.names(.))
) %>%
  data.table::rbindlist(., id = 'phenophase')

ggplot (broom::tidy(model.aic.backward_sos, conf.int = TRUE) %>%
          #           #rename(phenophase = name) %>%
          #           filter(term !='(Intercept)'),
          #         aes(xmin=conf.low, xmax = conf.high, y=term))+#, color = spec))+
          #   geom_vline(aes(xintercept =0))+
          #   geom_point(aes(x =estimate))+
          #   geom_linerange(position=position_dodge(width=0.5))
          
          tst<-all_mods %>%
          mutate(ci_lower = Estimate - 2*Std..Error,
                 ci_upper = Estimate + 2*Std..Error)
        
        
        ggplot (all_mods %>%
                  filter(param %in%c('GDD', 'moist_7d', 'temp_7d',
                                     'days_since_snow'))%>%
                  mutate(phenophase =
                           case_when(phenophase == 'sos' ~ 'greenup',
                                     phenophase == 'pop' ~ 'peak',
                                     phenophase == 'eos' ~ 'senescence'))%>%
                  mutate(ci_lower = Estimate - 2*Std..Error,
                         ci_upper = Estimate + 2*Std..Error),
                aes(x=phenophase,
                    y=Estimate, color = phenophase))+
          geom_line(aes(group =1), color = 'black')+
          geom_linerange(aes(
            ymin = ci_lower,
            ymax= ci_upper))+
          geom_point()+
          geom_hline(aes(yintercept =0), linetype = 'dotted')+
          facet_wrap(~param, scales = 'free_y')+
          ylab ('B \n delay <-> advance')+
          scale_colour_manual(values = c(greenup = "green", peak = "blue",
                                         senescence = "brown"))+
          theme_classic()
        
        

# all junk afte here ------------------------------------------------------

        # eos_effects <- plot_conditional_effect(dataset =phen_data_sos,
        #                                        model =sos2_fitted,
        #                                        all_rmse = all_rmse_sos,
        #                                        soil_or_air = 'soil', linecolor = 'green')
        
        
        
        #model.aic.backward_sos_2 <- step(sos_all_2, direction = "backward", trace = 1)
        #summary (model.aic.backward_sos_2)
        #sos_all_2
        
        
        #sos2_fitted
        
        #this is a start
        # if we want to do it
        # though I have not independently
        # recalcuated the runnigng means
        #forgot to make the cumsums
        
        mod = update (model, data = phen_data_fit_i %>%
                        filter(days_since_snow>0))
        
        dataset = phen_data_sos_pred
        surv  = 
          phen_data_sos_pred %>%
          ungroup %>%
          filter(days_since_snow>0) %>%
          mutate(hazard = c(predict(sos2_fitted, newdata =.,
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
          ungroup()%>%
          summarize(pred = mean (doy_event_pred ))
        
        surv_early_snow = phen_data_sos_pred %>%
          ungroup %>%
          mutate (days_since_snow = days_since_snow+18,
                  snowmelt_doy_infilled = snowmelt_doy_infilled-18,
                  sf = if_else (yday>snowmelt_doy_infilled,1,0))%>%
          #use air for soil when not sf
          mutate(soiltemp_5cm_avg_ave = if_else(sf ==1, soiltemp_5cm_avg_ave, airtemp_avg_ave))%>%
          arrange(sensornode, year, yday)%>%
          group_by(sensornode, year) %>%
          mutate(DD = ifelse(soiltemp_5cm_avg_ave>0, soiltemp_5cm_avg_ave, 0),
                 GDD = cumsum(DD))%>%
          ungroup() %>%
          filter(days_since_snow>0) %>%
          mutate(hazard = c(predict(sos2_fitted, newdata =.,
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
          ungroup()%>%
          summarize(pred = mean (doy_event_pred ))
        
        
        
        surv_late_snow = phen_data_sos_pred %>%
          ungroup %>%
          #use 0 for soil when snow covered
          mutate(soiltemp_5cm_avg_ave = if_else((yday > snowmelt_doy_infilled & yday<=snowmelt_doy_infilled+18),
                                                0, soiltemp_5cm_avg_ave))%>%
          mutate (days_since_snow = days_since_snow-18,
                  snowmelt_doy_infilled = snowmelt_doy_infilled+18,
                  sf = if_else (yday>snowmelt_doy_infilled,1,0))%>%
          arrange(sensornode, year, yday)%>%
          group_by(sensornode, year) %>%
          mutate(DD = ifelse(soiltemp_5cm_avg_ave>0, soiltemp_5cm_avg_ave, 0),
                 GDD = cumsum(DD))%>%
          ungroup() %>%
          filter(days_since_snow>0) %>%
          mutate(hazard = c(predict(sos2_fitted, newdata =.,
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
          ungroup()%>%
          summarize(pred = mean (doy_event_pred ))
        
        
        
        
        
        
        
        
        
        
#I think this is all junk below here?

library (gridExtra)
library(grid)
y.grob <- textGrob("daily probability of senescence", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)

grid.arrange(arrangeGrob(cowplot::plot_grid(plotlist = eos_effects), 
                         left =y.grob))

#add to plot

jpg('pop_effect.jpg')
pop_effects = grid.arrange(arrangeGrob(cowplot::plot_grid(plotlist =pop_effects_panels), 
                                       left =y.grob))








pop_effects_panels <- plot_conditional_effect(dataset =phen_data_eos,
                                              phenometric = 'peak', model =pop_fitted)

cowplot::plot_grid(plotlist =pop_effects)

library (gridExtra)
library(grid)
y.grob <- textGrob("daily probability of greenup", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)

x.grob <- textGrob("Common X", 
                   gp=gpar(fontface="bold", col="blue", fontsize=15))

#add to plot
grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))

jpg('pop_effect.jpg')
pop_effects = grid.arrange(arrangeGrob(cowplot::plot_grid(plotlist =pop_effects_panels), 
                                       left =y.grob))




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
all_rmse_sos_delta = all_rmse_sos %>%
  mutate(full_spatial = all_rmse_sos$spatial[all_rmse_sos$model == 'sos_full_rmse'],
         full_temporal = all_rmse_sos$temporal[all_rmse_sos$model == 'sos_full_rmse'])%>%
  mutate(spatial = spatial - full_spatial,
         temporal = temporal - full_temporal) %>%
  select(model, spatial, temporal) %>%
  pivot_longer(cols = c(spatial, temporal),
               names_to = 'type',
               values_to = 'delta_rmse')

inset.plot<-ggplot(all_rmse_sos_delta  %>%
                     filter(model == 'sos_drop_GDD_rmse'), aes(fill=type, y=delta_rmse, x=model))+
  # geom_hline(aes(yintercept = all_rmse_sos$rmse[all_rmse_sos$mod=='best model (soil temperature)' & all_rmse_sos$type == 'spatial']),
  #            color = 'black', alpha = 0.5)+
  # geom_hline(aes(yintercept = all_rmse_sos$rmse[all_rmse_sos$mod=='best model (soil temperature)' & all_rmse_sos$type == 'temporal']),
  #            color = 'grey', alpha = 0.5) + 
  geom_bar(position="dodge", stat="identity", color = 'black')+
  theme(axis.text.x = element_text(angle = 90)) +
  #scale_fill_manual(values=c("black", "grey"))+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  #theme(axis.text.x = element_text(angle = 45, hjust =1))+
  #theme(axis.text.x = NULL)+
  scale_fill_grey()+ #guides(fill=guide_legend(title="test set"))+
  #ggtitle ('sos_drivers_soil')+
  ylim(-1,1)+
  theme(legend.position = 'none')



main.plot <- 
  ggplot(data = mpg, aes(x = cty, y = hwy, colour = factor(cyl))) + 
  geom_point(size = 2.5)

inset.plot <- main.plot + theme(legend.position = "none")



plot.with.inset <-
  ggdraw() +
  draw_plot(sos_effects[[1]]) +
  draw_plot(inset.plot, x = 0.5, y = .7, width = .5, height = .3)


cowplot::plot_grid(plotlist =list(plot.with.inset, plot.with.inset2))


#dataset = phen_data_sos
#min_event = min_sos
#max_event = max_sos
#model = sos_fitted
#term = 'GDD'
sos_effects <- plot_conditional_effect(dataset =phen_data_sos,
                                       phenometric = 'greenup', model =sos_fitted,
                                       linecolor = 'green')

cowplot::plot_grid(plotlist =list(plot.with.inset, plot.with.inset2))

eos_effects <- plot_conditional_effect(dataset =phen_data_eos,
                                       phenometric = 'senescence', model =eos_fitted)

cowplot::plot_grid(plotlist =eos_effects)

pop_effects_panels <- plot_conditional_effect(dataset =phen_data_eos,
                                              phenometric = 'peak', model =pop_fitted)

cowplot::plot_grid(plotlist =pop_effects)




# avg_values_sos = phen_data_sos%>%
#   filter(yday>=min_sos&yday<max_sos)%>%
#   summarize(GDD = mean(GDD), days_since_snow =mean(days_since_snow),
#             temp_7d = mean(temp_7d),
#             yday = mean (yday))

#

avg_values_sos = phen_data_sos%>%
  filter(yday>=min_sos&yday<max_sos)%>%
  summarize(GDD = mean(GDD), days_since_snow =mean(days_since_snow),
            temp_7d = mean(temp_7d),
            yday = mean (yday))

GDD_preds = data.frame(GDD = 
                         phen_data_sos%>%
                         filter(yday>=min_sos&yday<max_sos)%>%
                         select(GDD),
                       yday = avg_values_sos$yday,
                       temp_7d = avg_values_sos$temp_7d,
                       days_since_snow = avg_values_sos$days_since_snow)

GDD_preds$hazard = predict(sos_fitted, GDD_preds, type = 'response')

ggplot (GDD_preds, aes(x=GDD, y=hazard))+
  geom_line()+
  ylab ('increase in daily probability of greenup')+
  theme_classic()+ylim(0,0.5)

days_since_snow_preds = data.frame(days_since_snow = 
                                     phen_data_sos%>%
                                     filter(yday>=min_sos&yday<max_sos)%>%
                                     select(days_since_snow),
                                   yday = avg_values_sos$yday,
                                   temp_7d = avg_values_sos$temp_7d,
                                   GDD = avg_values_sos$GDD)

days_since_snow_preds$hazard = predict(sos_fitted, days_since_snow_preds, type = 'response')

#https://cran.r-project.org/web/packages/ggeffects/ggeffects.pdf
ggplot (days_since_snow_preds, aes(x=days_since_snow, y=hazard))+
  geom_line()+
  ylab ('increase in daily probability of greenup')+
  theme_classic()+ylim(0,0.5)

plot(effects::allEffects(sos_fitted,
                         fixed.predictors =
                           list(given.values=
                                  c(yday=mean(phen_data_sos$yday[phen_data_sos$yday>=min_sos&phen_data_sos$yday<=max_sos]), 
                                    GDD = mean(phen_data_sos$GDD[phen_data_sos$yday>=min_sos&phen_data_sos$yday<=max_sos]))
                           )),
     #xlevels = dlist,
     confint = FALSE,
     type = 'response',
     rug = FALSE,
     axes=list(
       grid=TRUE,
       y=list(rotate=90,
              lim = c(0,0.8),
              lab=c("hazard")
              
       ),
       x=list(
         rotate=30,
         yday=list(lab="day of year",
                   lim=c(min_sos, max_sos)),
         GDD = list(lab = 'cumulative GDD',
                    lim=c(min(phen_data_sos_pred$GDD[phen_data_sos_pred$yday>=min_sos&phen_data_sos_pred$yday<=max_sos]),
                          max(phen_data_sos_pred$GDD[phen_data_sos_pred$yday>=min_sos&phen_data_sos_pred$yday<=max_sos]))))
     )
)

ggplot ()


###end repeat for sos

## run for sos


#end repeat for sos

#sce should compare these to 2 null models,
# first is GDD_event ~ 1
# 2nd is yday_event ~ 1

simplest_mod = lm (yday ~ 1,
                   data = phen_data_pop %>%
                     ungroup %>%
                     filter(yday == doy_event))

#both improve a lot over the null model of yday
# t2 = model_cv_simplest(phen_data_sos, phen_data_sos_pred, simplest_mod)

#try to get things onto laptop to use lasso?
# tst =phen_data_sos %>%
#   ungroup %>%
#   filter(yday == doy_event)
# 
# write.csv(tst, "tst.csv")

model_cv_simplest = function(
    phen_data_fit, phen_data_pred, model){
  out = list()
  for (i in unique (phen_data_fit$sensornode)){
    phen_data_fit_i = phen_data_fit %>%
      filter(sensornode != i) %>%
      ungroup %>%
      filter(yday == doy_event)
    phen_data_pred_i = phen_data_pred %>%
      filter(sensornode == i)%>%
      ungroup %>%
      filter(yday == doy_event)
    mod = update (model, data = phen_data_fit_i)
    # surv  = predict(mod, newdata =phen_data_pred_i,
    #                 type = 'response')
    surv = phen_data_pred_i %>%
      mutate(doy_event_pred = c(predict(mod, newdata =.,
                                        type = 'response'))) %>%
      select(subject, doy_event_pred, year) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, doy_event)%>%
                  distinct(),by = c("subject", "year"))
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
      filter(year == i)%>%
      ungroup %>%
      filter(yday == doy_event)
    mod = update (model, data = phen_data_fit_i)
    # surv  = predict(mod, newdata =phen_data_pred_i,
    #                 type = 'response')
    surv = phen_data_pred_i %>%
      mutate(doy_event_pred = c(predict(mod, newdata =.,
                                        type = 'response'))) %>%
      select(subject, doy_event_pred, year) %>%
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

#both improve a lot over the null model of yday
model_cv_simplest(phen_data_sos, phen_data_sos_pred, simplest_mod)

simple_mod = lm (GDD ~ 1,
                 data = phen_data_sos %>%
                   ungroup %>%
                   filter(yday == doy_event))

t3 =model_cv_simple(phen_data_sos, phen_data_sos_pred, simple_mod)

#t1 gets the temporal better than t3, though at a cost of misspecifying the
# spatial a little bit?

#t2 sucks in any regard
# t2

model_cv_simple = function(
    phen_data_fit, phen_data_pred, model){
  out = list()
  for (i in unique (phen_data_fit$sensornode)){
    phen_data_fit_i = phen_data_fit %>%
      filter(sensornode != i) %>%
      ungroup %>%
      filter(yday == doy_event)
    phen_data_pred_i = phen_data_pred %>%
      filter(sensornode == i) #%>%
    #select(subject, year) %>%
    #distinct()
    mod = update (model, data = phen_data_fit_i)
    # surv  = predict(mod, newdata =phen_data_pred_i,
    #                 type = 'response')
    surv = phen_data_pred_i %>%
      select(subject, year) %>%
      distinct() %>%
      mutate(GDD_pred = c(predict(mod, newdata =.,
                                  type = 'response'))) %>%
      select(subject, GDD_pred, year) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, GDD, doy_event, yday)%>%
                  distinct(),
                by = c("subject", "year")) %>%
      filter(GDD>GDD_pred) %>%
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
      filter(year == i) #%>%
    #select(subject, year) %>%
    #distinct()
    mod = update (model, data = phen_data_fit_i)
    # surv  = predict(mod, newdata =phen_data_pred_i,
    #                 type = 'response')
    surv = phen_data_pred_i %>%
      select(subject, year) %>%
      distinct() %>%
      mutate(GDD_pred = c(predict(mod, newdata =.,
                                  type = 'response'))) %>%
      select(subject, GDD_pred, year) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, GDD, doy_event, yday)%>%
                  distinct(),
                by = c("subject", "year")) %>%
      filter(GDD>GDD_pred) %>%
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


#can remake these fits
# for the others, too

simple_mod = lm (GDD ~ 1,
                 data = phen_data_pop %>%
                   ungroup %>%
                   filter(yday == doy_event))

pop_fitted = glm(formula = status ~ GDD + yday + temp_7d,
                 family = binomial(link = "logit"),
                 data = phen_data_pop%>%
                   ungroup %>%
                   filter(days_since_snow>0))

pop_fitted_2 = glm(formula = status ~ GDD + I(GDD^2) + yday + temp_7d,
                   family = binomial(link = "logit"),
                   data = phen_data_pop%>%
                     ungroup %>%
                     filter(days_since_snow>0))

#these help in both cross valid
t3 =model_cv_simple(phen_data_pop, phen_data_pop_pred, simple_mod)
t2 =model_cv_simplest(phen_data_pop, phen_data_pop_pred, simplest_mod)
t1 =model_cv(phen_data_pop, phen_data_pop_pred, pop_fitted)
t0 =model_cv(phen_data_pop, phen_data_pop_pred, pop_fitted_2)

plot(effects::allEffects(pop_fitted_2), type = 'response',
     rug = FALSE,
     axes=list(
       grid=TRUE,
       y=list(rotate=30,
              lim = c(0,0.4)
       )))

plot(effects::allEffects(pop_fitted), type = 'response',
     rug = FALSE,
     axes=list(
       grid=TRUE,
       y=list(rotate=30,
              lim = c(0,0.4)
       )))


#I think this fixes
plot(effects::allEffects(Gompertz_Model_Baseline5_eos_wtemp_3,
                         fixed.predictors = list(given.values=c(yday=200, 
                                                                days_since_snow = 60))),
     type = 'response',
     rug = FALSE,
     axes=list(
       grid=TRUE,
       y=list(rotate=30,
              lim = c(0,0.4)
       )))



eos_fitted = glm(formula = status ~ days_since_snow + yday + temp_7d + moist_7d,
                 family = binomial(link = "logit"),
                 data = phen_data_eos%>%
                   ungroup %>%
                   filter(days_since_snow>0))

eos_fitted_2 = glm(formula = status ~ I(yday^2) +I(temp_7d^2) +
                     days_since_snow + yday + temp_7d + I(moist_7d^2),
                   family = binomial(link = "logit"),
                   data = phen_data_eos%>%
                     ungroup %>%
                     filter(days_since_snow>0))

#if you care to look at the 2 fitteds,
# the 2nd one makes the temp effect come on more steeply
# the yday - I dunno looks similar
plot(effects::allEffects(eos_fitted_2), type = 'response')
plot(effects::allEffects(eos_fitted), type = 'response')

,
rug = FALSE,
axes=list(
  grid=TRUE,
  y=list(rotate=30,
         lim = c(0,0.4)
  )))


# eos_fitted_2 = glm(formula = status ~ GDD + yday + temp_7d,
#                  family = binomial(link = "logit"),
#                  data = phen_data_eos%>%
#                    ungroup %>%
#                    filter(yday>180))


t_gdd =model_cv_simple(phen_data_eos, phen_data_eos_pred, simple_mod)

#nope, something weird about the eos model
# that is overfitting vis-a-vis the
# start of season things.
# maybe bc it's starting at snowmelt?

#nope, it's weirdnesses
t_yday =model_cv_simplest(phen_data_eos, phen_data_eos_pred, simplest_mod)

#best mod - hooray
t_surv =model_cv(phen_data_eos, phen_data_eos_pred, eos_fitted)

#adding in the quadratics don't really help - just checking for thresholds?
eos_fitted_2 = glm(formula = status ~ days_since_snow + yday + I(yday^2) + temp_7d + I(temp_7d^2)+ moist_7d + I(moist_7d^2),
                   family = binomial(link = "logit"),
                   data = phen_data_eos%>%
                     ungroup %>%
                     filter(days_since_snow>0))


## could mess with nonlinear functs if desired but
## actually can't get this to load
#https://rdrr.io/cran/gnlm/man/bnlr.html


simple_mod = lm (GDD ~ 1,
                 data = phen_data_sos %>%
                   ungroup %>%
                   filter(yday == doy_event))

sos_fitted = glm(formula = status ~ GDD + yday,
                 family = binomial(link = "logit"),
                 data = phen_data_sos %>%
                   ungroup %>%
                   filter(days_since_snow>0))
phen_data_fit = phen_data_sos
phen_data_pred = phen_data_sos_pred
model = sos_fitted
model_cv = function(
    phen_data_fit, phen_data_pred, space_or_time, model
)
  
  
  pred_trans = PredSurv %>%
  group_by(subject) %>%
  arrange(yday) %>%
  filter(failure<0.5) %>%
  slice(1) %>%
  select(subject, yday)

sos_fitted = glm(formula = status ~ GDD + yday,
                 family = binomial(link = "logit"),
                 data = phen_data_sos %>%
                   ungroup %>%
                   filter(days_since_snow>0))




tst = predict(sos_fitted, new_data = phen_data_sos_pred %>%
                ungroup %>%
                filter(days_since_snow>0),
              type = 'response')

tst  = 
  phen_data_sos_pred %>%
  ungroup %>%
  filter(days_since_snow>0) %>%
  mutate(hazard = c(predict(sos_fitted, newdata =.,
                            type = 'response'))) %>%
  mutate(haz_compl = 1-hazard) %>%
  group_by(subject) %>%
  arrange(yday) %>%
  mutate(surv = cumprod(haz_compl))


#could make a drop-1-site
# and drop-1-year bakeoff?




#these kind of get it?
#though easier in bayesian land to get the CIs
# on the 50 Surv date which we don't have
# here we only have point est?
ggplot (tst,
        aes(x=yday, y=surv, group = subject, color = factor(year)))+
  geom_line()+
  geom_vline(aes(xintercept = doy_event, color = factor(year)))+
  facet_wrap(~sensornode)

OneMinusPredHaz <- 1 - predict(Gompertz_Model_Baseline5_eos_wtemp_3, 
                               newdata = phen_data_eos %>%
                                 ungroup %>%
                                 filter(days_since_snow>0), 
                               type = "response")
#calculate individual predicted survival curves
PredSurv <- aggregate(formula = OneMinusPredHaz ~ subject, 
                      data = phen_data_eos %>%
                        ungroup %>%
                        filter(days_since_snow>0), 
                      FUN=cumprod)


#as you add in pairs, having both GDD and yday in there
# is actually the best
summary (sos_yday)

car::vif(eos_all)

Gompertz_Model_Baseline <- glm(formula = status ~ GDD + days_since_snow,
                               family = binomial(link = "cloglog"),
                               data = phen_data_sos %>%
                                 ungroup %>%
                                 filter(days_since_snow>0))

Gompertz_Model_Baseline_eos <- glm(formula = status ~ GDD + days_since_snow,
                                   family = binomial(link = "cloglog"),
                                   data = phen_data_eos %>%
                                     ungroup %>%
                                     filter(days_since_snow>0))
summary (Gompertz_Model_Baseline_eos) #warmer and earlier snowmelt advance fall?


#order agnostic
Gompertz_Model_Baseline2 <- glm(formula = status ~ days_since_snow + GDD,
                                family = binomial(link = "cloglog"),
                                data = phen_data_sos %>%
                                  ungroup %>%
                                  filter(days_since_snow>0))

Gompertz_Model_Baseline3 <- glm(formula = status ~ yday + GDD,
                                family = binomial(link = "cloglog"),
                                data = phen_data_sos %>%
                                  ungroup %>%
                                  filter(days_since_snow>0))

Gompertz_Model_Baseline3_eos <- glm(formula = status ~ yday + GDD,
                                    family = binomial(link = "cloglog"),
                                    data = phen_data_eos %>%
                                      ungroup %>%
                                      filter(days_since_snow>0))

summary (Gompertz_Model_Baseline3_eos) #no, it's just that GDD is counting time
# and hazard increases into fall


#extra info added by yday, partially substitutive
summary (Gompertz_Model_Baseline3)

Gompertz_Model_Baseline4 <- glm(formula = status ~ yday * GDD,
                                family = binomial(link = "cloglog"),
                                data = phen_data_sos %>%
                                  ungroup %>%
                                  filter(days_since_snow>0))
#additive effects, not interactive
summary (Gompertz_Model_Baseline4)

Gompertz_Model_Baseline5<- glm(formula = status ~ yday + GDD + moist_7d,
                               family = binomial(link = "cloglog"),
                               data = phen_data_sos %>%
                                 ungroup %>%
                                 filter(days_since_snow>0))

#dryer it gets the more the hazard increases - this makes sense!
#this one is actually better than the thresholded one
# for whatever reason...
Gompertz_Model_Baseline5_eos<- glm(formula = status ~ yday + moist_7d,
                                   family = binomial(link = "cloglog"),
                                   data = phen_data_eos %>%
                                     ungroup %>%
                                     filter(days_since_snow>0))
summary (Gompertz_Model_Baseline5_eos)

Gompertz_Model_Baseline5_eos_chai<- glm(formula = status ~ yday + chai_moist_7d,
                                        family = binomial(link = "cloglog"),
                                        data = phen_data_eos %>%
                                          ungroup %>%
                                          filter(days_since_snow>0))
summary (Gompertz_Model_Baseline5_eos_chai)

#and even better if you put in temp
# basically the coldness and the drying co-control fall
Gompertz_Model_Baseline5_eos_wtemp<- glm(formula = status ~ yday + moist_7d + temp_7d,
                                         family = binomial(link = "cloglog"),
                                         data = phen_data_eos %>%
                                           ungroup %>%
                                           filter(days_since_snow>0))
summary (Gompertz_Model_Baseline5_eos_wtemp)


Gompertz_Model_Baseline5_eos_wtemp_2<- glm(formula = status ~ days_since_snow + moist_7d + temp_7d,
                                           family = binomial(link = "cloglog"),
                                           data = phen_data_eos %>%
                                             ungroup %>%
                                             filter(days_since_snow>0))
summary (Gompertz_Model_Baseline5_eos_wtemp_2)
summary (Gompertz_Model_Baseline5_eos_wtemp)

#probably cannot put both days since snow and yday in the same model?
Gompertz_Model_Baseline5_eos_wtemp_3<- glm(formula = status ~ days_since_snow + 
                                             yday+ moist_7d + temp_7d,
                                           family = binomial(link = "cloglog"),
                                           data = phen_data_eos %>%
                                             ungroup %>%
                                             filter(days_since_snow>0))
summary (Gompertz_Model_Baseline5_eos_wtemp_3)
#summary (Gompertz_Model_Baseline5_eos_wtemp)

#make the brms version
#probably cannot put both days since snow and yday in the same model?
Gompertz_Model_Baseline5_eos_wtemp_3<- glm(formula = status ~ days_since_snow + 
                                             yday+ moist_7d + temp_7d,
                                           family = binomial(link = "cloglog"),
                                           data = phen_data_eos %>%
                                             ungroup %>%
                                             filter(days_since_snow>0))
summary (Gompertz_Model_Baseline5_eos_wtemp_3)


Gompertz_Model_Baseline5_sos_wtemp_3<- glm(formula = status ~ days_since_snow + 
                                             yday+ moist_7d + temp_7d,
                                           family = binomial(link = "cloglog"),
                                           data = phen_data_sos %>%
                                             ungroup %>%
                                             filter(days_since_snow>0))

plot(effects::allEffects(Gompertz_Model_Baseline5_eos_wtemp_3), type = 'response',
     rug = FALSE,
     axes=list(
       grid=TRUE,
       y=list(rotate=30,
              lim = c(0,0.4)
       )))

#I think this fixes
plot(effects::allEffects(Gompertz_Model_Baseline5_eos_wtemp_3,
                         fixed.predictors = list(given.values=c(yday=200, 
                                                                days_since_snow = 60))),
     type = 'response',
     rug = FALSE,
     axes=list(
       grid=TRUE,
       y=list(rotate=30,
              lim = c(0,0.4)
       )))

Data_DevResid <- tibble(Pred_Haz = predict(Gompertz_Model_Baseline5_eos_wtemp_3, type = "response"),
                        Event = pull(phen_data_eos %>%
                                       ungroup %>%
                                       filter(days_since_snow>0), status),
                        ID = pull(phen_data_eos %>%
                                    ungroup %>%
                                    filter(days_since_snow>0), subject)) %>%
  mutate(DevRes = if_else(Event == 0, 
                          -sqrt(-2*log(1-Pred_Haz)),
                          sqrt(-2*log(Pred_Haz))))
Data_DevResid %>%
  ggplot(aes(x = ID, y = DevRes)) +
  geom_point()

# a few deviance resids >3
# so some poor fit to events, better
# fit to nonevents


#calculate 1 - predicted hazards
OneMinusPredHaz <- 1 - predict(Gompertz_Model_Baseline5_eos_wtemp_3, 
                               newdata = phen_data_eos %>%
                                 ungroup %>%
                                 filter(days_since_snow>0), 
                               type = "response")
#calculate individual predicted survival curves
PredSurv <- aggregate(formula = OneMinusPredHaz ~ subject, 
                      data = phen_data_eos %>%
                        ungroup %>%
                        filter(days_since_snow>0), 
                      FUN=cumprod)

PredErrCurve <- predErrDiscShort(timepoints = 1:100, 
                                 estSurvList = PredSurv[[2]], #survival curves
                                 newTime = phen_data_eos %>%
                                   ungroup %>%
                                   filter(days_since_snow>0)$yday, #time points in the test set
                                 newEvent = phen_data_eos %>%
                                   ungroup %>%
                                   filter(days_since_snow>0)$status, #event in the test set
                                 trainTime = Scania_Person_Train$exit, #time points in the training set
                                 trainEvent = Scania_Person_Train$event) #event in the training set



plot(effects::allEffects(Gompertz_Model_Baseline5_sos_wtemp_3), type = 'response',
     rug = FALSE,
     axes=list(
       grid=TRUE,
       y=list(rotate=30,
              lim = c(0,0.4)
       )))

#estimates mean for every one deg increase in temp
# pro of fall drops ~40%
jtools::summ(Gompertz_Model_Baseline5_eos_wtemp_3, exp = T, digits =8)

#scaled values - the times (calendar day, days_since snow have strongest effects)
# i.e. there is some fixedy in the timing of fall unrelated to the above?
jtools::summ(Gompertz_Model_Baseline5_eos_wtemp_3, exp = T, digits =8, scale = TRUE)

#more days_since_snow increases the hazard, as do later dates
jtools::plot_summs(Gompertz_Model_Baseline5_eos_wtemp_3, exp = T, scale = T)
plot(effects::allEffects(Gompertz_Model_Baseline5_eos_wtemp_3), type = 'response')

#probability of fall increases w days_since_snow, yday, 
#decreases with moisture, and rolling mean temp in the last 5d
plot(effects::allEffects(Gompertz_Model_Baseline5_eos_wtemp_3, 
                         axes=list(y=list(type="response", grid=TRUE,
                                          lab="probabilty scale, probability labels"))))

plot(effects::allEffects(Gompertz_Model_Baseline5_sos_wtemp_3, 
                         axes=list(y=list(type="response", grid=TRUE,
                                          lab="probabilty scale, probability labels"))))

#try the bayes mod
library (brms)


set_inits <- function(seed = 1) {
  set.seed(seed)
  list(
    Intercept = rnorm(n = 1, mean = 0, sd = 0.001),
    b  = rnorm(n = 4, mean = 0, sd = 0.001)
  )
}

set_inits_cheat <- function(seed = 1) {
  set.seed(seed)
  list(
    Intercept = rnorm(n = 1, mean = -13, sd = 1),
    b  = c(
      rnorm(n = 1, mean = 0.02, sd = 0.001),
      rnorm(n = 1, mean = 0.07, sd = 0.001),
      rnorm(n = 1, mean = -12.48, sd = 0.001),
      rnorm(n = 1, mean = -0.44, sd = 0.001)))
}

# try it out
set_inits(seed = 0)

my_second_list_of_inits <- list(
  # different seed values will return different results
  set_inits(seed = 1),
  set_inits(seed = 2)
)

my_second_list_of_inits <- list(
  # different seed values will return different results
  set_inits_cheat(seed = 1),
  set_inits_cheat(seed = 2)
)



#have to set the inits at 0 to make it go
Bayes_Model <- brm(formula = status ~ days_since_snow + 
                     yday+ moist_7d + temp_7d,  
                   data = phen_data_eos %>%
                     ungroup %>%
                     filter(days_since_snow>0), 
                   family = bernoulli(link = "cloglog"),
                   warmup = 500, 
                   iter = 2000, 
                   chains = 2, 
                   inits= "0", 
                   cores=2)#,
#seed = 123)
#sce is totally unclear why the 2nd chain cant get inits
Bayes_Modelwo_snowmelt <- update(Bayes_Model, formula. = ~ . - days_since_snow,
                                 seed = 123)
#seed = 123)




Bayes_Model_logit <- brm(formula = status ~ days_since_snow + 
                           yday+ moist_7d + temp_7d,  
                         data = phen_data_eos %>%
                           ungroup %>%
                           filter(days_since_snow>0), 
                         family = bernoulli(link = "logit"),
                         warmup = 500, 
                         iter = 2000, 
                         chains = 2, 
                         inits= "0", 
                         cores=2)#,
#apparently logit more numerically stable
#https://discourse.mc-stan.org/t/initialization-error-try-specifying-initial-values-reducing-ranges-of-constrained-values-or-reparameterizing-the-model/4401/2
Bayes_Modelwo_snowmelt_logit <- update(Bayes_Model_logit, formula. = ~ . - days_since_snow,
                                       seed = 123)



Bayes_Model_logit_sos <- brm(formula = status ~ days_since_snow + 
                               yday+ moist_7d + temp_7d,  
                             data = phen_data_sos %>%
                               ungroup %>%
                               filter(days_since_snow>0), 
                             family = bernoulli(link = "logit"),
                             warmup = 500, 
                             iter = 2000, 
                             chains = 2, 
                             inits= "0", 
                             cores=2)#,


Bayes_Model_logit_pop <- brm(formula = status ~ days_since_snow + 
                               yday+ moist_7d + temp_7d,  
                             data = phen_data_pop %>%
                               ungroup %>%
                               filter(days_since_snow>0), 
                             family = bernoulli(link = "logit"),
                             warmup = 500, 
                             iter = 2000, 
                             chains = 2, 
                             inits= "0", 
                             cores=2)#,





#seed = 123)

#comparing w loo and AIC
#https://bookdown.org/content/4253/extending-the-discrete-time-hazard-model.html
Bayes_Model  <- add_criterion(Bayes_Model , c("loo", "waic"))
Bayes_Model_logit  <- add_criterion(Bayes_Model_logit , c("loo", "waic"))
Bayes_Modelwo_snowmelt_logit  <- add_criterion(Bayes_Modelwo_snowmelt_logit , c("loo", "waic"))

#basically the same probit or logit
#it's a little worse if you ditch out days since snowmelt but not sure if 3 is really that much?
loo_compare(Bayes_Model, Bayes_Model_logit, Bayes_Modelwo_snowmelt_logit, criterion = "loo") %>% print(simplify = F)

#https://discourse.mc-stan.org/t/interpreting-elpd-diff-loo-package/1628/4
#need to compare these

library (brms)
#if you do random inits it is unhappy

#problem seemed to be that we didn't have things scaled
Bayes_Model3 <- brm(formula = status ~ scale(days_since_snow) + scale(GDD) + 
                      scale(yday)+ scale(moist_7d) + scale(temp_7d) + (1|sensorid),  
                    data = phen_data_sos %>%
                      ungroup %>%
                      filter(days_since_snow>0),
                    family = bernoulli(link = "logit"),
                    warmup = 500, 
                    iter = 2000, 
                    chains = 2, 
                    inits= "random", 
                    cores=2)

#adding in the ranefs here
Bayes_Model3_pop <- brm(formula = status ~ scale(days_since_snow) + scale(GDD) + 
                          scale(yday)+ scale(moist_7d) + scale(temp_7d) + (1|sensorid),  
                        data = phen_data_pop %>%
                          ungroup %>%
                          filter(days_since_snow>0),
                        family = bernoulli(link = "logit"),
                        warmup = 500, 
                        iter = 2000, 
                        chains = 2, 
                        inits= "random", 
                        cores=2)

Bayes_Model3_eos <- brm(formula = status ~ scale(days_since_snow) + scale(GDD) + 
                          scale(yday)+ scale(moist_7d) + scale(temp_7d) + (1|sensorid),  
                        data = phen_data_pop %>%
                          ungroup %>%
                          filter(days_since_snow>0),
                        family = bernoulli(link = "logit"),
                        warmup = 500, 
                        iter = 2000, 
                        chains = 2, 
                        inits= "random", 
                        cores=2)






#try setting the intis to the posterior
# can we get it to start?
# sometimes...
#this is with the b list

#this seems to work just fine without
# the ranefs, we can add them in as nec
Bayes_Model3 <- brm(formula = status ~ days_since_snow + 
                      yday+ moist_7d + temp_7d,  
                    data = phen_data_eos %>%
                      ungroup %>%
                      filter(days_since_snow>0), 
                    family = bernoulli(link = "cloglog"),
                    warmup = 500, 
                    iter = 2000, 
                    chains = 2, 
                    inits= my_second_list_of_inits , 
                    cores=2)

#

tst = stancode(Bayes_Model)
summary(Gompertz_Model_Baseline5_eos_wtemp_3)

#https://bookdown.org/ajkurz/Statistical_Rethinking_recoded/counting-and-classification.html
b10.1 <- add_criterion(Bayes_Model, "waic")

#https://bookdown.org/content/4253/extending-the-discrete-time-hazard-model.html
# Get the plots generated by conditional_effects()
#  p1 is a list containing four elements. Each element is a plot.
p1 = plot(conditional_effects(Bayes_Model_logit), plot = F)

library (patchwork)
patchwork::wrap_plots(p1) + 
  plot_annotation(title="daily 'risk' of fall as a function of envt drivers and time", 
                  theme=theme(plot.title=element_text(hjust=0.5)))

p1 = plot(conditional_effects(Bayes_Model_logit_sos), plot = F)

patchwork::wrap_plots(p1) + 
  plot_annotation(title="daily 'risk' of spring as a function of envt drivers and time", 
                  theme=theme(plot.title=element_text(hjust=0.5)))

p1 = plot(conditional_effects(Bayes_Model_logit_pop), plot = F)

patchwork::wrap_plots(p1) + 
  plot_annotation(title="daily 'risk' of peak season as a function of envt drivers and time", 
                  theme=theme(plot.title=element_text(hjust=0.5)))


plot(conditional_effects(Bayes_Model_logit), plot = F)[[1]] + ggtitle("daily 'risk' of fall as a function of envt drivers and time")

plot(conditional_effects(Bayes_Model_logit)$`days_since_snow`)

plot(conditional_effects(Bayes_Model_logit), ask = FALSE, ncol =2)

#can't do it via brms bc we dropped subject already
#kfold1 <- kfold(Bayes_Model, chains = 1, group=subject)
subjs=phen_data_eos %>%
  ungroup %>%
  filter(days_since_snow>0)%>%dplyr::pull(subject)
folds1 <- loo::kfold_split_stratified(K = length(unique(subjs)), x = subjs)

kf_rtm0 <- brms::kfold(Bayes_Model, folds = folds1)

kf_rtm0 <- brms::kfold(Bayes_Model3, folds = folds1)

#can just do the regular loo?
#but that is not right bc there arent 3792 values in the loglike
uvsdt_loo <- loo(Bayes_Model,cores = 4)

kfold1 <- brms::kfold(Bayes_Model, chains = 1, reloo = TRUE) #doesn't know what reloo is
kfold1 <- brms::kfold(Bayes_Model3, chains = 1, seed=12345)
kfold1 <- brms::kfold(Bayes_Model3, chains = 1)

model_1_split_1 <- kfold(Bayes_Model, K=10)

kf<-brms::kfold(Bayes_Model,K=5,save_fits=TRUE,compare=TRUE)

#basically gives the same estimates as the glm
#see if we can loo cv it

#first make predictions
preds = predict(Bayes_Model)



formula = status ~ days_since_snow + 
  yday+ moist_7d + temp_7d,
family = binomial(link = "cloglog"),
data = phen_data_sos %>%
  ungroup %>%
  filter(days_since_snow>0))


#make predictions from the model
# or try it the other way
make_preds <- function(newdata = phen_data_eos_pred, model = Gompertz_Model_Baseline5_eos_wtemp){
  PredSurv = data.frame(
    cbind(newdata %>%
            ungroup %>%
            filter(days_since_snow>0),
          hazard = 
            predict(model, 
                    newdata = newdata %>%
                      ungroup %>%
                      filter(days_since_snow>0), 
                    type = "response"))) %>%
    mutate (surv = 1-hazard) %>%
    arrange(subject, days_since_snow) %>%
    group_by(subject) %>%
    mutate(failure = cumprod(surv)) %>%
    mutate(p_occur = hazard * lag(failure, 1))
  pred_trans = PredSurv %>%
    group_by(subject) %>%
    arrange(yday) %>%
    filter(failure<0.5) %>%
    slice(1) %>%
    select(subject, yday)
  pred_v_act = (newdata %>%
                  #filter(subject == "sn_1_2019") %>%
                  select(subject, doy_event, year, sensornode) %>%distinct()) %>%
    left_join(pred_trans %>%
                rename(pred_doy_event = yday))
  rmse = Metrics::rmse(pred_v_act$doy_event, pred_v_act$pred_doy_event)
  naive_model = pred_v_act %>%
    group_by(sensornode) %>%
    summarize(doy_event_mean = mean(doy_event))
  pred_v_act = pred_v_act %>%
    left_join(naive_model)
  rmse_naive = Metrics::rmse(pred_v_act$doy_event, pred_v_act$doy_event_mean)
  
  return(list(PredSurv=PredSurv, pred_trans = pred_trans, pred_v_act = pred_v_act,
              rmse = rmse,
              rmse_naive = rmse_naive))
}

#down to 9d rmse for fall, which isn't *that* bad all things considering
#but then note that the SD (naive model) is only 12d so perhaps this
# is not an amazing accomplishment
preds = make_preds(newdata = phen_data_eos_pred, model = Gompertz_Model_Baseline5_eos_wtemp)

#can drop nearly another day of rmse by
# adding in the snowmelt date
# which seems like it should be totally colinear
# w the yday but somehow not so bad that it won't fit?
preds2 = make_preds(newdata = phen_data_eos_pred, model = Gompertz_Model_Baseline5_eos_wtemp_3)

preds3 = make_preds(newdata = phen_data_eos, model = Gompertz_Model_Baseline5_eos_wtemp_3)

head (preds3$PredSurv)

tst = preds3$PredSurv
tst = tst %>%
  rowwise() %>%
  mutate(ll = dbinom(status, size =1, prob = hazard, log = TRUE ))


#sum (tst$ll)
#[1] -162.9099

# ok have calc-ed the log like correctly
# so could probably figure out how to get this into the loo_compare
# if we go brms route
logLik(Gompertz_Model_Baseline5_eos_wtemp_3)

loglike <- dbinom(X, size = n, prob = p, log = TRUE)

pluginlogscore_3 <- sum(dnorm(kidiq$kid_score, fitted(fit_3), sigma(fit_3), log = TRUE))
round(pluginlogscore_3, 1)


ggplot (preds$PredSurv%>%
          filter(yday<300),
        #filter(subject == 'sn_1_2019'),
        aes(x=yday, y=p_occur))+
  geom_point()+
  geom_vline(data = phen_data_eos_pred %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, doy_event) %>%distinct(),
             aes(xintercept = doy_event))+
  geom_vline(data = preds$pred_trans %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, yday) %>%distinct(),
             aes(xintercept = yday), color = 'red') +
  facet_wrap(~subject)

ggplot (preds$PredSurv%>%
          filter(yday<300),
        #filter(subject == 'sn_1_2019'),
        aes(x=yday, y=p_occur))+
  geom_point()+
  geom_vline(data = phen_data_eos_pred %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, doy_event) %>%distinct(),
             aes(xintercept = doy_event))+
  geom_vline(data = preds$pred_trans %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, yday) %>%distinct(),
             aes(xintercept = yday, color = 'red')) +
  facet_wrap(~subject)

ggplot (preds$pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(year)))+
  geom_point()+
  geom_abline(slope = 1, intercept =0)

ggplot (preds$pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(year)))+
  geom_point()+
  geom_abline(slope = 1, intercept =0)


#model of the mean, spatial rank is kind of preserved
# but the year preds are way off
ggplot (preds2$pred_v_act, 
        aes(x = doy_event, y=doy_event_mean, color = factor(year)))+
  geom_point()+
  geom_abline(slope = 1, intercept =0)


# this partially fixes it, i.e. there are points both above and below
# the abline for each yr
ggplot (preds2$pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(year)))+
  geom_point()+
  geom_abline(slope = 1, intercept =0)





#naive model - take the multiyear average


#slightly underpredicting the tails, i.e. some
# events as early as 215 but model doesnt think they'd happen until 220/230
ggplot (preds$pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(sensornode)))+
  geom_point()+
  geom_abline(intercept =0, slope =1)

#plot the hazard as a function of temp, moist
# or just interpret directly
# or just give an example of making it 5% dryer

preds2_dry = make_preds(newdata = phen_data_eos_pred %>%
                          mutate(moist_7d = moist_7d -0.05),
                        model = Gompertz_Model_Baseline5_eos_wtemp_3)

preds2_warm = make_preds(newdata = phen_data_eos_pred %>%
                           mutate(temp_7d = temp_7d + 1),
                         model = Gompertz_Model_Baseline5_eos_wtemp_3)

#5% dryer leads to on average a 4.5 d earlier drydown
tst = preds2_dry$pred_trans %>%
  rename(pred_yday_dry = yday)%>%
  left_join(preds2$pred_trans)

tst = preds2_warm$pred_trans %>%
  rename(pred_yday_warm = yday)%>%
  left_join(preds2$pred_trans)

#warming by 1d will make fall 3.5 d later.
mean (tst$pred_yday_warm -tst$yday)


%>%
  ggplot(aes(gdpPercap,lifeExp, color=year)) +
  geom_point(aes(fill=year),size=3) +
  scale_x_log10()+
  geom_line(aes(group = paired),
            color="grey",
            arrow = arrow(type = "closed",
                          length=unit(0.075, "inches")))
ggsave("customizing_scatterplot_connecting_paired_points_with_lines_arrows_ggplot2.png")






#think interestingly about how we can get this into a space-time-subst framework?        

ggplot(make_preds(newdata = phen_data_eos_pred, model = Gompertz_Model_Baseline5_eos_wtemp))

pred_trans = PredSurv_pop %>%
  group_by(subject) %>%
  arrange(yday) %>%
  filter(failure<0.5) %>%
  slice(1) %>%
  select(subject, yday)

#looks ok ish maybe?
# 6 always occurs way earlier than expected,
# rest not too TOO bad?
ggplot (PredSurv_pop %>%
          filter(yday<300),
        #filter(subject == 'sn_1_2019'),
        aes(x=yday, y=p_occur))+
  geom_point()+
  geom_vline(data = phen_data_pop_pred %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, doy_event) %>%distinct(),
             aes(xintercept = doy_event))+
  geom_vline(data = pred_trans %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, yday) %>%distinct(),
             aes(xintercept = yday, color = 'red')) +
  facet_wrap(~subject)

pred_v_act = (phen_data_pop_pred %>%
                #filter(subject == "sn_1_2019") %>%
                select(subject, doy_event, year, sensornode) %>%distinct()) %>%
  left_join(pred_trans %>%
              rename(pred_doy_event = yday))

#this basically has the yrs right - no sig bias
ggplot (pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(year)))+
  geom_point()

#sn 6 is definitely biased late
ggplot (pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(sensornode)))+
  geom_point()





##end functiondef

#here, increasing moisture actually
# decreases the hazard (but that's probably bc
# we are counting from the beginning of the year
# where the moisture is high and it's not fall?)
summary(Gompertz_Model_Baseline5_eos)


#dryer it gets the more the hazard increases - this makes sense!
#also colder falls increase the hazard
# also 
Gompertz_Model_Baseline5_eos<- glm(formula = status ~ yday + moist_7d +temp_7d,
                                   family = binomial(link = "cloglog"),
                                   data = phen_data_eos %>%
                                     ungroup %>%
                                     filter(days_since_snow>0))
#here, increasing moisture actually
# decreases the hazard (but that's probably bc
# we are counting from the beginning of the year
# where the moisture is high and it's not fall?)
summary(Gompertz_Model_Baseline5_eos)

# or try it the other way
PredSurv_eos = data.frame(
  cbind(phen_data_eos_pred %>%
          ungroup %>%
          filter(days_since_snow>0),
        hazard = 
          predict(Gompertz_Model_Baseline5_eos, 
                  newdata = phen_data_eos_pred %>%
                    ungroup %>%
                    filter(days_since_snow>0), 
                  type = "response"))) %>%
  mutate (surv = 1-hazard) %>%
  arrange(subject, days_since_snow) %>%
  group_by(subject) %>%
  mutate(failure = cumprod(surv)) %>%
  mutate(p_occur = hazard * lag(failure, 1)) 

pred_trans = PredSurv_pop %>%
  group_by(subject) %>%
  arrange(yday) %>%
  filter(failure<0.5) %>%
  slice(1) %>%
  select(subject, yday)

#looks ok ish maybe?
# 6 always occurs way earlier than expected,
# rest not too TOO bad?
ggplot (PredSurv_pop %>%
          filter(yday<300),
        #filter(subject == 'sn_1_2019'),
        aes(x=yday, y=p_occur))+
  geom_point()+
  geom_vline(data = phen_data_pop_pred %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, doy_event) %>%distinct(),
             aes(xintercept = doy_event))+
  geom_vline(data = pred_trans %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, yday) %>%distinct(),
             aes(xintercept = yday, color = 'red')) +
  facet_wrap(~subject)

pred_v_act = (phen_data_pop_pred %>%
                #filter(subject == "sn_1_2019") %>%
                select(subject, doy_event, year, sensornode) %>%distinct()) %>%
  left_join(pred_trans %>%
              rename(pred_doy_event = yday))

#this basically has the yrs right - no sig bias
ggplot (pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(year)))+
  geom_point()

#sn 6 is definitely biased late
ggplot (pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(sensornode)))+
  geom_point()





Gompertz_Model_Baseline5_eos_2<- glm(formula = status ~ yday * moist_7d,
                                     family = binomial(link = "cloglog"),
                                     data = phen_data_eos %>%
                                       ungroup %>%
                                       filter(days_since_snow>0))
summary(Gompertz_Model_Baseline5_eos_2)



Gompertz_Model_Baseline5_eos_4<- glm(formula = status ~ GDD + moist_7d,
                                     family = binomial(link = "cloglog"),
                                     data = phen_data_eos %>%
                                       ungroup %>%
                                       filter(yday>200))
summary(Gompertz_Model_Baseline5_eos_4)
anova(Gompertz_Model_Baseline5_eos_3, Gompertz_Model_Baseline5_eos_4)




moist_pop<- glm(formula = status ~ yday + GDD + moist_7d,
                family = binomial(link = "cloglog"),
                data = phen_data_pop %>%
                  ungroup %>%
                  filter(days_since_snow>0))

temp_pop<- glm(formula = status ~ yday + temp_7d + days_since_snow,
               family = binomial(link = "cloglog"),
               data = phen_data_pop %>%
                 ungroup %>%
                 filter(days_since_snow>0))

temp_pop2<- glm(formula = status ~ yday + GDD + days_since_snow,
                family = binomial(link = "cloglog"),
                data = phen_data_pop %>%
                  ungroup %>%
                  filter(days_since_snow>0))
summary (temp_pop2)

# Gompertz_Model_Baseline2_pop <- glm(formula = status ~ days_since_snow + GDD,
#                                     family = binomial(link = "cloglog"),
#                                     data = phen_data_pop %>%
#                                       ungroup %>%
#                                       filter(days_since_snow>0))


anova(Gompertz_Model_Baseline3, Gompertz_Model_Baseline4, test ="Chisq")

OneMinusPredHaz <- 1 - predict(Gompertz_Model_Baseline3, 
                               newdata = phen_data_sos %>%
                                 ungroup %>%
                                 filter(days_since_snow>0), 
                               type = "response")

#calculate individual predicted survival curves
#sce notes this is pretty close but need to cumprod
# them in a different way and also for dates
# future to see how close it is
PredSurv <- aggregate(formula = OneMinusPredHaz ~ subject, 
                      data = phen_data_sos %>%
                        ungroup %>%
                        filter(days_since_snow>0), 
                      FUN=cumprod)

# or try it the other way
PredSurv = data.frame(
  cbind(phen_data_sos_pred %>%
          ungroup %>%
          filter(days_since_snow>0),
        hazard = 
          predict(Gompertz_Model_Baseline3, 
                  newdata = phen_data_sos_pred %>%
                    ungroup %>%
                    filter(days_since_snow>0), 
                  type = "response"))) %>%
  mutate (surv = 1-hazard) %>%
  arrange(subject, days_since_snow) %>%
  group_by(subject) %>%
  mutate(failure = cumprod(surv)) %>%
  mutate(p_occur = hazard * lag(failure, 1))



# or try it the other way
PredSurv_pop = data.frame(
  cbind(phen_data_pop_pred %>%
          ungroup %>%
          filter(days_since_snow>0),
        hazard = 
          predict(Gompertz_Model_Baseline3_pop, 
                  newdata = phen_data_pop_pred %>%
                    ungroup %>%
                    filter(days_since_snow>0), 
                  type = "response"))) %>%
  mutate (surv = 1-hazard) %>%
  arrange(subject, days_since_snow) %>%
  group_by(subject) %>%
  mutate(failure = cumprod(surv)) %>%
  mutate(p_occur = hazard * lag(failure, 1)) 

pred_trans = PredSurv_pop %>%
  group_by(subject) %>%
  arrange(yday) %>%
  filter(failure<0.5) %>%
  slice(1) %>%
  select(subject, yday)

#looks ok ish maybe?
# 6 always occurs way earlier than expected,
# rest not too TOO bad?
ggplot (PredSurv_pop %>%
          filter(yday<300),
        #filter(subject == 'sn_1_2019'),
        aes(x=yday, y=p_occur))+
  geom_point()+
  geom_vline(data = phen_data_pop_pred %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, doy_event) %>%distinct(),
             aes(xintercept = doy_event))+
  geom_vline(data = pred_trans %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, yday) %>%distinct(),
             aes(xintercept = yday, color = 'red')) +
  facet_wrap(~subject)

pred_v_act = (phen_data_pop_pred %>%
                #filter(subject == "sn_1_2019") %>%
                select(subject, doy_event, year, sensornode) %>%distinct()) %>%
  left_join(pred_trans %>%
              rename(pred_doy_event = yday))

#this basically has the yrs right - no sig bias
ggplot (pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(year)))+
  geom_point()

#sn 6 is definitely biased late
ggplot (pred_v_act, 
        aes(x = doy_event, y=pred_doy_event, color = factor(sensornode)))+
  geom_point()


#can get peak doy correct to 8 days, which isn't too to bad
Metrics::rmse (pred_v_act$doy_event, pred_v_act$pred_doy_event)



#or ggplot vs 50% predictions
PredSurv_pop %>%
  group_by(subject) %>%
  arrange(yday) 


#looks ok ish maybe?
# 6 always occurs way earlier than expected,
# rest not too TOO bad?
ggplot (PredSurv %>%
          filter(yday<250),
        #filter(subject == 'sn_1_2019'),
        aes(x=yday, y=p_occur))+
  geom_point()+
  geom_vline(data = phen_data_sos_pred %>%
               #filter(subject == "sn_1_2019") %>%
               select(subject, doy_event) %>%distinct(),
             aes(xintercept = doy_event))+
  facet_wrap(~subject)

%>%
  cbine


dif(PredSurv[[2]][1]
    
    
    
    library (jtools) # won't load just yet, filled too much memory need to crash sesson
    
    
    mutate(ct = dplyr::n()) %>%
      pivot_wider(names_from = status, values_from =ct), values_fill =0)



#see if we can get sensible predictions from this?
#sce attempt
preds = phen_data_sos %>%
  filter(days_since_snow>0) %>%
  mutate(log_haz = (1/(1+exp((2.2/39.33) * (158.04-GDD))))) %>%
  mutate(survival = exp(-log_haz))

library (tidybayes) # for fitted draws
#predict posterior means. course notes say automatically adds in the offset when it predicts the posterior. 
#but sce needs to figure out what exactly fitted_draws is doing for nonlinear formulas
PP_factor <- fitted_draws(fit4a4, newdata = phen_data_sos_pred, scale = 'linear') 
PP_factor_2 <- fitted_draws(fit4a4, newdata = phen_data_sos_pred, scale = 'response') 

surv_factor <- PP_factor %>% 
  #filter(.draw == 1, x=='Maintained') %>% 
  select(-.chain, -.iteration) %>% 
  #nope this is not right
  #group_by(.draw, subject, yday) %>% 
  group_by(.draw, subject, yday) %>% 
  mutate(S = exp(-cumsum(exp(.value)))) 

surv_factor_2 <- PP_factor %>% 
  #filter(.draw == 1, x=='Maintained') %>% 
  select(-.chain, -.iteration) %>% 
  #nope this is not right
  #group_by(.draw, subject, yday) %>% 
  group_by(.draw, subject, yday) %>% 
  mutate(S = exp(-cumsum(.value)))


# haz_ordered <- PP_factor %>% 
#   #filter(.draw == 1, x=='Maintained') %>% 
#   select(-.chain, -.iteration) %>% 
#   #nope this is not right
#   #group_by(.draw, subject, yday) %>% 
#   group_by(.draw, subject) %>%
#   arrange(yday) %>% 
#   mutate(S = exp(-cumsum(exp(.value)))) 

#cumulative hazard = -log(S(t))
#or Survival =  exp (-cumsum(lambda(t)))

surv_factor_ordered <- PP_factor %>% 
  #filter(.draw == 1, x=='Maintained') %>% 
  select(-.chain, -.iteration) %>% 
  #nope this is not right
  #group_by(.draw, subject, yday) %>% 
  group_by(.draw, subject) %>%
  arrange(yday) %>% 
  mutate(S = exp(-cumsum(exp(.value)))) 

surv_factor_mean <- surv_factor_ordered %>% 
  group_by(subject, year, yday) %>% 
  summarise(S = mean(S))

haz_factor_mean <- PP_factor %>% 
  #filter(.draw == 1, x=='Maintained') %>% 
  select(-.chain, -.iteration) %>% 
  #nope this is not right
  #group_by(.draw, subject, yday) %>% 
  group_by(subject, year,yday) %>%
  arrange(yday) %>% 
  summarise(h = mean(.value))

#sce see if we can smartly plot the hazards?
#plot hazards
# the haz seems right as a function of GDDs
# but the surv curve indicates the predicted day of transition
# should be way earlier than it is
p_haz_brms_factor <- ggplot() + #filter(surv_factor %>%
  #filter(subject == "sn_1_2019"), .draw %in% 1:1000), aes(yday, S)) +
  #geom_jitter(width = .1, height = 0, alpha = .01) +
  geom_line(data = haz_factor_mean%>%
              filter(subject == "sn_1_2019"), aes(yday, h))+
  geom_vline(data = phen_data_sos_pred %>%
               filter(subject == "sn_1_2019") %>%
               select(subject, doy_event) %>%distinct(),
             aes(xintercept = doy_event))+
  #geom_line(data = surv_obs, aes(tstart, surv, group = x), alpha = .7)+
  labs(x = 'Time', y = 'hazard', title = 'Poisson, time = factor', color = '') +
  theme(legend.position = 'bottom')+
  facet_wrap(~subject)+
  geom_line(data = surv_factor_mean%>%
              filter(subject == "sn_1_2019"), aes(yday, S), color = "blue")



#plot survival curves
p_surv_brms_factor <- ggplot(filter(surv_factor %>%
                                      filter(subject == "sn_1_2019"), .draw %in% 1:1000), aes(yday, S)) +
  geom_jitter(width = .1, height = 0, alpha = .01) +
  geom_line(data = surv_factor_mean%>%
              filter(subject == "sn_1_2019"), aes(yday, S))+
  #geom_line(data = surv_obs, aes(tstart, surv, group = x), alpha = .7)+
  labs(x = 'Time', y = 'Probability of survival', title = 'Poisson, time = factor', color = '') +
  theme(legend.position = 'bottom')+
  facet_wrap(~subject)

p_surv_brms_factor <- ggplot(filter(surv_factor %>%
                                      filter(subject == "sn_1_2019"), .draw %in% 1:1000), aes(yday, S)) +
  geom_jitter(width = .1, height = 0, alpha = .01) +
  geom_line(data = surv_factor_mean%>%
              filter(subject == "sn_1_2019"), aes(yday, S))+
  #geom_line(data = surv_obs, aes(tstart, surv, group = x), alpha = .7)+
  labs(x = 'Time', y = 'Probability of survival', title = 'Poisson, time = factor', color = '') +
  theme(legend.position = 'bottom')+
  facet_wrap(~subject)




fit4a4_snowmelt <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow)))),
                          p2 + p3  ~ 1, nl = TRUE),
                       data = phen_data_sos %>%
                         filter(days_since_snow>0), family = c("poisson", "identity"),
                       prior = prior4_4,
                       cores =2, control = list(adapt_delta = 0.999))


fit4a4_snowmelt_GDD <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp,
                              p2 + p3  ~ 1, nl = TRUE),
                           data = phen_data_sos %>%
                             filter(days_since_snow>0), family = c("poisson", "identity"),
                           prior = prior4_4,
                           cores =2, control = list(adapt_delta = 0.999))

#estimates a surv curve that has most change around 157 GDD

#this fits w somewhat reasonable values?

forcing = function(p1, p2, p3, x){
  return(p1/(1+exp((2.2/p3) * (p2-x))))
}

#p1 - max
# p2 = point of inflection
# p3 = time between 0.1 and 0.5 and beween 0.5 and 0.9

tst = phen_data_sos %>%
  mutate(pred = forcing(GDD, p1 =1, p2 =157.97, p3 =39.31))

#this plots GDD accumulation until event
#or hazard accumulation as a function of doy
ggplot (tst)+
  geom_line(aes(x=GDD,y=pred, group=subject), color = 'black')+
  geom_line(aes(x=yday, y=pred, group = subject), color = 'red')+
  facet_wrap(~subject)

#different forcings by snowmelt doy
#we are ending up with something where there is none
# this makes no sense
ggplot (tst, aes (x=x, y=pred))+
  geom_line()+
  geom_line(aes(y=pred2), color = 'red')+
  geom_line(aes(y=pred3), color = 'blue')






prior4_5 <- 
  prior(uniform(0, 1), lb = 0, ub =1, nlpar = "p1") + 
  prior(uniform(0, 300), lb = 0, ub =300, nlpar = "p2") + 
  prior(uniform(0, 300), lb = 0, ub =300, nlpar = "p3")


#try to see if ading back in the P1 param helps? 
fit4_5 <- brm(bf(status ~ (p1/(1+exp((2.2/p3) * (p2-GDD)))),
                 p1+p2 + p3  ~ 1, nl = TRUE),
              data = phen_data_sos %>%
                filter(days_since_snow>0), family = "poisson", prior = prior4_5,
              cores =2, control = list(adapt_delta = 0.999))





prior_binom.v3=prior(uniform(-10, 3), lb=-10, ub=3, nlpar = "aSite") +
  #prior is that background recruitment when you scarify also pretty low
  prior(uniform(-10, 3), lb=-10, ub=3, nlpar = "bScar")+
  prior(normal(0, 10), nlpar = "bSeed")+
  prior(normal(0, 5), nlpar = "bSeedScar")

#e.g. something like this
fitbinom_hierarch_odds_allzones.germ.y1.v5 <- brm(bf(germ.y1~ #mu, 
                                                       #mu = 
                                                       #seedT*(exp(b1+log(numSeeded)))+exp(b2),
                                                       #no seeding or scar
                                                       exp(aSite)*plot_size*(seedT-1)*(scarT-1)+
                                                       #scarification only
                                                       exp(bScar)*plot_size*scarT*(-1*(seedT-1))+
                                                       #seedT only
                                                       #(numSeeded*(exp(bSeed)/(1+exp(bSeed)))+exp(aSite))*seedT*(-1*(scarT-1))+
                                                       (viable.seeds.added*(exp(bSeed)/(1+exp(bSeed)))+exp(aSite)*plot_size)*seedT*(-1*(scarT-1))+
                                                       
                                                       #seedT and scarT
                                                       #(numSeeded*(exp(bSeedScar)/(1+exp(bSeedScar)))+exp(aSite+bScar))*seedT*scarT,
                                                       (viable.seeds.added*(exp(bSeed+bSeedScar)/(1+exp(bSeed+bSeedScar)))+exp(bScar)*plot_size)*seedT*scarT,
                                                     aSite ~ 0 + zone + (zone||site),
                                                     bScar ~ 0 + zone + (zone||site),
                                                     bSeed ~ 0 + zone + (zone||site),
                                                     bSeedScar ~ 0 + zone + (zone||site), nl=TRUE), family=c("poisson", "identity"),
                                                  data=df, iter=2000, warmup=1000, prior = prior_binom.v3,
                                                  control = list(adapt_delta = 0.99, max_treedepth = 15), cores=2)



# and one more without the GDD
#make a simpler one where the max is just 1
prior5 <- 
  prior(normal(3, 1000), nlpar = "p2") + 
  prior(normal(50, 1000), nlpar = "p3") 

fit5 <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow)))),
               p2 + p3  ~ 1, nl = TRUE),
            data = phen_data_sos, family = "poisson", prior = prior5,
            cores =2)

# won't solve, has divergent transitions after warmup but rha - 1.38
fit5a <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow)))),
                p2 + p3  ~ 1, nl = TRUE),
             data = phen_data_sos %>%
               filter(days_since_snow>0), family = "poisson", prior = prior5,
             cores =2)

summary (fit3)
summary (fit3a)
summary (fit4)
summary (fit4a)
summary (fit5)
summary (fit5a)

#sce also might just want to try this as days_since_snow
# seems dumb otherwise

phen_data_pop = phen_data_pop %>%
  unite(subject, sensornode, year, remove = FALSE) %>%
  #note that all the clim should be shifted back one
  # doy for forcings but ok for now
  mutate(status = ifelse(yday == doy_event, 1, 0),
         last_day = yday-1,
         days_since_snow = yday - snowmelt_doy_infilled) %>%
  mutate(days_since_snow = ifelse(days_since_snow<0, 0, days_since_snow),
         sensorid = as.numeric(as.factor(sensornode))) %>%
  arrange(sensornode, year, yday)%>%
  group_by(sensornode, year) %>%
  mutate(DD = ifelse(soiltemp_5cm_avg_ave>0&snow==1, soiltemp_5cm_avg_ave, 0),
         GDD = cumsum(DD)) %>%
  ungroup()


phen_data_eos = phen_data_eos %>%
  unite(subject, sensornode, year, remove = FALSE) %>%
  #note that all the clim should be shifted back one
  # doy for forcings but ok for now
  mutate(status = ifelse(yday == doy_event, 1, 0),
         last_day = yday-1,
         days_since_snow = yday - snowmelt_doy_infilled) %>%
  mutate(days_since_snow = ifelse(days_since_snow<0, 0, days_since_snow),
         sensorid = as.numeric(as.factor(sensornode))) %>%
  arrange(sensornode, year, yday)%>%
  group_by(sensornode, year) %>%
  mutate(DD = ifelse(soiltemp_5cm_avg_ave>0&snow==1, soiltemp_5cm_avg_ave, 0),
         GDD = cumsum(DD)) %>%
  ungroup()


#try to do an AFT on days since snow?
# fit <- coxph(Surv(last_day, yday, status) ~ days_since_snow,
#              data=phen_data %>%
#                filter(days_since_snow>0))

#this is not quite
# I think what I actually want to do
# is make the surv process on days_since_snow directly
#but this does run
mod_tvc_sce <- stan_surv(
  formula = Surv(last_day, yday, status) ~ days_since_snow,
  basehaz = 'exp',
  data = phen_data_sos %>%
    filter(days_since_snow>0),
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)

plot (mod_tvc_sce)

saveRDS(mod_tvc_sce, 'sce-tst/stan_temp/mod_tvc_sce.rds')

#view the weibull aft model by group
#nd <- data.frame(sensornode = unique (phen_data$sensornode))
nd <- data.frame(days_since_snow=seq(min(phen_data_sos$days_since_snow),max(phen_data_sos$days_since_snow)))


#sce need to sort out the warnings here
ps_sos_dss <- posterior_survfit(mod_tvc_sce,
                                newdata = nd,
                                times = 0,
                                type = c('haz'),
                                extrapolate = TRUE,
                                control = list(edist = 1000))

tst <-ps_sos_dss %>%
  as.data.frame(.) %>%
  rename(med = median) %>%
  select(-time) %>%
  distinct()

#this has hazard increasing as a function of days-since-snow
# so that seems rightish?
ggplot (tst, aes(x=id, y=med))+
  geom_line()

#but to convert it into probabilities of transitions
#we'd need to go back to set
# days_since_snow to 0 for all the prior
# events and then predict forward? Would that make sense?

#nevertheless let's see if we can get GDD in there too
#nope, can't put both in the same model
# and get it to converge
# mod_tvc_sce_GDD <- stan_surv(
#   formula = Surv(last_day, yday, status) ~ days_since_snow + GDD,
#   basehaz = 'exp',
#   data = phen_data_sos %>%
#     filter(days_since_snow>0),
#   chains = 2,
#   cores = 2,
#   seed = 1234,
#   iter = 10000)

#model doy as a function of GDD,
# but filtered to days_since_snow>0
# shows a significant effect of GDD on surv
# e.g. higher GDD is later surv - can that possibly
# be right?
# I think (though am not sure)
# that it is on the hazard ratio framework
# so for every one unit increase in GDD you get a 2% increase
# in the instantanous rate of the event?
# or will get there 2% faster
mod_tvc_sce_GDD <- stan_surv(
  formula = Surv(last_day, yday, status) ~ GDD,
  basehaz = 'exp',
  data = phen_data_sos %>%
    filter(days_since_snow>0),
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)

summary (mod_tvc_sce_GDD, digits = 4)
saveRDS(mod_tvc_sce, 'sce-tst/stan_temp/mod_tvc_sce.rds')

#try generating predictions
#but this is weird bc it's
# making all the combinations of times and GDDs
# not time-varying covariate predictions



#ok now have them in discrete time ish
# haxard is the instantaneous rate of occurrence
ps_haz <- posterior_survfit(mod_tvc_sce_GDD,
                            newdata = phen_data_sos %>%
                              filter(days_since_snow>0),
                            type = 'haz',
                            times = 0,
                            extrapolate = TRUE,
                            control = list(epoints = 10, edist = 9))

#cumhaz
# integral of hazard from t=0 -> time t
ps_cumhaz <- posterior_survfit(mod_tvc_sce_GDD,
                               newdata = phen_data_sos %>%
                                 filter(days_since_snow>0),
                               type = 'cumhaz',
                               times = 0,
                               extrapolate = TRUE,
                               control = list(epoints = 10, edist = 9))

surv_function = function(x){
  return (ps_haz$median[1] * x)
}

surv_function_log = function(x){
  return (log(ps_haz$median[1]) * x)
}

surv_function_exp = function(x){
  return (exp(ps_haz$median[1]) * x)
}


integrate(surv_function, 0,1)
integrate(surv_function_exp, 0,1)
integrate(surv_function, 1,2)
exp(-1*(integrate(surv_function, 0,2)$value))

#ok, confirming that exp(H) is surv
# now just need to get from little h to big H
#seems right up to rounding
#now I just need to figure out how to do it with time-varying hazards
calcs = ps_cumhaz %>%
  mutate(sce_surv_est_from_cumhaz = exp(-median)) %>%
  left_join(., ps_surv %>% rename(surv = median) %>%
              select(id, time, surv)) %>%
  left_join(.,
            ps_haz %>%
              mutate(sce_surv_est_from_haz = exp(-1*median*time))%>%
              select(id, time, sce_surv_est_from_haz))

#cumhaze incrementes by .0063 every time step
# the hazard is constant at 0.006079925


#surv = esp (-cumhaz)
ps_surv <- posterior_survfit(mod_tvc_sce_GDD,
                             newdata = phen_data_sos %>%
                               filter(days_since_snow>0),
                             type = 'surv',
                             times = 0,
                             extrapolate = TRUE,
                             control = list(epoints = 10, edist = 9))

tst = ps_cumhaz %>%
  as.data.frame() %>%
  mutate(difhaz = c(NA, diff(median)))

#so need to figure out how to make the function
# that integrates instantanous hazard rate: 0.0061
# to cumhaz ->0.006277116

#this is the instantaneous hazard rate
#-0.006079925

1-exp(-0.006079925)

#big H above should be integral of little H above that

ps2 <- posterior_survfit(mod_tvc_sce_GDD,
                         newdata = phen_data_sos %>%
                           filter(days_since_snow>0),
                         extrapolate = FALSE)

pps <- plot(ps) +
  facet_grid(~ id, labeller = labeller(id = panel_labels))

#simple comparison - did not deal w time varying loo yet
compare_models(loo(mod_tvc_sce),
               loo(mod_tvc_sce_GDD))

#fit with all the days, not just the days past snowmelt
mod_tvc_sce_GDD_2 <- stan_surv(
  formula = Surv(last_day, yday, status) ~ GDD,
  basehaz = 'exp',
  data = phen_data_sos,
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)

mod_tvc_sce_2 <- stan_surv(
  formula = Surv(last_day, yday, status) ~ days_since_snow,
  basehaz = 'exp',
  data = phen_data_sos,
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)

#then this would
loo_compare(loo(mod_tvc_sce_2),
            loo(mod_tvc_sce_GDD_2))

#other things to try when time is there
#ie will this fit?
#no this will not fit
# this is fitting an accelleration based on days_since_snow AND then that differs by sensorid
#apparently takes foreer so letting it run then think
# about what it means
mod_tvc_sce_2_ranef <- stan_surv(
  formula = Surv(last_day, yday, status) ~ days_since_snow + (days_since_snow|sensorid),
  basehaz = 'exp',
  data = phen_data_sos,
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)

#no this will not fit either
mod_tvc_sce_2_GDD_ranef <- stan_surv(
  formula = Surv(last_day, yday, status) ~ GDD + (GDD|sensorid),
  basehaz = 'exp',
  data = phen_data_sos,
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)







#gdd model is better than just the days
# past snowmelt (which we kind of knew/expected?)
loo_compare(loo(mod_tvc_sce),
            loo(mod_tvc_sce_GDD))

#need to LOO over individuals not meas within, see p 31
#sce start HERE, assuming the above fit
R> ids <- dat$id
R> ll <- log_lik(mod_tvc)
R> ll <- apply(ll, 1L, function(row) tapply(row, ids, sum))
R> ll <- t(ll)
R> loo(ll)





#this is the prob of peaking in terms of days
# past snowmelt per each node
ggplot (ps_pop, aes(x=time, y=haz_2))+
  geom_line()+
  facet_wrap(~ id,
             labeller = labeller(id = panel_labels))

#summary (mod_tvc_sce)
#sce is unclear on why so many (97.8pct)
# show up as right censored?
#maybe those just the ones that didn't occur ON
# that day though weirdly we know them
# so this is kind of confusing
# as a metric

mod_tvc_sce_2 <- stan_surv(
  formula = Surv(prev_day, days_since_snow, status) ~ 1,
  basehaz = 'exp',
  data = phen_data %>%
    filter(days_since_snow>0) %>%
    mutate(prev_day = days_since_snow-1),
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)

#accellerated failure time model
mod_tvc_sce_3 <- stan_surv(
  formula = Surv(prev_day, days_since_snow, status) ~ 1,
  basehaz = 'exp-aft',
  data = phen_data %>%
    filter(days_since_snow>0) %>%
    mutate(prev_day = days_since_snow-1),
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)


#accellerated failure time model (exp) w frailty terms by node
mod_tvc_sce_4 <- stan_surv(
  formula = Surv(prev_day, days_since_snow, status) ~ 1+ (1|sensorid),
  basehaz = 'exp-aft',
  data = phen_data %>%
    filter(days_since_snow>0) %>%
    mutate(prev_day = days_since_snow-1),
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)


#bsplines frailty w random effects of node
mod_tvc_sce_5 <- stan_surv(
  formula = Surv(prev_day, days_since_snow, status) ~ 1+ (1|sensorid),
  basehaz = 'ms',
  data = phen_data %>%
    filter(days_since_snow>0) %>%
    mutate(prev_day = days_since_snow-1),
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10000)

#seem all to be the same?
plot(posterior_survfit(mod_tvc_sce_3, newdata = nd, type = "haz"))
plot(posterior_survfit(mod_tvc_sce_4, newdata = nd, type = "haz"))


m_aft_pop <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ sensornode,
                       data = phen_data_pop %>%
                         filter(days_since_snow>0) %>%
                         mutate(prev_day = days_since_snow-1),
                       basehaz = "weibull-aft",
                       chains = 2,
                       cores = 2,
                       seed = 1233,
                       iter = 10000)


m_aft_sos <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ sensornode,
                       data = phen_data_sos %>%
                         filter(days_since_snow>0) %>%
                         mutate(prev_day = days_since_snow-1),
                       basehaz = "weibull-aft",
                       chains = 2,
                       cores = 2,
                       seed = 1233,
                       iter = 10000)

m_aft_eos <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ sensornode,
                       data = phen_data_eos %>%
                         filter(days_since_snow>0) %>%
                         mutate(prev_day = days_since_snow-1),
                       basehaz = "weibull-aft",
                       chains = 2,
                       cores = 2,
                       seed = 1233,
                       iter = 10000)

####try all 3 as a function of GDD
m_aft_pop_GDD <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ GDD + (1|sensornode),
                           data = phen_data_pop %>%
                             filter(days_since_snow>0) %>%
                             mutate(prev_day = days_since_snow-1),
                           basehaz = "weibull-aft",
                           chains = 2,
                           cores = 2,
                           seed = 1233,
                           iter = 10000)

#also does not converge
# 1: There were 3963 divergent transitions after warmup. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
# http://mc-stan.org/misc/warnings.html#bfmi-low 
# 3: Examine the pairs() plot to diagnose sampling problems
# 
# 4: The largest R-hat is 2.14, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#r-hat 
# 5: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#bulk-ess 
# 6: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#tail-ess 
# 7: Markov chains did not converge! Do not analyze results! 




m_aft_pop_GDD_2 <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ GDD + (GDD|sensornode),
                             data = phen_data_pop %>%
                               filter(days_since_snow>0) %>%
                               mutate(prev_day = days_since_snow-1),
                             basehaz = "weibull-aft",
                             chains = 2,
                             cores = 2,
                             seed = 1233,
                             iter = 10000)

#the model above does not converge
# Warning messages:
#   1: There were 487 divergent transitions after warmup. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
# http://mc-stan.org/misc/warnings.html#bfmi-low 
# 3: Examine the pairs() plot to diagnose sampling problems
# 
# 4: The largest R-hat is 2.23, indicating chains have not mixed.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#r-hat 
# 5: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#bulk-ess 
# 6: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#tail-ess 
# 7: Markov chains did not converge! Do not analyze results! 



m_aft_sos_GDD <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ GDD + (1|sensornode),
                           data = phen_data_sos %>%
                             filter(days_since_snow>0) %>%
                             mutate(prev_day = days_since_snow-1),
                           basehaz = "weibull-aft",
                           chains = 2,
                           cores = 2,
                           seed = 1233,
                           iter = 10000)

#get this warning for the model above
# Warning messages:
#   1: There were 1422 divergent transitions after warmup. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: There were 1834 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 15. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems
# 
# 4: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#bulk-ess 
# 5: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
# Running the chains for more iterations may help. See
# http://mc-stan.org/misc/warnings.html#tail-ess 

#killed this one bc I don't think it will converge
# m_aft_sos_GDD_2 <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ GDD + (GDD|sensornode),
#                            data = phen_data_sos %>%
#                              filter(days_since_snow>0) %>%
#                              mutate(prev_day = days_since_snow-1),
#                            basehaz = "weibull-aft",
#                            chains = 2,
#                            cores = 2,
#                            seed = 1233,
#                            iter = 10000)

#trying this one next
#AFT where measuring survival time in terms of days since snow
# with the accelleration based on AGDD since snow
m_aft_sos_GDD_simple <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ GDD,
                                  data = phen_data_sos %>%
                                    filter(days_since_snow>0) %>%
                                    mutate(prev_day = days_since_snow-1),
                                  basehaz = "weibull-aft",
                                  chains = 2,
                                  cores = 2,
                                  seed = 1233,
                                  iter = 10000)

#then try this one
m_aft_eos_GDD_simple <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ GDD,
                                  data = phen_data_eos %>%
                                    filter(days_since_snow>0) %>%
                                    mutate(prev_day = days_since_snow-1),
                                  basehaz = "weibull-aft",
                                  chains = 2,
                                  cores = 2,
                                  seed = 1233,
                                  iter = 10000)


m_aft_pop_GDD_simple <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ GDD,
                                  data = phen_data_pop %>%
                                    filter(days_since_snow>0) %>%
                                    mutate(prev_day = days_since_snow-1),
                                  basehaz = "weibull-aft",
                                  chains = 2,
                                  cores = 2,
                                  seed = 1233,
                                  iter = 10000)



#did not bother trying these, pretty sure they won't fit
# m_aft_eos_GDD <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ GDD + (1|sensornode),
#                        data = phen_data_eos %>%
#                          filter(days_since_snow>0) %>%
#                          mutate(prev_day = days_since_snow-1),
#                        basehaz = "weibull-aft",
#                        chains = 2,
#                        cores = 2,
#                        seed = 1233,
#                        iter = 10000)
# 
# m_aft_eos_GDD_2 <- stan_surv(formula = Surv(prev_day, days_since_snow, status) ~ GDD + (GDD|sensornode),
#                            data = phen_data_eos %>%
#                              filter(days_since_snow>0) %>%
#                              mutate(prev_day = days_since_snow-1),
#                            basehaz = "weibull-aft",
#                            chains = 2,
#                            cores = 2,
#                            seed = 1233,
#                            iter = 10000)




####end sce experiment 7/25/2021


#view the weibull aft model by group
nd <- data.frame(sensornode = unique (phen_data$sensornode))
# ps_pop <- posterior_survfit(m_aft,
#                         newdata = nd,
#                         times = 0,
#                         extrapolate = TRUE,
#                         control = list(edist = 100))
# 
# 
# plot(m_aft)
# plot (ps)

#sce need to sort out the warnings here
ps_pop <- posterior_survfit(m_aft_pop,
                            newdata = nd,
                            times = 0,
                            type = c('cdf'),
                            extrapolate = TRUE,
                            control = list(edist = 1000))

#see if we can get the posteriors of the expected date
#I thought it would be the difference in the cumulative haz
# but clearly that is not right
#instead it's the diffs in the cdf (estimated failure probability - I think?)
# kind of looks ok but I think the ranefs just having 
# different intercepts isn't quite right
# we'd want the whole baseline hazard to vary by subject
ps_pop <-ps_pop %>%
  as.data.frame(.) %>%
  rename(med = median) %>%
  group_by(id) %>%
  mutate(haz_2 = c(NA,diff(med)))

#this is the prob of peaking in terms of days
# past snowmelt per each node
ggplot (ps_pop, aes(x=time, y=haz_2))+
  geom_line()+
  facet_wrap(~ id,
             labeller = labeller(id = panel_labels))


#sce need to sort out the warnings here
ps_sos<- posterior_survfit(m_aft_sos,
                           newdata = nd,
                           times = 0,
                           type = c('cdf'),
                           extrapolate = TRUE,
                           control = list(edist = 1000))

#see if we can get the posteriors of the expected date
#I thought it would be the difference in the cumulative haz
# but clearly that is not right
#instead it's the diffs in the cdf (estimated failure probability - I think?)
# kind of looks ok but I think the ranefs just having 
# different intercepts isn't quite right
# we'd want the whole baseline hazard to vary by subject
ps_sos <-ps_sos %>%
  as.data.frame(.) %>%
  rename(med = median) %>%
  group_by(id) %>%
  mutate(haz_2 = c(NA,diff(med)))

#sce need to sort out the warnings here
ps_eos<- posterior_survfit(m_aft_eos,
                           newdata = nd,
                           times = 0,
                           type = c('cdf'),
                           extrapolate = TRUE,
                           control = list(edist = 1000))

#at some point we'll have to unbake the code bc we do
# have an estimate of the baseline survival
# after the last point (assuming exp or weibull, something parametric)

#see if we can get the posteriors of the expected date
#I thought it would be the difference in the cumulative haz
# but clearly that is not right
#instead it's the diffs in the cdf (estimated failure probability - I think?)
# kind of looks ok but I think the ranefs just having 
# different intercepts isn't quite right
# we'd want the whole baseline hazard to vary by subject
ps_eos <-ps_eos %>%
  as.data.frame(.) %>%
  rename(med = median) %>%
  group_by(id) %>%
  mutate(haz_2 = c(NA,diff(med)))

#this is the prob of peaking in terms of days
# past snowmelt per each node
ggplot (ps_sos, aes(x=time, y=haz_2))+
  geom_line(color = "green")+
  geom_line(data = ps_pop, color = "blue")+
  geom_line(data = ps_eos, color = "brown")+
  facet_wrap(~ id,
             labeller = labeller(id = panel_labels))

#other factors we could add in - GDD?
# also need to split these somehow
# sensibly into test-train or something for model eval

#could also wrap it into calendar time?
# doing that comparison you will almost certainly see a wide discrep
# bc of the differences among yrs

#Prediction data for 'posterior_survfit' cannot include delayed entry.
tst = posterior_survfit(mod_tvc_sce_5) 

#this is just the basic plot
#showing the increase in hazard rate as days since snowmelt gets onger
plot (mod_tvc_sce_5)
#I feel like 4 should show aft but it doesn't seem to?
plot (mod_tvc_sce_4)

#nd <- data.frame(sensornode = as.factor(unique (phen_data$sensornode)))
nd <- data.frame(sensorid = unique (phen_data$sensorid))

#sce need to sort out the warnings here
ps <- posterior_survfit(mod_tvc_sce_5,
                        newdata = nd,
                        times = 0,
                        type = c('cumhaz'),
                        extrapolate = TRUE,
                        control = list(edist = 1000))

#see if we can get the posteriors of the expected date
#I thought it would be the difference in the cumulative haz
# but clearly that is not right
ps <-ps %>%
  as.data.frame(.) %>%
  rename(med = median) %>%
  group_by(id) %>%
  mutate(haz_2 = c(NA,diff(med)))

ggplot (ps, aes(x=time, y=haz_2))+
  geom_line()+
  facet_wrap(~id)

#sce need to sort out the warnings here
ps <- posterior_survfit(mod_tvc_sce_5,
                        newdata = nd,
                        times = 0,
                        type = c('cdf'),
                        extrapolate = TRUE,
                        control = list(edist = 1000))

#see if we can get the posteriors of the expected date
#I thought it would be the difference in the cumulative haz
# but clearly that is not right
#instead it's the diffs in the cdf (estimated failure probability - I think?)
# kind of looks ok but I think the ranefs just having 
# different intercepts isn't quite right
# we'd want the whole baseline hazard to vary by subject
ps <-ps %>%
  as.data.frame(.) %>%
  rename(med = median) %>%
  group_by(id) %>%
  mutate(haz_2 = c(NA,diff(med)))

ggplot (ps, aes(x=time, y=haz_2))+
  geom_line()+
  facet_wrap(~ id,
             labeller = labeller(id = panel_labels))




plot(posterior_survfit(mod_tvc_sce_5,
                       newdata = nd,
                       times = 0,
                       type = 'haz',
                       extrapolate = TRUE,
                       control = list(edist = 100)))+
  facet_wrap(~ id,
             labeller = labeller(id = panel_labels))

plot(posterior_survfit(mod_tvc_sce_5,
                       newdata = nd,
                       times = 0,
                       type = 'cumhaz',
                       extrapolate = TRUE,
                       control = list(edist = 100)))+
  facet_wrap(~ id,
             labeller = labeller(id = panel_labels))

plot(posterior_survfit(mod_tvc_sce_5,
                       newdata = nd,
                       times = 0,
                       type = 'cdf',
                       extrapolate = TRUE,
                       control = list(edist = 100)))+
  facet_wrap(~ id,
             labeller = labeller(id = panel_labels))



phen_data %>%
  select(sensornode, sensorid)%>%
  distinct()
panel_labels <- c('1' = "sn_01",
                  '2' = "sn_10",
                  '3' = "sn_11",
                  '4' = "sn_12",
                  '5' = "sn_13",
                  '6' = "sn_14",
                  '7' = "sn_15",
                  '8' = "sn_16",
                  '9' = "sn_17",
                  '10' = "sn_19",
                  '11' = "sn_20",
                  '12' = "sn_21",
                  '13' = "sn_6",
                  '14' = "sn_7",
                  '15' = "sn_8",
                  '16' = "sn_9")

#sce needs to deal with the levels which
# arent getting backed out..
#this is normalized over all of them
plot(posterior_survfit(mod_tvc_sce_5,
                       newdata = nd,
                       times = 0,
                       type = 'cdf',
                       standardise = TRUE,
                       extrapolate = TRUE,
                       control = list(edist = 100)))

#ph <- posterior_survfit(mod1, newdata = nd, type = "haz")


tst <- phen_data %>% select(year, sensornode, snowmelt_doy_infilled, yday, doy_event, days_since_snow) %>%
  dplyr::filter(yday == doy_event) 
#my test
# maybe I want it to just have the
# last dates for the surv?
# but then I don't know how to deal with the non-events
#earlier properly?
tst_2 = phen_data %>%
  filter(days_since_snow>0) %>%
  mutate(prev_day = days_since_snow-1)


#brms diez hazard farm - can't fit
#install.packages("rstanarm", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library (rstanarm)
library (survival)
dat <- survival::pbc
dat <- dat[dat$id <= 312, ]
ts_dat <- survival::pbcseq

#status = 0 = censored, 1 = transplant, 2 = dead
dat2 <- tmerge(dat, dat, id = id,
               death = event(time, as.numeric(status == 2)))

#sce understands this but thinks 
# that we need to expand so we have
# all the doys where it's not yet there
# yet in the model which could make it go slow
dat3 <- tmerge(dat2, ts_dat, id = id,
               ascites = tdc(day, ascites),
               bili = tdc(day, bili),
               albumin = tdc(day, albumin),
               protime = tdc(day, protime),
               alk.phos = tdc(day, alk.phos))

#example, works
mod_tvc <- stan_surv(
  formula = Surv(tstart, tstop, death) ~ log(bili) + log(protime),
  data = dat3,
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10)

#sce shift to constant baseline hazard
mod_tvc <- stan_surv(
  formula = Surv(tstart, tstop, death) ~ log(bili) + log(protime),
  basehaz = 'exp',
  data = dat3,
  chains = 2,
  cores = 2,
  seed = 1234,
  iter = 10)

#or could also do aft one of these two ways
#"exp-aft": an exponential distribution for the event times.
#"weibull-aft": a Weibull distribution for the event times.

#alternatively, to get the days-past-snowmelt in there
# we could to b-splines with
# time functioned as days-past-snowmelt and look for
# sharp peaks in the function [e.g. p 23 example bottom 3 panels]
# but then add in how the baseline hazard is modified 
# by GDD, etc.
# basically I think I want the baseline-hazard to be
# subject (node-specific) effect but need to sort
# how the random effects do that


#this seems to be running and if you google for random effects 
# apparently you can do this but unsure if its
# just on the frailty term of the exp/weibull or could 
# be put into something thresholded

#copied from cox_ph attempt
#library (tidyverse)
#ok try this one with the jeffdiez poisson model
# but fit in brms to make faster

#https://cran.r-project.org/web/packages/brms/vignettes/brms_nonlinear.html

# prior1 <- prior(normal(1, 2), nlpar = "b1") +
#   prior(normal(0, 2), nlpar = "b2")
# fit1 <- brm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE),
#             data = dat1, prior = prior1)

library (brms)

#this will run, not super exciting
fit1 <- brm(bf(status ~ days_since_snow),
            data = phen_data_sos, family = "poisson")

#p1 - max
# p2 = point of inflection
# p3 = time between 0.1 and 0.5 and beween 0.5 and 0.9

#sce could consider fixing p1 to 1 if can't get identifiability here


#just try to see if stan will run it just
# having the snowfree function in there
# I think this might actually need an identity link
# not the one I have but whatevs
prior1 <- prior(normal(0, 1000), nlpar = "p1") +
  prior(normal(3, 1000), nlpar = "p2") + 
  prior(normal(50, 1000), nlpar = "p3")

fit2 <- brm(bf(status ~ p1/(1+exp((2.2/p3) * (p2-days_since_snow))),
               p1 + p2 + p3 ~ 1, nl = TRUE),
            data = phen_data_sos, family = "poisson", prior = prior1)

#see if we can just
# fit the generic model of
# fitting the days-since-snow model
# and see if we can plot the expected values

#make uniform 0-100
prior2 <- prior(uniform(0, 100), nlpar = "p1") +
  prior(uniform(0, 100), nlpar = "p2") + 
  prior(uniform(0, 100), nlpar = "p3")

fit2a <- brm(bf(status ~ p1/(1+exp((2.2/p3) * (p2-days_since_snow))),
                p1 + p2 + p3 ~ 1, nl = TRUE),
             data = phen_data_sos %>%
               filter(days_since_snow>0), family = "poisson", prior = prior1)

fit2a2 <- brm(bf(status ~ p1/(1+exp((2.2/p3) * (p2-days_since_snow))),
                 p1 + p2 + p3 ~ 1, nl = TRUE),
              data = phen_data_sos %>%
                filter(days_since_snow>0), family = "poisson", prior = prior2)


summary (fit2a2)

forcing = function(p1, p2, p3, x){
  return(p1/(1+exp((2.2/p3) * (p2-x))))
}

#p1 - max
# p2 = point of inflection
# p3 = time between 0.1 and 0.5 and beween 0.5 and 0.9

tst = data.frame(x=seq(1:100)) %>%
  mutate(pred = forcing(x, p1 =1, p2 =30, p3 =20))%>%
  mutate(pred2 = forcing(x, p1 =1, p2 =30, p3 =40))%>%
  mutate(pred3 = forcing(x, p1 =0, p2 = 0.86, p3 =1.09))

#different forcings by snowmelt doy
#we are ending up with something where there is none
# this makes no sense
ggplot (tst, aes (x=x, y=pred))+
  geom_line()+
  geom_line(aes(y=pred2), color = 'red')+
  geom_line(aes(y=pred3), color = 'blue')


summary (fit2) # totally doesn't converge
#then the 3rd one would be modifying
# this

prior3 <- prior(normal(0, 1000), nlpar = "p1") +
  prior(normal(3, 1000), nlpar = "p2") + 
  prior(normal(50, 1000), nlpar = "p3") + 
  prior(normal(0, 1000), nlpar = "beta1") 



#beta0[pp]~ dnorm(0,.001) 

#sce should play more with this and then also
# work on integrating the survival function?
# using the code from the peeps
fit3 <- brm(bf(status ~ (p1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
               p1 + p2 + p3 + beta1 ~ 1, nl = TRUE),
            data = phen_data_sos, family = "poisson", prior = prior3,
            cores =2)

fit3a <- brm(bf(status ~ (p1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
                p1 + p2 + p3 + beta1 ~ 1, nl = TRUE),
             data = phen_data_sos %>%
               filter(days_since_snow>0), family = "poisson", prior = prior3,
             cores =2)


#make a simpler one where the max is just 1
prior4 <- 
  prior(normal(3, 1000), nlpar = "p2") + 
  prior(normal(50, 1000), nlpar = "p3") + 
  prior(normal(0, 1000), nlpar = "beta1") 

prior4a <- 
  prior(uniform(0, 50), lb = 0, ub =50, nlpar = "p2") + 
  prior(uniform(0, 50), lb = 0, ub =50,nlpar = "p3") + 
  prior(uniform(-2, 2), lb = -2, ub =2, nlpar = "beta1") 


#no this does not fit
fit4 <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
               p2 + p3 + beta1 ~ 1, nl = TRUE),
            data = phen_data_sos, family = "poisson", prior = prior4,
            cores =2, chains =2)

#switchting to identity link which is I think(?) what jags is doing
#and then adding better priors to the above
#this is close to Jeff's model but bc p1 is set to 1
# it's equally weighting the snow hazard and the GDD hazard (I think?)
# but what I really need to do (maybe?) - is make hazard not on days_since_snow
# but just on yday (maybe? this would make it so that)
prior4a_better <- 
  prior(uniform(0, 300), lb = 0, ub =300, nlpar = "p2") + 
  prior(uniform(0, 40), lb = 0, ub =40,nlpar = "p3") + 
  prior(uniform(-2, 2), lb = -2, ub =2, nlpar = "beta1") 

#sweet - this does fit and you can see
# the positive forcing effect of GDDs
fit4_better_priors <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
                             p2 + p3 + beta1 ~ 1, nl = TRUE),
                          data = phen_data_sos, family = c("poisson", "identity"), prior = prior4a_better,
                          cores =2)

#try getting p1 in there?
prior4a_better_wp1 <- 
  prior(uniform(0, 1), lb = 0, ub =1, nlpar = "p1") + 
  prior(uniform(0, 300), lb = 0, ub =300, nlpar = "p2") + 
  prior(uniform(0, 40), lb = 0, ub =40,nlpar = "p3") + 
  prior(uniform(-2, 2), lb = -2, ub =2, nlpar = "beta1") 

#sweet - this does fit and you can see
# the positive forcing effect of GDDs
# consider adding in the running 7d mean moisture prior
# and or adding in the yday
#get a few divergent transitions but not too bad
#but notably if you do this you can see that the est
# for p1 is basically 0
fit4_better_priors_wp1 <- brm(bf(status ~ (p1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
                                 p1+p2 + p3 + beta1 ~ 1, nl = TRUE),
                              data = phen_data_sos, family = c("poisson", "identity"), prior = prior4a_better_wp1,
                              cores =2, chains =2)


forcing = function(p1, p2, p3, x){
  return(p1/(1+exp((2.2/p3) * (p2-x))))
}

#p1 - max
# p2 = point of inflection
# p3 = time between 0.1 and 0.5 and beween 0.5 and 0.9

tst = phen_data_sos %>%
  mutate(pred = forcing(days_since_snow, p1 =0.01, p2 =6.67, p3 =2.03))%>%
  mutate(pred1 = pred*exp(0.01*GDD))

#this plots GDD accumulation until event
#or hazard accumulation as a function of doy
ggplot (tst)+
  geom_line(aes(x=days_since_snow,y=pred, group=subject), color = 'black')+
  geom_line(aes(x=days_since_snow,y=pred1, group=subject), color = 'red')+
  #geom_line(aes(x=yday, y=pred, group = subject), color = 'red')+
  facet_wrap(~subject)

ggplot (tst %>%
          filter(days_since_snow>0))+
  geom_line(aes(x=yday,y=pred, group=subject), color = 'black')+
  geom_line(aes(x=yday,y=pred1, group=subject), color = 'red')+
  #geom_line(aes(x=yday, y=pred, group = subject), color = 'red')+
  facet_wrap(~subject)

#try to get yday in there too?
#try getting p1 in there?
prior4a_better_wp1_w_yday <- 
  prior(uniform(0, 1), lb = 0, ub =1, nlpar = "p1") + 
  prior(uniform(0, 300), lb = 0, ub =300, nlpar = "p2") + 
  prior(uniform(0, 40), lb = 0, ub =40,nlpar = "p3") + 
  prior(uniform(-2, 2), lb = -2, ub =2, nlpar = "beta1") +
  prior(uniform(-2, 2), lb = -2, ub =2, nlpar = "beta2") 

#this one doesn't converge, maybe try it after filtering just to days after snowmelt
fit4_better_priors_wp2_yday <- brm(bf(status ~ (p1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD+ beta2*yday),
                                      p1+p2 + p3 + beta1 + beta2 ~ 1, nl = TRUE),
                                   data = phen_data_sos, family = c("poisson", "identity"), prior = prior4a_better_wp1_w_yday,
                                   cores =2, chains =2,
                                   control = list(adapt_delta = 0.99))

#could do by daylength 173 is solstice tho
#this doesn't converge, though
summary(fit4_better_priors_wp2_yday) 


#this has 3 divergent transitions so may or may not converge in the end
# also wide uncertainty on the p parameters
# which don't really converge so I think it can't fit that much
fit4_better_priors_wp2_yday_only_after_snow <- brm(bf(status ~ (p1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD+ beta2*yday),
                                                      p1+p2 + p3 + beta1 + beta2 ~ 1, nl = TRUE),
                                                   data = phen_data_sos %>% filter(days_since_snow>0),
                                                   family = c("poisson", "identity"), prior = prior4a_better_wp1_w_yday,
                                                   cores =2, chains =2,
                                                   control = list(adapt_delta = 0.99))

summary(fit4_better_priors_wp2_yday_only_after_snow) 
pairs(fit4_better_priors_wp2_yday_only_after_snow)



#temp_7d
#this runs and we can see the
# temp and yday effects both
# but the days since snow
# params aren't solving that well
fit4_better_priors_wp2_yday_only_after_snow_rollingtemp <- brm(bf(status ~ (p1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*temp_7d+ beta2*yday),
                                                                  p1+p2 + p3 + beta1 + beta2 ~ 1, nl = TRUE),
                                                               data = phen_data_sos %>% filter(days_since_snow>0),
                                                               family = c("poisson", "identity"), prior = prior4a_better_wp1_w_yday,
                                                               cores =2, chains =2,
                                                               control = list(adapt_delta = 0.99))

summary (fit4_better_priors_wp2_yday_only_after_snow_rollingtemp)
pairs(fit4_better_priors_wp2_yday_only_after_snow_rollingtemp)

#sce thinks this is right but need to double check
# the math against here
#https://data.princeton.edu/wws509/r/c7s1
tst = phen_data_sos %>%
  filter(days_since_snow>0) %>%
  mutate(pred = forcing(days_since_snow, p1 =0.49, p2 =152.08, p3 =34.13 ))%>%
  mutate(pred1 = pred*exp(0.29*temp_7d + 0.03*yday)) %>%
  group_by(subject, yday) %>%
  arrange(yday) %>%
  summarize(cumhaz = cumsum(pred1)) %>%
  mutate(surv = exp(-cumhaz)) %>%
  arrange(subject, yday) %>%
  mutate(risk = 1-surv) %>% #risk = complement of surv
  mutate(pevent = c(NA, diff(risk, 1)))

ggplot (tst, aes(x=yday, y=pevent))+
  geom_line()+
  facet_wrap(~subject)

tst = phen_data_sos %>%
  filter(days_since_snow>0) %>%
  #mutate(pred = forcing(days_since_snow, p1 =0.49, p2 =152.08, p3 =34.13 ))%>%
  mutate(pred1 = exp(0.29*temp_7d + 0.03*yday))

#maybe kind of doing the right thing? Not sure
# it is making the hazard go up steeper when it's further along in the season
ggplot (tst, aes(x=yday, y=pred1))+
  geom_line()+
  facet_wrap(~subject)

#this plots GDD accumulation until event
#or hazard accumulation as a function of doy
ggplot (tst)+
  geom_line(aes(x=days_since_snow,y=pred, group=subject), color = 'black')+
  #geom_line(aes(x=days_since_snow,y=pred1, group=subject), color = 'red')+
  #geom_line(aes(x=yday, y=pred, group = subject), color = 'red')+
  facet_wrap(~subject)

#need to figure out the actual preds


ggplot (tst %>%
          filter(days_since_snow>0))+
  geom_line(aes(x=yday,y=pred, group=subject), color = 'black')+
  geom_line(aes(x=yday,y=pred1, group=subject), color = 'red')+
  #geom_line(aes(x=yday, y=pred, group = subject), color = 'red')+
  facet_wrap(~subject)




#could do by daylength 173 is solstice tho
#this doesn't converge, though
summary(fit4_better_priors_wp2_yday)

#try something where the counting process just starts at snowmelt
#and accelerates faster if (1) it's warmer or (2) it's later in the season

prior4_aftersnow_rolling_temp <- 
  prior(uniform(-2, 2), lb = -2, ub =2, nlpar = "beta1") +
  prior(uniform(-2, 2), lb = -2, ub =2, nlpar = "beta2") 
# 
fit4_after_snow_rollingtemp <- brm(bf(status ~ exp(beta1*temp_7d+ beta2*yday),
                                      beta1 + beta2 ~ 1, nl = TRUE),
                                   data = phen_data_sos %>% filter(days_since_snow>0),
                                   family = c("poisson", "identity"), prior = prior4_aftersnow_rolling_temp ,
                                   cores =2, chains =2,
                                   control = list(adapt_delta = 0.99))

summary (fit4_after_snow_rollingtemp, digits =4)
pairs(prior4_aftersnow_rolling_temp)
plot (conditional_effects(fit4_after_snow_rollingtemp))

pp_check(fit4_after_snow_rollingtemp)

tst = predict(fit4_after_snow_rollingtemp, type = response) %>%
  cbind(., phen_data_sos %>% filter(days_since_snow>0)) %>%
  mutate(pred = exp(0.38*temp_7d+ -0.03*yday))

# I expect this is not quite right
# bc we need to do it invidual-wise (summing log-likelihood)
#over indivs but maybe ok just to check
#smaller looic is the better model, so 
# here it like the more complex model (maybe?)
loo(fit4_after_snow_rollingtemp, fit4_better_priors_wp2_yday_only_after_snow_rollingtemp)

#need a model that describes the fact that
# 1. things that melt out later have generally later phenology
# 2. but telescoping such that days between snowmelt and event are


fit4a <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
                p2 + p3 + beta1 ~ 1, nl = TRUE),
             data = phen_data_sos %>%
               filter(days_since_snow>0), family = "poisson", prior = prior4,
             cores =2, control = list(adapt_delta = 0.99))

#try just seeing if you make it the case that the snowfree can
# only
#sce should look at the code for the tree stuff
# to try to remember how to make these as
# a function of node, or year, or whatever.
fit4a2 <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
                 p2 + p3 + beta1 ~ 1, nl = TRUE),
              data = phen_data_sos %>%
                filter(days_since_snow>0), family = "poisson", prior = prior4a,
              cores =2, control = list(adapt_delta = 0.99))

fit4a2 <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
                 p2 + p3 + beta1 ~ 1, nl = TRUE),
              data = phen_data_sos %>%
                filter(days_since_snow>0), family = "poisson", prior = prior4a,
              cores =2, control = list(adapt_delta = 0.99))

#trying this with p2 and p3 (rate and place it maximizes)
# as a function of sensornode
# probably won't work...
#if you do this without the identity link
# then doesn't fit
# also doesn't fit with the ranefs but let's just try without
fit4a3 <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
                 p2 ~ (1|sensornode),
                 p3 ~ (1|sensornode),
                 beta1 ~ 1, nl = TRUE),
              data = phen_data_sos %>%
                filter(days_since_snow>0), family = c("poisson", "identity"), prior = prior4a,
              cores =2, control = list(adapt_delta = 0.999))

fit4a3_simple <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-days_since_snow))))*exp(beta1*GDD),
                        p2 + p3 + beta1 ~ 1, nl = TRUE),
                     data = phen_data_sos %>%
                       filter(days_since_snow>0), family = c("poisson", "identity"), prior = prior4a,
                     cores =2, control = list(adapt_delta = 0.999))


# other option is just to make it a saturating function of GDD
# instead of days-since-snow? would that work?
prior4_4 <- 
  prior(uniform(0, 300), lb = 0, ub =300, nlpar = "p2") + 
  prior(uniform(0, 40), lb = 0, ub =40, nlpar = "p3")


#this is down to just 70 divergent trans w 0,300 for the bounds
# for each
#still busted need fewer widths

#sce repeat - does this work (maybe)
#does it make sense (not sure>)
fit4a4 <- brm(bf(status ~ (1/(1+exp((2.2/p3) * (p2-GDD)))),
                 p2 + p3  ~ 1, nl = TRUE),
              data = phen_data_sos %>%
                filter(days_since_snow>0), family = c("poisson", "identity"),
              prior = prior4_4,
              cores =2, control = list(adapt_delta = 0.999))


#explorations with segmented - don't
# really seem to make sense
library (segmented)
eos_all_seg = glm(formula = status ~ GDD + 
                    days_since_snow + 
                    moist_7d + 
                    yday + 
                    temp_7d,
                  family = binomial(link = "logit"),
                  data = phen_data_eos %>%
                    ungroup %>%
                    filter(days_since_snow>0))

#try to fit each with a number of breakpoints
o<-segmented(eos_all_seg ,seg.Z=~GDD + 
               days_since_snow + 
               moist_7d + 
               yday + 
               temp_7d)

plot.segmented(o, "GDD")
plot.segmented(o, "days_since_snow")
plot.segmented(o, "moist_7d")
plot.segmented(o, "yday")

slope(o)

surv  = 
  phen_data_pred %>%
  ungroup %>%
  filter(days_since_snow>0) %>%
  mutate(hazard = c(predict(o, newdata =.,
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


summary(o,short=TRUE)

#is the yday just keeping it from being too early?
range (phen_data_eos$doy_event)
plot.segmented(o, "x", .coef=fixef(o), conf.level=.95)


#library(mgcViz) would be good but won't install on ruderalis
sos_all_gam = mgcv::gam(formula = status ~ s(GDD) + s(days_since_snow) +
                          s(moist_7d) + s(yday) + s(temp_7d),
                        family = binomial(link = "logit"),
                        data = phen_data_sos %>%
                          ungroup %>%
                          filter(days_since_snow>0))

#can see more or less that moisture is the only thing that
# has a nonlinear (threshold) effect
eos_all_gam = mgcv::gam(formula = status ~ s(GDD) + s(days_since_snow) +
                          s(moist_7d) + s(yday) + s(temp_7d),
                        family = binomial(link = "logit"),
                        data = phen_data_eos %>%
                          ungroup %>%
                          filter(days_since_snow>0))

eos_all_gam3 = mgcv::gam(formula = status ~ GDD + days_since_snow +
                           s(moist_7d) + yday + temp_7d,
                         family = binomial(link = "logit"),
                         data = phen_data_eos %>%
                           ungroup %>%
                           filter(days_since_snow>0))


##gompertz model - trial and error is very similar to logit so
# used logit as more numerically stable

Gompertz_Model_Baseline <- glm(formula = status ~ GDD + days_since_snow,
                               family = binomial(link = "cloglog"),
                               data = phen_data_sos %>%
                                 ungroup %>%
                                 filter(days_since_snow>0))

Gompertz_Model_Baseline_eos <- glm(formula = status ~ GDD + days_since_snow,
                                   family = binomial(link = "cloglog"),
                                   data = phen_data_eos %>%
                                     ungroup %>%
                                     filter(days_since_snow>0))


#could try them one at a time but it's slow
sos_GDD = glm(formula = status ~ GDD,
              family = binomial(link = "logit"),
              data = phen_data_sos %>%
                ungroup %>%
                filter(days_since_snow>0))

sos_GDD_yday = glm(formula = status ~ GDD + yday,
                   family = binomial(link = "logit"),
                   data = phen_data_sos %>%
                     ungroup %>%
                     filter(days_since_snow>0))

#can't really put GDD and days_since_snow sensibly in the same model
sos_GDD_snowmelt = glm(formula = status ~ GDD + days_since_snow,
                       family = binomial(link = "logit"),
                       data = phen_data_sos %>%
                         ungroup %>%
                         filter(days_since_snow>0))


sos_yday = glm(formula = status ~ yday,
               family = binomial(link = "logit"),
               data = phen_data_sos %>%
                 ungroup %>%
                 filter(days_since_snow>0))

sos_yday_temp7 = glm(formula = status ~ yday + temp_7d,
                     family = binomial(link = "logit"),
                     data = phen_data_sos %>%
                       ungroup %>%
                       filter(days_since_snow>0))

sos_yday_snowmelt = glm(formula = status ~ yday + days_since_snow,
                        family = binomial(link = "logit"),
                        data = phen_data_sos %>%
                          ungroup %>%
                          filter(days_since_snow>0))


sos_GDD_yday = glm(formula = status ~ GDD + yday,
                   family = binomial(link = "logit"),
                   data = phen_data_sos %>%
                     ungroup %>%
                     filter(days_since_snow>0))

sos_snowmelt = glm(formula = status ~ days_since_snow,
                   family = binomial(link = "logit"),
                   data = phen_data_sos %>%
                     ungroup %>%
                     filter(days_since_snow>0))

sos_temp7 = glm(formula = status ~ temp_7d,
                family = binomial(link = "logit"),
                data = phen_data_sos %>%
                  ungroup %>%
                  filter(days_since_snow>0))

sos_moist7 = glm(formula = status ~ moist_7d,
                 family = binomial(link = "logit"),
                 data = phen_data_sos %>%
                   ungroup %>%
                   filter(days_since_snow>0))

#univariate, GDD is the best
AIC (sos_GDD, sos_yday, sos_snowmelt, sos_temp7, sos_moist7,
     sos_GDD_yday, sos_GDD_snowmelt, sos_yday_temp7,
     sos_yday_snowmelt)

lm.full <- lm(bwt ~ age + lwt + race.cat + smoke + preterm + ht + ui + ftv.cat, data = lbw)


#if you check carefully you actually want the moist to be
# squared/threshold, but ignoring that for now
eos_all_2 = glm(formula = status ~ GDD + I(GDD^2) +
                  days_since_snow + I(days_since_snow^2) +
                  moist_7d + I(moist_7d^2) +
                  yday + I(yday^2) +
                  temp_7d +I(temp_7d^2) ,
                family = binomial(link = "logit"),
                data = phen_data_eos %>%
                  ungroup %>%
                  filter(days_since_snow>0))






summary (eos_all_2)
model.eos.aic.backward_2 <- step(eos_all_2, direction = "backward", trace = 1)
# gets reduced to just GDD + yday




#define function to cv null model of doy
model_cv_simplest = function(
    phen_data_fit, phen_data_pred, model){
  out = list()
  for (i in unique (phen_data_fit$sensornode)){
    phen_data_fit_i = phen_data_fit %>%
      filter(sensornode != i) %>%
      ungroup %>%
      filter(yday == doy_event)
    phen_data_pred_i = phen_data_pred %>%
      filter(sensornode == i)%>%
      ungroup %>%
      filter(yday == doy_event)
    mod = update (model, data = phen_data_fit_i)
    # surv  = predict(mod, newdata =phen_data_pred_i,
    #                 type = 'response')
    surv = phen_data_pred_i %>%
      mutate(doy_event_pred = c(predict(mod, newdata =.,
                                        type = 'response'))) %>%
      select(subject, doy_event_pred, year) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, doy_event)%>%
                  distinct(),by = c("subject", "year"))
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
      filter(year == i)%>%
      ungroup %>%
      filter(yday == doy_event)
    mod = update (model, data = phen_data_fit_i)
    # surv  = predict(mod, newdata =phen_data_pred_i,
    #                 type = 'response')
    surv = phen_data_pred_i %>%
      mutate(doy_event_pred = c(predict(mod, newdata =.,
                                        type = 'response'))) %>%
      select(subject, doy_event_pred, year) %>%
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

# simple_mod = lm (GDD ~ 1,
#                  data = phen_data_sos %>%
#                    ungroup %>%
#                    filter(yday == doy_event))

#define function just based on GDD past snowmelt
# simple_mod = lm (GDD ~ 1,
#                  data = phen_data_sos %>%
#                    ungroup %>%
#                    filter(yday == doy_event))

#could rewrite this so it's any threshold
model_cv_simple = function(
    phen_data_fit, phen_data_pred, model){
  out = list()
  for (i in unique (phen_data_fit$sensornode)){
    phen_data_fit_i = phen_data_fit %>%
      filter(sensornode != i) %>%
      ungroup %>%
      filter(yday == doy_event)
    phen_data_pred_i = phen_data_pred %>%
      filter(sensornode == i) #%>%
    #select(subject, year) %>%
    #distinct()
    mod = update (model, data = phen_data_fit_i)
    # surv  = predict(mod, newdata =phen_data_pred_i,
    #                 type = 'response')
    surv = phen_data_pred_i %>%
      select(subject, year) %>%
      distinct() %>%
      mutate(GDD_pred = c(predict(mod, newdata =.,
                                  type = 'response'))) %>%
      select(subject, GDD_pred, year) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, GDD, doy_event, yday)%>%
                  distinct(),
                by = c("subject", "year")) %>%
      filter(GDD>GDD_pred) %>%
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
      filter(year == i) #%>%
    #select(subject, year) %>%
    #distinct()
    mod = update (model, data = phen_data_fit_i)
    # surv  = predict(mod, newdata =phen_data_pred_i,
    #                 type = 'response')
    surv = phen_data_pred_i %>%
      select(subject, year) %>%
      distinct() %>%
      mutate(GDD_pred = c(predict(mod, newdata =.,
                                  type = 'response'))) %>%
      select(subject, GDD_pred, year) %>%
      left_join(., phen_data_pred_i %>%
                  select(subject, year, GDD, doy_event, yday)%>%
                  distinct(),
                by = c("subject", "year")) %>%
      filter(GDD>GDD_pred) %>%
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

as.formula(paste(thresh, paste(1), sep = '~'))


#could rewrite this so it's any threshold

# model_cv_simple_2(phen_data_sos, phen_data_sos_pred, GDD, 'GDD')
# 
# model_cv_simple_2(phen_data_sos, phen_data_sos_pred, 'GDD')
# 
# model_cv_simple_2(phen_data_sos, phen_data_sos_pred, 'days_since_snow')




#look at synchrony among years
# all_inc = sumry %>%
#   filter(year>2017) %>%
#   group_by(sensornode) %>%
#   dplyr::tally() %>%
#   filter(n==4) %>%
#   left_join(sumry) %>%
#   filter(year>2017) %>%
#   group_by(year) %>%
#   summarise (mean_sos = mean(sos),
#              sd_sos = sd (sos),
#              mean_pop = mean(pop),
#              sd_pop = sd (pop),
#              mean_eos= mean(eos),
#              sd_eos = sd (eos))
# 
# 
# ggplot (all_inc,
#         aes(x=mean_sos, y=sd_sos))+
#   geom_point()
# 
# ggplot (all_inc,
#         aes(x=mean_pop, y=sd_pop))+
#   geom_point(aes(color = factor(year)))
# 
# ggplot (all_inc,
#         aes(x=mean_eos, y=sd_eos))+
#   geom_point(aes(color = factor(year)))
# 
# 
# 
# ggplot (all_inc,
#         aes(x=year, y=sd_sos))+
#   geom_point()








####sce original if we need it



####Read in data and munge a bit to proper format####
# met data comes from script 1b_infill_met.R
met_data <-read.csv("data/met_data_infilled.csv")


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



