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

#coxph
library (survival)

fit_sos <- coxph(Surv(last_day, yday, status) ~
                   GDD + days_since_snow+
                   moist_7d + temp_7d,
             data=phen_data_sos %>%
               ungroup %>%
               filter(days_since_snow>0),
             id = subject)

x <- survfit(fit_sos, newdata = phen_data_sos_pred %>%
               ungroup %>%
               filter(days_since_snow>0),
             id = subject)

#actual predictions
preds = data.frame(pred = summary(x)$table[,"median"])
preds$subject = row.names (preds)
preds = preds %>%
  left_join(., phen_data_sos_pred %>%
              select(subject, doy_event, snowmelt_doy_infilled) %>%
              distinct()) %>%
  mutate(doy_event_pred = snowmelt_doy_infilled+pred)

#6.5 RMSE using coxPH for SOS w forcing off doy
Metrics::rmse(preds$doy_event, preds$doy_event_pred)

#try second method
fit_sos_2 <- coxph(Surv(last_day, days_since_snow, status) ~
                   GDD + yday+
                   moist_7d + temp_7d,
                 data=phen_data_sos %>%
                   ungroup %>%
                   filter(days_since_snow>0) %>%
                   mutate(last_day = days_since_snow-1),
                 id = subject)

x_2 <- survfit(fit_sos_2, newdata = phen_data_sos_pred %>%
               ungroup %>%
               filter(days_since_snow>0) %>%
               mutate(last_day = days_since_snow-1),
             id = subject)

#actual predictions
preds_2 = data.frame(pred = summary(x_2)$table[,"median"])
preds_2$subject = row.names (preds_2)
preds_2 = preds_2 %>%
  left_join(., phen_data_sos_pred %>%
              select(subject, doy_event, snowmelt_doy_infilled) %>%
              distinct()) %>%
  mutate(doy_event_pred = snowmelt_doy_infilled+pred)

#6.52 RMSE using coxPH for SOS w forcing off days_since_snow
Metrics::rmse(preds_2$doy_event, preds_2$doy_event_pred)

##pop
fit_pop <- coxph(Surv(last_day, yday, status) ~
                   GDD + days_since_snow+
                   moist_7d + temp_7d,
                 data=phen_data_pop %>%
                   ungroup %>%
                   filter(days_since_snow>0),
                 id = subject)

x <- survfit(fit_pop, newdata = phen_data_pop_pred %>%
               ungroup %>%
               filter(days_since_snow>0),
             id = subject)

#actual predictions
preds = data.frame(pred = summary(x)$table[,"median"])
preds$subject = row.names (preds)
preds = preds %>%
  left_join(., phen_data_pop_pred %>%
              select(subject, doy_event, snowmelt_doy_infilled) %>%
              distinct()) %>%
  mutate(doy_event_pred = snowmelt_doy_infilled+pred)

#7.62 RMSE using coxPH for pop w forcing off doy
Metrics::rmse(preds$doy_event, preds$doy_event_pred)

#try second method
fit_pop_2 <- coxph(Surv(last_day, days_since_snow, status) ~
                     GDD + yday+
                     moist_7d + temp_7d,
                   data=phen_data_pop %>%
                     ungroup %>%
                     filter(days_since_snow>0) %>%
                     mutate(last_day = days_since_snow-1),
                   id = subject)

x_2 <- survfit(fit_pop_2, newdata = phen_data_pop_pred %>%
                 ungroup %>%
                 filter(days_since_snow>0) %>%
                 mutate(last_day = days_since_snow-1),
               id = subject)

#actual predictions
preds_2 = data.frame(pred = summary(x_2)$table[,"median"])
preds_2$subject = row.names (preds_2)
preds_2 = preds_2 %>%
  left_join(., phen_data_pop_pred %>%
              select(subject, doy_event, snowmelt_doy_infilled) %>%
              distinct()) %>%
  mutate(doy_event_pred = snowmelt_doy_infilled+pred)

#7.278617 RMSE using coxPH for pop w forcing off days_since_snow
Metrics::rmse(preds_2$doy_event, preds_2$doy_event_pred)

##popair

fit_pop_air <- coxph(Surv(last_day, yday, status) ~
                       GDD_air + days_since_snow+
                   moist_7d + temp_7d_air,
                 data=phen_data_pop %>%
                   ungroup %>%
                   filter(days_since_snow>0),
                 id = subject)

x <- survfit(fit_pop_air, newdata = phen_data_pop_pred %>%
               ungroup %>%
               filter(days_since_snow>0),
             id = subject)

#actual predictions
preds = data.frame(pred = summary(x)$table[,"median"])
preds$subject = row.names (preds)
preds = preds %>%
  left_join(., phen_data_pop_pred %>%
              select(subject, doy_event, snowmelt_doy_infilled) %>%
              distinct()) %>%
  mutate(doy_event_pred = snowmelt_doy_infilled+pred)

#6.820685 RMSE using coxPH for pop w forcing off doy w air temps
#so consistent that for pop air is a little better than soil?
Metrics::rmse(preds$doy_event, preds$doy_event_pred)

#try second method
fit_pop_2_air <- coxph(Surv(last_day, days_since_snow, status) ~
                     GDD_air + yday+
                     moist_7d + temp_7d_air +
                       strata(sensornode),
                   data=phen_data_pop %>%
                     ungroup %>%
                     filter(days_since_snow>0) %>%
                     mutate(last_day = days_since_snow-1),
                   id = subject)

x_2 <- survfit(fit_pop_2_air, newdata = phen_data_pop_pred %>%
                 ungroup %>%
                 filter(days_since_snow>0) %>%
                 mutate(last_day = days_since_snow-1),
               id = subject)

zp1 <- cox.zph(fit_pop_2_air)

#actual predictions
preds_2 = data.frame(pred = summary(x_2)$table[,"median"])
preds_2$subject = row.names (preds_2)
preds_2 = preds_2 %>%
  left_join(., phen_data_pop_pred %>%
              select(subject, doy_event, snowmelt_doy_infilled) %>%
              distinct()) %>%
  mutate(doy_event_pred = snowmelt_doy_infilled+pred)

#7.278617 RMSE using coxPH for pop w forcing off days_since_snow w soil temps
# can get down to 6.27 using air temps so that actually
# works oddly better for the peak
Metrics::rmse(preds_2$doy_event, preds_2$doy_event_pred)



####

#eos
fit_eos <- coxph(Surv(last_day, yday, status) ~
                   GDD + days_since_snow+
                   moist_7d + temp_7d,
                 data=phen_data_eos %>%
                   ungroup %>%
                   filter(days_since_snow>0),
                 id = subject)

x <- survfit(fit_eos, newdata = phen_data_eos_pred %>%
               ungroup %>%
               filter(days_since_snow>0),
             id = subject)

#actual predictions
preds = data.frame(pred = summary(x)$table[,"median"])
preds$subject = row.names (preds)
preds = preds %>%
  left_join(., phen_data_eos_pred %>%
              select(subject, doy_event, snowmelt_doy_infilled) %>%
              distinct()) %>%
  mutate(doy_event_pred = snowmelt_doy_infilled+pred)

#8.722784 RMSE using coxPH for eos w forcing off doy
Metrics::rmse(preds$doy_event, preds$doy_event_pred)

#try second method
fit_eos_2 <- coxph(Surv(last_day, days_since_snow, status) ~
                     GDD + yday+
                     moist_7d + temp_7d,
                   data=phen_data_eos %>%
                     ungroup %>%
                     filter(days_since_snow>0) %>%
                     mutate(last_day = days_since_snow-1),
                   id = subject)

x_2 <- survfit(fit_eos_2, newdata = phen_data_eos_pred %>%
                 ungroup %>%
                 filter(days_since_snow>0) %>%
                 mutate(last_day = days_since_snow-1),
               id = subject)

#actual predictions
preds_2 = data.frame(pred = summary(x_2)$table[,"median"])
preds_2$subject = row.names (preds_2)
preds_2 = preds_2 %>%
  left_join(., phen_data_eos_pred %>%
              select(subject, doy_event, snowmelt_doy_infilled) %>%
              distinct()) %>%
  mutate(doy_event_pred = snowmelt_doy_infilled+pred)

#9.143779 RMSE using coxPH for eos w forcing off days_since_snow
Metrics::rmse(preds_2$doy_event, preds_2$doy_event_pred)




tst = phen_data_sos %>%
  ungroup %>%
  filter(days_since_snow>0),
id = subject)
tst$pred = summary(x)$table[,"median"]



fit_pop <- coxph(Surv(last_day, yday, status) ~
                   GDD + days_since_snow+
                   moist_7d + temp_7d,
                 data=phen_data_pop %>%
                   ungroup %>%
                   filter(days_since_snow>0),
                 id = subject)

fit_eos <- coxph(Surv(last_day, yday, status) ~
                   GDD + days_since_snow+
                   moist_7d + temp_7d,
                 data=phen_data_eos %>%
                   ungroup %>%
                   filter(days_since_snow>0),
                 id = subject)

anova (fit_eos)
