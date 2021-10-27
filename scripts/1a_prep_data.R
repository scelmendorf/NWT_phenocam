# script to pre-process phenocam and met data
# sce 10 Nov 2020
# creates file 'data/phen_clim_all.csv'
# which has phenocam data together with daily climate data
# where it exists

# Setup ########################################################################
#my R version
library (tidyverse)
library (magrittr)

#whether to remake a zillion plots
plot_data = FALSE
path2data = 'data/'


#read and process met data #####################################################

#define function to remove all flagged values
rm_flag_values = function(df){
  for (col in grep('flag', names(df), value = TRUE)){
    col2edit=gsub('^flag_', '', col)
    #reset values below flag to NA
    df[[col2edit]][!is.na(df[[col]])]<-NA
  }
  return(df)
}

#using the fill by rel values method
all_met = list()
nodes = c('01', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15',
          '16', '17', '19', '20', '21')

#set warns to fail in case of parsing errors so you def read them all
options(warn = 2)

#process each node
#first, remove all flagged values
#then calculate daily mean over all variables
#as well as 90th quantile of dts sample (for snow depth)
#as well as diurnal temperature range at 5cm soil depth, which I have defined
#as soiltemp_dtr = 90th quantile - 10th quantile of 5cm temperature, per day

#note if we wanted to use an alternate approach we could follow Jabis 2020
# who did days with >0.5C diel soil temp variations. 

#using no freeze fill
all_met = list()
nodes = c('01', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15',
          '16', '17', '19', '20', '21')

#set warns to fail in case of parsing errors so you def read them all
options(warn = 2)

#process each node
#first, remove all flagged values
#then calculate daily mean over all variables
#as well as 90th quantile of dts sample (for snow depth)
#as well as diurnal temperature range at 5cm soil depth, which I have defined
#as soiltemp_dtr = 90th quantile - 10th quantile of 5cm temperature, per day

#note if we wanted to use an alternate approach we could follow Jabis 2020
# who did days with >0.5C diel soil temp variations. 

for (node in nodes){
  df = readr::read_csv(
    paste0(path2data, 'sn_', node, '_tenminute.jm.data.csv'),
    guess_max = 1000000, na = c('', 'NaN'))
  #remove flagged values
  df = rm_flag_values(df)
  
  all_met[[node]] = df %>%
    mutate(ymd = lubridate::date(date)) %>%
    select(-contains('flag')) %>%
    select(-LTER_site, -local_site) %>%
    #then avg by node and day
    group_by(sensornode, ymd) %>%
    summarise_at(vars(contains('avg')),
                 .funs = list(ave = ~mean(., na.rm = TRUE)),
                 .groups = 'drop') %>%
    full_join(., df %>%
                mutate(ymd = lubridate::date(date)) %>%
                select(-contains('flag')) %>%
                #even though dts is just a sample not an ave, rename so the rest
                #of the logic cors
                #rename(dts_avg = dts_sample) %>%
                select(-LTER_site, -local_site) %>%
                group_by(sensornode, ymd) %>%
                summarise(snow_90 = quantile(dts_sample, 0.9, na.rm = TRUE),
                          .groups = 'drop')) %>%
    full_join(., df %>%
                mutate(ymd = lubridate::date(date)) %>%
                select(-contains('flag')) %>%
                #even though dts is just a sample not an ave, rename so the rest
                #of the logic cors
                #rename(dts_avg = dts_sample) %>%
                select(-LTER_site, -local_site) %>%
                group_by(sensornode, ymd) %>%
                summarise(max_t = quantile(soiltemp_5cm_avg, 0.9, na.rm = TRUE),
                          min_t = quantile(soiltemp_5cm_avg, 0.1, na.rm = TRUE),
                          .groups = 'drop') %>%
                mutate(soiltemp_dtr = max_t-min_t))
}

all_met = all_met %>%
  data.table::rbindlist(.)

#set options back ok for warnings after here
options(warn = 1)

#sn 15 the soil_temp_30 is bad
all_met = all_met %>%
  mutate(soiltemp_30cm_avg_ave = ifelse(sensornode == 15, NA,
                                        soiltemp_30cm_avg_ave))

write.csv(all_met, paste0(path2data, 'all_met_no_freezefill.csv'), row.names = FALSE)

#read and process phenocam data ################################################
#read phenocam data
all_phen = readr::read_csv(
  paste0(path2data, 'nwt_sennet_phenocams_filtered.ke.data.csv'),
  guess_max = 100000, na = c('', 'NaN')) %>%
  rename(ymd = date, sensornode = node)

#read in the known snowmelt dates as visually assessed from phenocams
snowmelt = readr::read_csv(
  paste0(path2phendata, "nwt_sennet_phenocams_snowmelt.ke.data.csv"),
  guess_max = 1000000, na = c('', 'NaN')) %>%
  mutate(snowmelt_date = dplyr::if_else(confidence_snow_off == 'high',
                         lubridate::ymd(snow_off), lubridate::ymd(NA))) %>%
  mutate(snowmelt_doy = lubridate::yday(snowmelt_date)) %>%
  rename(sensornode = node)

#plot phenocam vs various sensor node met parameters ###########################
if (plot_data) {

  #plot temperature vs phenocam GCC
  #figure out scaling factor in order to overplot gcc overtop temperature
  var = 'airtemp_avg_ave'
  var_mean = mean (all_met[[var]], na.rm = TRUE)
  var_sd = sd(all_met[[var]], na.rm = TRUE)

  phen_var = 'max_filtered'
  phen_mean = mean (all_phen[[phen_var]], na.rm = TRUE)
  phen_var_sd = sd (all_phen[[phen_var]], na.rm = TRUE)

  #if we were to scale both, we would subtract each one's mean and add to
  # each sd.
  # so for gcc - need to subtract the mean gcc, divide by the gcc sd
  # then multiply by the met sd and add in the met mean
  coeff = mean(all_phen[["max_filtered"]], na.rm = TRUE) / mean(all_met[[var]],
    na.rm = TRUE)

  # it looks to track both air and soilT pretty well,
  # but not as jaggedy - might better track
  # the 5d running mean of soilT or something?
  # kind of looks somewhat like the moisture
  # areas, spring tracks the increas in soil temperature
  # whereas in the moist meadows it lags?
  # nodes 6(MM), 7(MM), 8(MM), 11(DM3), 12,(MM),
  # 13(WM), 14(DM1), 15(DM3), 16(DM2), 19(SA) spring tracks soil temp
  # nodes 9(DM1), 10(DM1), 21(DM2), 21(DM2), spring lags increase in soil temp.
  (ggplot(all_met, aes(x = ymd)) +
    geom_line(aes(y = airtemp_avg_ave), color = "red") +
    geom_line(aes(y = soiltemp_5cm_avg_ave), color = "orange") +
    geom_line(aes(y = soiltemp_30cm_avg_ave), color = "green") +
    geom_line(
      data = all_phen,
      aes(y = (((max_filtered - phen_mean) / phen_var_sd) * var_sd) + var_mean),
      color = "blue"
    ) +
    scale_y_continuous(
      # Features of the first axis
      name = "airtemp (red)/soiltemp_5(orange)/soiltemp_30(green)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~ (((. - var_mean) / var_sd) * phen_var_sd) + phen_mean,
        name = "max_filtered (blue)"
      )
    ) + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
    facet_wrap(~sensornode)) %>% 
    ggsave(file='plots/raw_data_plots/gcc_vs_all_temps.jpg', ., width = 5, height = 5)
  
  
  (ggplot(all_met, aes(x = ymd)) +
      #geom_line(aes(y = airtemp_avg_ave), color = "red") +
      #geom_line(aes(y = soiltemp_5cm_avg_ave), color = "orange") +
      #geom_line(aes(y = soiltemp_30cm_avg_ave), color = "green") +
      geom_point(aes(y = airtemp_avg_ave, color = commtype),size = 0.5)+
      geom_point(
        data = all_phen,
        aes(y = (((max_filtered - phen_mean) / phen_var_sd) * var_sd) + var_mean),
        color = "blue", size = 0.5
      ) +
      scale_y_continuous(
        # Features of the first axis
        name = "airtemp (red)/soiltemp_5(orange)/soiltemp_30(green)",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~ (((. - var_mean) / var_sd) * phen_var_sd) + phen_mean,
                            name = "max_filtered (blue)"
        )
      ) + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
      facet_wrap(~sensornode))
  
  
  
  
  (ggplot(all_met, aes(x = ymd))+
      geom_point(aes(y = airtemp_avg_ave, color = commtype),size = 0.5)+ #color = "purple", size = 0.5)+
      #geom_point(aes(y = soilmoisture_b_30cm_avg_ave), color = "pink", size = 0.5)+
      geom_point(data = all_phen,
                 aes(y = (((max_filtered-phen_mean)/phen_var_sd)*var_moist_sd)+
                       var_moist_mean),
                 color = "blue", size=0.5)+
      #geom_hline(aes(yintercept = 0.3))+
      scale_y_continuous(
        # Features of the first axis
        name = "soilmoist_5",
        # Add a second axis and specify its features
        sec.axis = sec_axis(~(((.-var_moist_mean)/var_moist_sd)*phen_var_sd)+
                              phen_mean, name="max_filtered (blue)")
      ) + 
      theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
      facet_wrap(~sensornode))
  
  
  
  

  # the lags behind air temp seem more similar among nodes -
  # spring they all lag the increase in airT by a hit
  (ggplot(all_met, aes(x = ymd)) +
    geom_line(aes(y = airtemp_avg_ave), color = "red") +
    # geom_line(aes(y=soiltemp_5cm_avg_ave), color="orange")+
    # geom_line(aes(y=soiltemp_30cm_avg_ave), color="green")+
    geom_line(data = all_phen, aes(y = (((max_filtered - phen_mean) / phen_var_sd) * var_sd) + var_mean), color = "blue") +
    scale_y_continuous(
      # Features of the first axis
      name = "airtemp (red)",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~ (((. - var_mean) / var_sd) * phen_var_sd) + phen_mean,
        name = "max_filtered (blue)"
      )
    ) + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + 
    facet_wrap(~sensornode)) %>%
    ggsave(file='plots/raw_data_plots/gcc_vs_air_temp.jpg', ., width = 5, height = 5)
  
  #Plot greenness vs moilsture
  var_moist = "soilmoisture_b_5cm_avg_ave"
  var_moist_mean = mean(all_met[[var_moist]], na.rm = TRUE)
  var_moist_sd = sd(all_met[[var_moist]], na.rm = TRUE)
  
  all_met = all_met %>%
    mutate(commtype =
             case_when(
               sensornode %in% c(9,10,14, 16, 20, 21, 11,15) ~ 'DRY',
               sensornode %in% c(6,7,8,12) ~ 'MOIST',
               sensornode %in% c(13,17) ~ 'WET',
               sensornode == 19 ~ 'SUBALPINE',
               TRUE ~ 'UNKNOWN'))
  
  (ggplot(all_met, aes(x = ymd))+
    geom_point(aes(y = soilmoisture_b_5cm_avg_ave, color = commtype),size = 0.5)+ #color = "purple", size = 0.5)+
    #geom_point(aes(y = soilmoisture_b_30cm_avg_ave), color = "pink", size = 0.5)+
    geom_point(data = all_phen,
              aes(y = (((max_filtered-phen_mean)/phen_var_sd)*var_moist_sd)+
                        var_moist_mean),
              color = "blue", size=0.5)+
    #geom_hline(aes(yintercept = 0.3))+
    scale_y_continuous(
      # Features of the first axis
      name = "soilmoist_5",
      # Add a second axis and specify its features
      sec.axis = sec_axis(~(((.-var_moist_mean)/var_moist_sd)*phen_var_sd)+
                            phen_mean, name="max_filtered (blue)")
    ) + 
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
    facet_wrap(~sensornode)) %>%
    ggsave(file='plots/raw_data_plots/gcc_vs_moist.jpg', ., width = 8, height = 5)
  
  #plot internnual variability
  # Plot variation over years in greenness
  (ggplot(all_phen %>% mutate(year = factor(year)), aes(x = doy))+
   geom_point(aes(y = max_filtered, group = year, color = year))+
   geom_vline(data = snowmelt %>% mutate(year = factor(year)),
              aes(xintercept = snowmelt_doy, group = year, color = year))+
   theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
   ggtitle ("max gcc (not infilled) for each node in each year \n vertical lines show snowmelt date")+
   facet_wrap(~sensornode)) %>%
    ggsave(file='plots/raw_data_plots/max_gcc_seasonality.jpg', ., width = 8, height = 5)
  
  # Plot variation in soil temperature
  # the soil temps aren't particularly strikingly different between the years
  #though they do warm up a bit earlier in 2018 than 2019 - it really
  # depends where you look, though clearly the daily temps warmer at snowmelt in 2018
  #than 2019 in at sn7, 16 21 etc
  (ggplot(all_met %>% 
         mutate(year = factor(lubridate::year(ymd)),
                        doy = lubridate::yday(ymd)),
    aes(x = doy))+
    geom_point(aes(y = soiltemp_5cm_avg_ave, group = year, color = year),
                size = 0.5, alpha = 0.1)+
    geom_smooth(aes(y = soiltemp_5cm_avg_ave, group = year, color = year),
                size = 0.5)+
    geom_vline(data = snowmelt %>% mutate(year = factor(year)),
              aes(xintercept = snowmelt_doy, group = year, color = year))+
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
    ggtitle ("5cm soil temp for each node in each year \n vertical lines show snowmelt date")+
    facet_wrap(~sensornode)) %>%
    ggsave(file='plots/raw_data_plots/soil_temp_seasonality.jpg', ., width = 8, height = 5)

  # Plot variation in air temperature
  # the air temps aren't particularly strikingly different between the years
  #though they do warm up a bit earlier in 2018 than 2019 - it really
  # depends where you look, though clearly the daily temps warmer at snowmelt in 2018
  #than 2019 in at sn7, 16 21 etc

  (ggplot(all_met %>% mutate(year = factor(lubridate::year(ymd)),
                          doy = lubridate::yday(ymd)), aes(x = doy))+
    geom_point(aes(y = airtemp_avg_ave, group = year, color = year),
              size = 0.5, alpha = 0.1)+
    geom_smooth(aes(y = airtemp_avg_ave, group = year, color = year), size = 0.5)+
    geom_vline(data = snowmelt %>% mutate(year = factor(year)),
              aes(xintercept=snowmelt_doy, group = year, color = year))+
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
    ggtitle ("air temp for each node in each year \n vertical lines show snowmelt date")+
    facet_wrap(~sensornode)) %>%
    ggsave(file='plots/raw_data_plots/air_temp_seasonality.jpg', ., width = 8, height = 5)

  # Plot variation in soil moisture
  # the soil moistures are more strikingly different among yrs
  #though 2017 largely missing which isn't so helpful
  #2018 the dry down definitely earlier than 2019.
  ((ggplot(all_met %>% mutate(year = factor(lubridate::year(ymd)),
                          doy = lubridate::yday(ymd)), aes(x = doy))+
    geom_point(aes(y = soilmoisture_b_5cm_avg_ave, group = year, color = year),
              size = 0.5, alpha = 0.5))+
    ggtitle ("soil moisture 5cm for each node in each year \n vertical lines show snowmelt date")+
    facet_wrap(~sensornode)) %>%
    ggsave(file='plots/raw_data_plots/soil_moist5_seasonality.jpg', ., width = 8, height = 5)

  (ggplot(all_met %>% mutate(year = factor(lubridate::year(ymd)),
                          doy = lubridate::yday(ymd)), aes(x = doy))+
    geom_point(aes(y = soilmoisture_b_30cm_avg_ave, group = year, color = year),
              size = 0.5, alpha = 0.5)+
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
    ggtitle ("soil moisture 30 cmfor each node in each year \n vertical lines show snowmelt date")+
    facet_wrap(~sensornode)) %>%
    ggsave(file='plots/raw_data_plots/soil_moist30_seasonality.jpg', ., width = 8, height = 5)

  # Plot variation in snow depth using snow depth sensor - doesn't work
  #the snowdepth sensor pretty much sucks all the time - can see a
  # bit in 15 a comparison among yrs - res there is hardly anything helpful
  # at all.
  (ggplot(all_met %>% mutate(year = factor(lubridate::year(ymd)),
                          doy = lubridate::yday(ymd)), aes(x = doy))+
    geom_point(aes(y = snow_90, group = year, color = year), size = 0.5,
              alpha = 0.5)+
    facet_wrap(~sensornode)+scale_y_continuous(trans = "reverse"))

  options(warn = 1)
  ggplot(all_met %>% mutate(year = factor(lubridate::year(ymd)),
                          doy = lubridate::yday(ymd)), aes(x = ymd))+
    geom_point(aes(y = snow_90, group = year, color = year), size = 0.5, alpha = 0.5)+
    facet_wrap(~sensornode)+scale_y_continuous(trans  =  "reverse")
  
} #end plot data loop

#estimate snowmelt off soil dtr where not observed #############################
#estimate doy of snowfree where not observed

#get some warnings where doys are missing in converting to numeric
# fine to ignore
est_sm = all_met %>%
  mutate(year = lubridate::year(ymd),
         doy = lubridate::yday(ymd)) %>%
  left_join(., snowmelt %>% select(sensornode, year, snowmelt_doy)) %>%
  group_by(sensornode) %>%
  mutate(soiltemp_dtr_ma_2 = slider::slide_dbl(soiltemp_dtr, mean, .before = 2,
                                               .after = 2)) %>%
  ungroup() %>%
  filter(!is.na(snowmelt_doy)) %>%
  group_by(sensornode) %>%
  #figure out mean moving average dtr when that node became snow free
  summarise(soiltemp_dtr_ma_2_avg = mean(soiltemp_dtr_ma_2, na.rm = TRUE)) %>%
  full_join(., all_met %>%
              group_by(sensornode) %>%
              mutate(soiltemp_dtr_ma_2 = slider::slide_dbl(soiltemp_dtr, mean, .before = 2, .after = 2))) %>%
  #assume snow-covered before that MA date
  mutate(snow_status_2 = ifelse(soiltemp_dtr_ma_2<soiltemp_dtr_ma_2_avg, 'snow_covered', 'snow_free')) %>%
  mutate(year = lubridate::year(ymd), doy = lubridate::yday(ymd)) %>%
  #find first day it is snow free after doy 100
  filter(snow_status_2 == 'snow_free'&doy>100) %>%
  group_by(sensornode, year) %>%
  summarise(doy = min(doy)) %>%
  #pull out a few that don't have reasonable data to estimate
  mutate(doy = 
           case_when(year == 2018&sensornode%in%c(6,17) ~ 'NA',
                     sensornode == 1 ~ 'NA',
                     year == 2017 ~ 'NA',
                     TRUE ~ as.character(doy))) %>%
  mutate(doy = as.numeric(doy)) %>%
  mutate(ymd = as.Date(doy-1, origin = paste0(year, '-01-01')))


#est snowfall by the same logic
#ie when 5d ma dtr drops below snow-covered value in fall (after sept 15)
est_sf = all_met %>%
  mutate(year = lubridate::year(ymd),
         doy = lubridate::yday(ymd)) %>%
  left_join(., snowmelt %>% select(sensornode, year, snowmelt_doy)) %>%
  group_by(sensornode) %>%
  mutate(soiltemp_dtr_ma_2 = slider::slide_dbl(soiltemp_dtr, mean, .before = 2, .after = 2)) %>%
  ungroup() %>%
  filter(!is.na(snowmelt_doy)) %>%
  group_by(sensornode) %>%
  summarise(soiltemp_dtr_ma_2_avg = mean(soiltemp_dtr_ma_2, na.rm = TRUE)) %>%
  full_join(., all_met %>%
              group_by(sensornode) %>%
              mutate(soiltemp_dtr_ma_2 = slider::slide_dbl(soiltemp_dtr, mean, .before = 2, .after = 2))) %>%
  mutate(snow_status_2 = ifelse(soiltemp_dtr_ma_2<soiltemp_dtr_ma_2_avg, 'snow_covered', 'snow_free')) %>%
  mutate(year = lubridate::year(ymd), doy = lubridate::yday(ymd)) %>%
  filter(snow_status_2 == 'snow_covered'&doy>258) %>% # first doy snow covered after sept 1
  group_by(sensornode, year) %>%
  summarise(doy = min(doy)) %>%
  mutate(doy = 
           case_when(#year == 2018&sensornode%in%c(6,17) ~ 'NA',
                     #year == 2019&sensornode%in%c(19) ~ 'NA',
                     #sensornode == 1 ~ 'NA',
                     year == 2017 ~ 'NA',
                     TRUE ~ as.character(doy))) %>%
  mutate(doy = as.numeric(doy)) %>%
  mutate(ymd = as.Date(doy-1, origin = paste0(year, '-01-01')))


snowmelt = snowmelt %>%
  left_join(., est_sm %>%rename(snowmelt_doy_2 = doy,
                            snowmelt_date_2 = ymd)) %>%
  mutate(snowmelt_doy_infilled = 
           ifelse(!is.na(snowmelt_date), snowmelt_doy, snowmelt_doy_2)) %>%
  mutate(snowmelt_date_infilled = as.Date(snowmelt_doy_infilled-1,
                                          origin = paste0(year, '-01-01'))) %>%
  left_join(., est_sf %>%rename(snowfall_doy_infilled = doy,
                         snowfall_date_infilled = ymd))

#plot the 5d MA diurnal temperate range agains the real (purple)
# and infilled (red) snowmelt dates
if (plot_data){
  #make df to plot
  dtr_img = all_met %>%
  mutate(year = lubridate::year(ymd),
       doy = lubridate::yday(ymd)) %>%
  left_join(., snowmelt %>%select(sensornode, year, snowmelt_doy)) %>%
  group_by(sensornode) %>%
  mutate(soiltemp_dtr_ma_2 = 
           slider::slide_dbl(soiltemp_dtr, mean, .before = 2, .after = 2)) %>%
  ungroup() %>%
  filter(!is.na(snowmelt_doy)) %>%
  group_by(sensornode) %>%
  summarise(soiltemp_dtr_ma_2_avg=mean(soiltemp_dtr_ma_2, na.rm=TRUE)) %>%
  full_join(., all_met %>%
              group_by(sensornode) %>%
              mutate(soiltemp_dtr_ma_2 =
                       slider::slide_dbl(soiltemp_dtr, mean, .before = 2,
                                         .after = 2))) %>%
  mutate(snow_status_2 = ifelse(soiltemp_dtr_ma_2 < soiltemp_dtr_ma_2_avg,
                                'snow_covered', 'snow_free')) %>%
  mutate(year=lubridate::year(ymd), doy=lubridate::yday(ymd))
  
  #here are the snowmelt dates we def know overplotted soiltemp_dtr
  # the nodes vary a lot in their summer dtrs at 5cm
  #not sure if this is insulation from the veg or differences in moisture
  #or differences in aspect
  #6 and 21 have really low daily temp fluctuations at 5cm
  #9, 12, 14, 10 have pretty high daily temp fluctuations
  # of course, this could be something artefact-y like
  # the sensors being shallower at one vs the other
  
  (ggplot(dtr_img, aes( x= ymd, y=soiltemp_dtr))+
  geom_line()+
  facet_wrap(~sensornode)+
  geom_vline(data = snowmelt, aes(xintercept = snowmelt_date_infilled),
           color="red")+
  geom_vline(data = snowmelt, aes(xintercept = snowfall_date_infilled),
                 color="blue")+    
  geom_vline(data = snowmelt, aes(xintercept = snowmelt_date),
           color = "purple")+
  ggtitle("5 day running mean soil diurnal temperature range \n
          observed snowmelt (purple) infilled snowment (red)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ) %>% 
  ggsave(file='plots/raw_data_plots/soil_dtr_with_infilled_snowmelt.jpg', ., width = 8, height = 5)

}

# trim all on first snowfall date
# by inspection of the phenocam data we can see that 2017 was doy 259
bounds <-snowmelt %>% group_by(year) %>%
  summarize(median_snowfall = median (snowfall_doy_infilled, na.rm = TRUE),
            min_snowmelt = min(snowmelt_doy_infilled, na.rm = TRUE))

#join all phen to snowmelt, remove any points before doy 118
all_phen = all_phen %>%
  left_join(., snowmelt) %>%
  #eyeballed snowfall doy for 2017 off greenness
  mutate(snowfall_doy_infilled = ifelse(year == 2017, 259, snowfall_doy_infilled)) %>%
  #for this one too
  mutate(snowfall_doy_infilled = ifelse(year == 2019&is.na(snowfall_doy_infilled), bounds$median_snowfall[bounds$year==2019], snowfall_doy_infilled)) %>%
  mutate(snowfall_doy_infilled = ifelse(year == 2018&is.na(snowfall_doy_infilled), bounds$median_snowfall[bounds$year==2018], snowfall_doy_infilled)) %>%
  mutate(snowfall_doy_infilled = ifelse(year == 2020&is.na(snowfall_doy_infilled), bounds$median_snowfall[bounds$year==2018], snowfall_doy_infilled)) %>%
  rowwise() %>%
  #2017 we are missing the doy but it definitely must have had snow on doy 118 which is the min over all other years
  mutate(spline_filtered = ifelse(((doy < snowmelt_doy_infilled&!is.na(snowmelt_doy_infilled))|doy < 118), NA, spline_filtered)) %>%
  mutate(spline_filtered = ifelse(doy >= snowfall_doy_infilled, NA, spline_filtered)) %>%
  mutate(max_filtered = ifelse(((doy < snowmelt_doy_infilled&!is.na(snowmelt_doy_infilled))|doy < 118), NA, max_filtered)) %>%
  mutate(max_filtered = ifelse(doy >= snowfall_doy_infilled, NA, max_filtered)) %>%
  select(year, sensornode, doy, max_filtered, spline_filtered,snowfall_doy_infilled, snowmelt_doy_infilled, ymd)


# and infilled (red) snowmelt dates
if (plot_data){
  (ggplot(all_phen, aes( x= doy, y=soiltemp_dtr))+
      geom_line()+
      facet_wrap(~sensornode)+
      geom_vline(data = all_phen, aes(xintercept = snowmelt_doy_infilled),
                 color="red")+
      geom_vline(data = all_phen, aes(xintercept = snowfall_doy_infilled),
                 color="blue")+ 
      ggtitle("5 day running mean soil diurnal temperature range \n
          observed snowmelt (purple) infilled snowment (red)")+
      theme(axis.text.x = element_text(angle = 90, hjust = 1))+
     facet_grid(sensornode ~ year)) %>% 
    ggsave(file='plots/raw_data_plots/soil_dtr_with_infilled_snowmelt_final.jpg', ., width = 8, height = 5)
}

#pad out doys
all_doy = expand.grid(year = unique (all_phen$year), doy = seq(1:365))
all_doy = all_doy %>% full_join(all_phen %>% select(sensornode, year,
                                                  snowmelt_doy_infilled,
                                                  snowfall_doy_infilled) %>% 
                                  distinct())
all_phen = full_join(all_phen, all_doy) %>%
  left_join(., all_met)

write.csv(all_phen, 'data/phen_clim_all.csv')



