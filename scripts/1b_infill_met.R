# script to infill the met data
# sce 22 Oct 2021
# creates file 'data/all_met_infilled_seq_5_no_freeze_fill.rds'
# which has all the met data post 2017 infilled
# basic logic is to fill temps using MICE
# and moisures using mtsdi
# because moisture shows more temporal autocorrelation

# use same node other sensors and/or other nodes as predictors
# start with nearby time windows to deal with seasonally varying correlations
# , expand windows as nec to get enough
# data to fit

# warning, this takes a LONG time to run.

#todo
#add progress bar
# benchmark to figure out if we could make the parallelizing faster
# it might not be set up quite right

#Setup####
library (tidyverse)
library (GGally)
library(doMC)
library (ggplot2)
registerDoMC(4)

#read in all met for infilling####
all_met = read.csv('data/all_met_no_freezefill.csv') %>%
    mutate(sensornode = paste0('sn_', sensornode))

all_doy = expand.grid(ymd = unique(all_met$ymd),
                      sensornode = unique(all_met$sensornode))

all_met = full_join(all_met, all_doy)

#remove leap yrs
all_met = all_met %>%
  mutate(yday = lubridate::ymd(ymd)%>%lubridate::yday(.))%>%
  filter(yday!=366) %>%
  select(-yday)


fill_func = function (df, missing){
  require (mice)
  df = df %>%
    mutate(ymd = lubridate::ymd(ymd))%>%
    mutate(yday = lubridate::yday(ymd))#%>%
  
  #do not gap fill the 2017 data, too much all missing to deal with
  if (is.null(missing)){
    missing = which (is.na(df$airtemp_avg_ave)&df$ymd>='2018-01-01')
  }
  
  r <- foreach(i=icount(length(missing)), .combine=c, .inorder=TRUE) %dopar%{
    preds = df %>%
      filter(
        ymd<=(df$ymd[missing[i]]+14)&ymd>=(df$ymd[missing[i]]-14))%>%
      pivot_wider(id_cols = ymd, names_from = sensornode,
                  values_from = airtemp_avg_ave) %>%
      left_join(., df %>% 
                  filter(sensornode ==df$sensornode[missing[i]])%>%
                  select(ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave ))
    #if you don't have at least 15d to infill off her, then use all yrs
    if (sum(!is.na(preds[[df$sensornode[missing[i]]]]))<15){
      #wrap around year end
      ydays = seq(df$yday[missing[i]]-14, df$yday[missing[i]]+14)
      ydays[ydays>365]=ydays[ydays>365]-365
      
      preds = df %>%
        filter( yday %in% ydays)%>%
        pivot_wider(id_cols = ymd, names_from = sensornode,
                    values_from = airtemp_avg_ave) %>%
        left_join(., df %>% 
                    filter(sensornode ==df$sensornode[missing[i]])%>%
                    select(ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave ))
      
    }
    
    #if you don't have at least 15d to infill off her, then use all yrs + 30d
    if (sum(!is.na(preds[[df$sensornode[missing[i]]]]))<15){
      #wrap around year end
      ydays = seq(df$yday[missing[i]]-30, df$yday[missing[i]]+30)
      ydays[ydays>365]=ydays[ydays>365]-365
      
      preds = df %>%
        #filter(sensornode !=df$sensornode[aT]&
        filter( yday %in% ydays)%>%
        pivot_wider(id_cols = ymd, names_from = sensornode,
                    values_from = airtemp_avg_ave) %>%
        left_join(., df %>% 
                    filter(sensornode ==df$sensornode[missing[i]])%>%
                    select(ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave ))
      
    }
    
    #if you don't have at least 20d to infill off her, then use all yrs + 45d
    if (sum(!is.na(preds[[df$sensornode[missing[i]]]]))<15){
      #wrap around year end
      ydays = seq(df$yday[missing[i]]-45, df$yday[missing[i]]+45)
      ydays[ydays>365]=ydays[ydays>365]-365
      preds = df %>%
        filter( yday %in% ydays)%>%
        pivot_wider(id_cols = ymd, names_from = sensornode,
                    values_from = airtemp_avg_ave) %>%
        left_join(., df %>% 
                    filter(sensornode ==df$sensornode[missing[i]])%>%
                    select(ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave ))
      
    }
    
    
    #if you don't have at least 20d to infill off her, then use all yrs + 60d
    if (sum(!is.na(preds[[df$sensornode[missing[i]]]]))<15){
      #wrap around year end
      ydays = seq(df$yday[missing[i]]-60, df$yday[missing[i]]+60)
      ydays[ydays>365]=ydays[ydays>365]-365
      preds = df %>%
        filter( yday %in% ydays)%>%
        pivot_wider(id_cols = ymd, names_from = sensornode,
                    values_from = airtemp_avg_ave) %>%
        left_join(., df %>% 
                    filter(sensornode ==df$sensornode[missing[i]])%>%
                    select(ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave ))
      
    }
    
    if (sum(!is.na(preds[[df$sensornode[missing[i]]]]))<15){
      return(NA)
    }else{
      #remove any variables from predictors that
      # are also missing on the doy of interest (otherwise can't predict that day)
      missing_vars = preds %>%
        filter(ymd == df$ymd[missing[i]])
      remove_preds = names (missing_vars)[is.na(missing_vars)]
      #obviously, the var of interest is missing so that's ok
      remove_preds = setdiff(remove_preds, df$sensornode[missing[i]])
      
      col_of_interest =df$sensornode[missing[i]]
      
      #also remove predictors that are missing in >75% of the rows
      # where the col of interest is present
      
      misspred = 
        preds %>%
        filter(!is.na(preds[[col_of_interest]]))%>%
        select(-ymd, -all_of(remove_preds))
      
      nrow_present =apply(misspred, 2, function(x) sum(!is.na(x)))
      frac_present= nrow_present/nrow_present[df$sensornode[missing[i]]]
      keep_vars = names(frac_present)[frac_present>0.25]
      
      mice_out <- mice(preds %>%
                         select(-ymd) %>%
                         select(all_of(keep_vars)),
                       m=5, maxit = 10,
                       method = 'pmm', seed = 500,
                       threshold=0.99999999)
     
      impdatlong <- complete(mice_out, action="long", include=F)
      
      #cannot figure out how to do this correctly without renaming
      infilled = impdatlong %>% 
        rename('col_of_interest' = !!(df$sensornode[missing[i]]))%>%
        group_by(.id)%>%
        summarize(out = mean (col_of_interest))%>%
        mutate (ymd = preds$ymd)
      
      return(infilled$out[infilled$ymd == df$ymd[missing[i]]])
    }
  }
  inf_values = cbind (df[missing,],r) %>%
    rename('aT_infilled' = r)
  df = df %>%left_join(., inf_values%>%select(sensornode, ymd, aT_infilled))
  return (df)
}

###infill soil temp 5cm
#infilled just off the other below-ground params
fill_func_soil_5 = function (df){
  require (mice)
  df = df %>%
    mutate(ymd = lubridate::ymd(ymd))%>%
    mutate(yday = lubridate::yday(ymd))#%>%
  #mutate(sensornode = paste0('sn_', sensornode))
  
  missing = which (is.na(df$soiltemp_5cm_avg_ave)&df$ymd>='2018-01-01')
  
  r <- foreach(i=icount(length(missing)), .combine=c, .inorder=TRUE,
               .packages=c('tidyverse', 'mice')) %dopar%{
                 preds = df %>%
                   filter(
                     ymd<=(df$ymd[missing[i]]+14)&ymd>=(df$ymd[missing[i]]-14)&
                       sensornode ==df$sensornode[missing[i]])%>%
                   select(sensornode, ymd, airtemp_avg_ave, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                          soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                          soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 #if you don't have at least 15d to infill off her, then use all yrs
                 if (sum(!is.na(preds$soiltemp_5cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-14, df$yday[missing[i]]+14)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, airtemp_avg_ave, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 
                 #if you don't have at least 15d to infill off her, then use all yrs + 30d
                 if (sum(!is.na(preds$soiltemp_5cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-30, df$yday[missing[i]]+30)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, airtemp_avg_ave, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 
                 #if you don't have at least 20d to infill off her, then use all yrs + 45d
                 if (sum(!is.na(preds$soiltemp_5cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-45, df$yday[missing[i]]+45)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd,airtemp_avg_ave, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 if (sum(!is.na(preds$soiltemp_5cm_avg_ave))<15){
                   return(NA)
                 }else{
                   missing_vars = preds %>%
                     filter(ymd == df$ymd[missing[i]])
                   remove_preds = names (missing_vars)[is.na(missing_vars)]
                   #obviously, the var of interest is missing so that's ok
                   remove_preds = setdiff(remove_preds, 'soiltemp_5cm_avg_ave')
                   
                   col_of_interest ='soiltemp_5cm_avg_ave'
                   
                   #also remove predictors that are missing in >75% of the rows
                   # where the col of interest is present
                   
                   misspred = 
                     preds %>%
                     filter(!is.na(preds[[col_of_interest]]))%>%
                     select(-ymd, sensornode, -all_of(remove_preds))
                   
                   nrow_present =apply(misspred, 2, function(x) sum(!is.na(x)))
                   frac_present= nrow_present/nrow_present[[col_of_interest]]
                   keep_vars = names(frac_present)[frac_present>0.25]

                   mice_out <- mice(preds %>%
                                      select(all_of(keep_vars)),
                                    m=5, maxit = 10,
                                    method = 'pmm', seed = 500,
                                    threshold=0.99999999)

                   impdatlong <- complete(mice_out, action="long", include=F)
                   
                   #cannot figure out how to do this correctly without renaming
                   infilled = impdatlong %>% 
                     group_by(.id)%>%
                     summarize(out = mean (soiltemp_5cm_avg_ave))%>%
                     mutate (ymd = preds$ymd)
                   
                   return(infilled$out[infilled$ymd == df$ymd[missing[i]]])
                 }
               }
  inf_values = cbind (df[missing,],r) %>%
    rename('soiltemp_5cm_avg_ave_infilled' = r)
  df = df %>%left_join(., inf_values%>%select(sensornode, ymd, soiltemp_5cm_avg_ave_infilled))
  return (df)
}

###infill soil temp 30cm, needs the stuff copied from 5 above
# re col of interest etc
fill_func_soil_30 = function (df){
  require (mice)
  df = df %>%
    mutate(ymd = lubridate::ymd(ymd))%>%
    mutate(yday = lubridate::yday(ymd))#%>%
    missing = which (is.na(df$soiltemp_30cm_avg_ave)&df$ymd>='2018-01-01')

  r <- foreach(i=icount(length(missing)), .combine=c, .inorder=TRUE,
               .packages=c('tidyverse', 'mice')) %dopar%{
                 
                 preds = df %>%
                   filter(
                     ymd<=(df$ymd[missing[i]]+14)&ymd>=(df$ymd[missing[i]]-14)&
                       sensornode ==df$sensornode[missing[i]])%>%
                   select(sensornode, ymd, airtemp_avg_ave, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                          soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                          soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 
                 #if you don't have at least 15d to infill off her, then use all yrs
                 if (sum(!is.na(preds$soiltemp_30cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-14, df$yday[missing[i]]+14)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, airtemp_avg_ave, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 
                 #if you don't have at least 15d to infill off her, then use all yrs + 30d
                 if (sum(!is.na(preds$soiltemp_30cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-30, df$yday[missing[i]]+30)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, airtemp_avg_ave, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 
                 #if you don't have at least 15d to infill off her, then use all yrs + 45d
                 if (sum(!is.na(preds$soiltemp_30cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-45, df$yday[missing[i]]+45)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd,airtemp_avg_ave, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 if (sum(!is.na(preds$soiltemp_30cm_avg_ave))<15){
                   return(NA)
                 }else{
                   missing_vars = preds %>%
                     filter(ymd == df$ymd[missing[i]])
                   remove_preds = names (missing_vars)[is.na(missing_vars)]
                   #obviously, the var of interest is missing so that's ok
                   remove_preds = setdiff(remove_preds, 'soiltemp_30cm_avg_ave')
                   
                   col_of_interest ='soiltemp_30cm_avg_ave'
                   
                   #also remove predictors that are missing in >75% of the rows
                   # where the col of interest is present
                   
                   misspred = 
                     preds %>%
                     filter(!is.na(preds[[col_of_interest]]))%>%
                     select(-ymd, sensornode, -all_of(remove_preds))
                   
                   nrow_present =apply(misspred, 2, function(x) sum(!is.na(x)))
                   frac_present= nrow_present/nrow_present[[col_of_interest]]
                   keep_vars = names(frac_present)[frac_present>0.25]
                   
                   mice_out <- mice(preds %>%
                                      select(all_of(keep_vars)),
                                    m=5, maxit = 10,
                                    method = 'pmm', seed = 500,
                                    threshold=0.99999999)

                   impdatlong <- complete(mice_out, action="long", include=F)
                   
                   #cannot figure out how to do this correctly without renaming
                   infilled = impdatlong %>% 
                     group_by(.id)%>%
                     summarize(out = mean (soiltemp_30cm_avg_ave))%>%
                     mutate (ymd = preds$ymd)
                   
                   return(infilled$out[infilled$ymd == df$ymd[missing[i]]])
                 }
               }
  inf_values = cbind (df[missing,],r) %>%
    rename('soiltemp_30cm_avg_ave_infilled' = r)
  df = df %>%left_join(., inf_values%>%select(sensornode, ymd, soiltemp_30cm_avg_ave_infilled))
  return (df)
}

####
#define function to infill moisture 5cm
fill_func_soil_5_moist = function (df){
  require (mice)
  df = df %>%
    mutate(ymd = lubridate::ymd(ymd))%>%
    mutate(yday = lubridate::yday(ymd))#%>%

  missing = which (is.na(df$soilmoisture_b_5cm_avg_ave))

  r <- foreach(i=icount(length(missing)), .combine=c, .inorder=TRUE,
               .packages=c('tidyverse', 'mice')) %dopar%{
                 
                preds = df %>%
                   filter(
                     ymd<=(df$ymd[missing[i]]+14)&ymd>=(df$ymd[missing[i]]-14)&
                       sensornode ==df$sensornode[missing[i]])%>%
                   select(sensornode, ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                          soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                          soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 
                 #if you don't have at least 15d to infill off her, then use all yrs
                 if (sum(!is.na(preds$soilmoisture_b_5cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-14, df$yday[missing[i]]+14)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 
                 #if you don't have at least 15d to infill off her, then use all yrs + 30d
                 if (sum(!is.na(preds$soilmoisture_b_5cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-30, df$yday[missing[i]]+30)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 
                 #if you don't have at least 15d to infill off her, then use all yrs + 45d
                 if (sum(!is.na(preds$soilmoisture_b_5cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-45, df$yday[missing[i]]+45)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 if (sum(!is.na(preds$soilmoisture_b_5cm_avg_ave))<15){
                   return(NA)
                 }else{
                   mice_out <- mice(preds %>%
                                      select(-ymd, -sensornode),
                                    m=5, maxit = 10,
                                    method = 'pmm', seed = 500,
                                    threshold=0.99999999)

                   impdatlong <- complete(mice_out, action="long", include=F)
                   
                   #cannot figure out how to do this correctly without renaming
                   infilled = impdatlong %>% 
                     group_by(.id)%>%
                     summarize(out = mean (soilmoisture_b_5cm_avg_ave))%>%
                     mutate (ymd = preds$ymd)
                   
                   return(infilled$out[infilled$ymd == df$ymd[missing[i]]])
                 }
               }
  inf_values = cbind (df[missing,],r) %>%
    rename('soilmoisture_b_5cm_avg_ave_infilled' = r)
  df = df %>%left_join(., inf_values%>%select(sensornode, ymd, soilmoisture_b_5cm_avg_ave_infilled))
  return (df)
}

#define function to infill moisture 30cm
fill_func_soil_30_moist = function (df){
  require (mice)
  df = df %>%
    mutate(ymd = lubridate::ymd(ymd))%>%
    mutate(yday = lubridate::yday(ymd))#%>%
  
  missing = which (is.na(df$soilmoisture_b_30cm_avg_ave))#[1:10]

  r <- foreach(i=icount(length(missing)), .combine=c,
               .inorder=TRUE,
               .packages=c('tidyverse', 'mice')) %dopar%{
                
                 preds = df %>%
                   filter(
                     ymd<=(df$ymd[missing[i]]+14)&ymd>=(df$ymd[missing[i]]-14)&
                       sensornode ==df$sensornode[missing[i]])%>%
                   select(sensornode, ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                          soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                          soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 
                 
                 #if you don't have at least 15d to infill off her, then use all yrs
                 if (sum(!is.na(preds$soilmoisture_b_30cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-14, df$yday[missing[i]]+14)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 
                 #if you don't have at least 15d to infill off her, then use all yrs + 30d
                 if (sum(!is.na(preds$soilmoisture_b_30cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-30, df$yday[missing[i]]+30)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 
                 #if you don't have at least 20d to infill off her, then use all yrs + 45d
                 if (sum(!is.na(preds$soilmoisture_b_30cm_avg_ave))<15){
                   #wrap around year end
                   ydays = seq(df$yday[missing[i]]-45, df$yday[missing[i]]+45)
                   ydays[ydays>365]=ydays[ydays>365]-365
                   preds = df %>%
                     filter( yday %in% ydays&
                               sensornode ==df$sensornode[missing[i]])%>%
                     select(sensornode, ymd, soiltemp_30cm_avg_ave, soiltemp_30cm_avg_ave,
                            soilmoisture_a_5cm_avg_ave,soilmoisture_a_30cm_avg_ave,soilmoisture_b_5cm_avg_ave, 
                            soilmoisture_b_30cm_avg_ave,soilmoisture_c_5cm_avg_ave,soilmoisture_c_30cm_avg_ave)
                 }
                 if (sum(!is.na(preds$soilmoisture_b_30cm_avg_ave))<15){
                   return(NA)
                 }else{
                   if (length(unique(preds$soilmoisture_b_30cm_avg_ave[!is.na(preds$soilmoisture_b_30cm_avg_ave)]))>1){
                     mice_out <- mice(preds %>%
                                        select(-ymd, -sensornode),
                                      m=5, maxit = 10,
                                      method = 'pmm', seed = 500,
                                      threshold=0.99999999)

                     impdatlong <- complete(mice_out, action="long", include=F)
                     
                     #cannot figure out how to do this correctly without renaming
                     infilled = impdatlong %>% 
                       group_by(.id)%>%
                       summarize(out = mean (soilmoisture_b_5cm_avg_ave))%>%
                       mutate (ymd = preds$ymd) 
                     return(infilled$out[infilled$ymd == df$ymd[missing[i]]])
                   }
                   else{return(unique(preds$soilmoisture_b_30cm_avg_ave[!is.na(preds$soilmoisture_b_30cm_avg_ave)]))}
                 }
               }
  inf_values = cbind (df[missing,],r) %>%
    rename('soilmoisture_b_30cm_avg_ave_infilled' = r)
  df = df %>%left_join(., inf_values%>%select(sensornode, ymd, soilmoisture_b_30cm_avg_ave_infilled))
  return (df)
}


#30 cm moisture
fill_func_sm_30 = function (df){
  require (mtsdi)
  df = df %>%
    mutate(ymd = lubridate::ymd(ymd))%>%
    mutate(yday = lubridate::yday(ymd))#%>%

  missing = which (is.na(df$soilmoisture_b_30cm_avg_ave))
  #do by node
  missing = unique (df$sensornode[missing])

  r <-     foreach(i= icount(length(missing)),.combine=rbind)%dopar%{
    preds = df %>%
      pivot_wider(id_cols = ymd, names_from = sensornode,
                  values_from = soilmoisture_b_30cm_avg_ave) %>%
      left_join(., df %>% 
                  filter(sensornode ==missing[i])%>%
                  select(ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                         soilmoisture_b_5cm_avg_ave, soilmoisture_a_5cm_avg_ave,
                         soilmoisture_a_30cm_avg_ave, soilmoisture_c_5cm_avg_ave,
                         soilmoisture_c_30cm_avg_ave))%>%
      arrange(ymd)
    
    allmisscols <- apply(preds,2, function(x)all(is.na(x)))  
    colswithallmiss <-names(allmisscols[allmisscols>0])  
    
    preds =preds%>% 
      select(-!!colswithallmiss)
    
    if(!missing[i]%in%names(preds)){return(NULL)}else{
      #try mnimputting it
      sp_control = 7
      pred_mtdsi <- mnimput(preds %>%
                              select(-ymd) %>%
                              as.data.frame(.),
                            data = preds%>%
                              select(-ymd) %>%
                              as.data.frame(.), 
                            maxit = 50, ts = TRUE,
                            sp.control =
                              list(df = rep(sp_control, ncol(preds - 1))))

      pred_mtdsi <- predict(pred_mtdsi) %>%
        as.data.frame(.) %>%
        mutate(ymd = preds %>% 
                 dplyr::pull(ymd),
        )
      
      pred_mtdsi$soilmoisture_b_30cm_avg_ave_infilled = unlist(pred_mtdsi[missing[i]])
      pred_mtdsi$sensornode = missing[i]
      pred_mtdsi = pred_mtdsi %>%select(sensornode, soilmoisture_b_30cm_avg_ave_infilled, ymd)
      return(pred_mtdsi)
    }
  }
  df = df %>%left_join(., r)
  df$soilmoisture_b_30cm_avg_ave_infilled[!is.na(df$soilmoisture_b_30cm_avg_ave)]<-NA
  return(df)
}

fill_func_sm_5 = function (df){
  require (mtsdi)
  df = df %>%
    mutate(ymd = lubridate::ymd(ymd))%>%
    mutate(yday = lubridate::yday(ymd))#%>%
    missing = which (is.na(df$soilmoisture_b_5cm_avg_ave))
  #do by node
  missing = unique (df$sensornode[missing])
  
  r <-     foreach(i= icount(length(missing)),.combine=rbind)%dopar%{
    preds = df %>%
      pivot_wider(id_cols = ymd, names_from = sensornode,
                  values_from = soilmoisture_b_5cm_avg_ave) %>%
      left_join(., df %>% 
                  filter(sensornode ==missing[i])%>%
                  select(ymd, soiltemp_5cm_avg_ave, soiltemp_30cm_avg_ave,
                         soilmoisture_b_30cm_avg_ave, soilmoisture_a_5cm_avg_ave,
                         soilmoisture_a_30cm_avg_ave, soilmoisture_c_5cm_avg_ave,
                         soilmoisture_c_30cm_avg_ave))%>%
      arrange(ymd)
    
    allmisscols <- apply(preds,2, function(x)all(is.na(x)))  
    colswithallmiss <-names(allmisscols[allmisscols>0])  
    
    preds =preds%>% 
      select(-!!colswithallmiss)
    
    #try mnimputting it
    sp_control = 7
    pred_mtdsi <- mnimput(preds %>%
                            select(-ymd) %>%
                            as.data.frame(.),
                          data = preds%>%
                            select(-ymd) %>%
                            as.data.frame(.),

                          maxit = 50, ts = TRUE,
                          sp.control =
                            list(df = rep(sp_control, ncol(preds - 1))))
    pred_mtdsi <- predict(pred_mtdsi) %>%
      as.data.frame(.) %>%
      mutate(ymd = preds %>% 
               dplyr::pull(ymd),
      )
    
    pred_mtdsi$soilmoisture_b_5cm_avg_ave_infilled = unlist(pred_mtdsi[missing[i]])
    pred_mtdsi$sensornode = missing[i]
    pred_mtdsi = pred_mtdsi %>%select(sensornode, soilmoisture_b_5cm_avg_ave_infilled, ymd)
    return(pred_mtdsi)
  }
  df = df %>%left_join(., r)
  df$soilmoisture_b_5cm_avg_ave_infilled[!is.na(df$soilmoisture_b_5cm_avg_ave)]<-NA
  return(df)
}

#process might be better in order, e.g first infill all the airs,
#use that to infill the 5cm soil (and then readjust the 5cm moisture)
#use that to infill the 30 cm soil (and then readjust the 30cm soil moisture)
#then do 5cm soil mositure
#then do 30cm soil moisture

#####
# infill all - mdtsi for the moistures, just mice for the temps
all_met_infilled_seq = fill_func(all_met, missing = NULL)


# if (freeze_fill){
#   saveRDS(all_met_infilled_seq, 'all_met_infilled_seq_1.rds')
# }else{
#   saveRDS(all_met_infilled_seq, 'all_met_infilled_seq_1_no_freeze_fill.rds')
# }


all_met_infilled_seq = all_met_infilled_seq  %>%
  mutate(aT_is_infilled = ifelse (is.na(aT_infilled), FALSE, TRUE))%>%
  mutate(airtemp_avg_ave = ifelse(is.na(airtemp_avg_ave), aT_infilled, airtemp_avg_ave))


#next 5cm soil temp
all_met_infilled_seq = fill_func_soil_5(all_met_infilled_seq)
# if (freeze_fill){
#   saveRDS(all_met_infilled_seq, 'all_met_infilled_seq_2.rds')
# }else{
#   saveRDS(all_met_infilled_seq, 'all_met_infilled_seq_2_no_freeze_fill.rds')
# }

#next 30cm soil temp
all_met_infilled_seq = fill_func_soil_30(all_met_infilled_seq)

all_met_infilled_seq = all_met_infilled_seq %>%
  mutate(soil30_is_infilled = ifelse (is.na(soiltemp_30cm_avg_ave_infilled), FALSE, TRUE))%>%
  mutate(soiltemp_30cm_avg_ave = ifelse(is.na(soiltemp_30cm_avg_ave),
                                        soiltemp_30cm_avg_ave_infilled, soiltemp_30cm_avg_ave))

#2nd to last 5cm moisture mtdsi
all_met_infilled_seq = fill_func_sm_5(all_met_infilled_seq)


all_met_infilled_seq = all_met_infilled_seq %>%
  mutate(soil5_moist_is_infilled = ifelse(is.na(soilmoisture_b_5cm_avg_ave_infilled), FALSE, TRUE))%>%
  mutate(soilmoisture_b_5cm_avg_ave= ifelse(is.na(soilmoisture_b_5cm_avg_ave), soilmoisture_b_5cm_avg_ave_infilled, soilmoisture_b_5cm_avg_ave))


#and 30
#last 5cm moisture mtdsi
all_met_infilled_seq = fill_func_sm_30(all_met_infilled_seq)

all_met_infilled_seq = all_met_infilled_seq %>%
  mutate(soil30_moist_is_infilled = ifelse(is.na(soilmoisture_b_30cm_avg_ave_infilled), FALSE, TRUE))%>%
  mutate(soilmoisture_b_30cm_avg_ave= ifelse(is.na(soilmoisture_b_30cm_avg_ave), soilmoisture_b_30cm_avg_ave_infilled,
                                             soilmoisture_b_30cm_avg_ave))

#sce backup just before fixing the rest
saveRDS(all_met_infilled_seq, "all_met_infilled_seq_5_no_freeze_fill.rds")

all_met_infilled_seq <- all_met_infilled_seq %>%
  mutate(soil30_moist_is_infilled = ifelse(is.na(soilmoisture_b_30cm_avg_ave_infilled), FALSE, TRUE)) %>%
  mutate(soilmoisture_b_30cm_avg_ave = ifelse(is.na(soilmoisture_b_30cm_avg_ave), soilmoisture_b_30cm_avg_ave_infilled,
                                              soilmoisture_b_30cm_avg_ave
  )) %>%
  mutate(soil5_temp_is_infilled = ifelse(is.na(soiltemp_5cm_avg_ave_infilled), FALSE, TRUE)) %>%
  mutate(soiltemp_5cm_avg_ave = ifelse(is.na(soiltemp_5cm_avg_ave), soiltemp_5cm_avg_ave_infilled,
                                       soiltemp_5cm_avg_ave
  )) %>%
  mutate(aT_is_infilled = ifelse(is.na(aT_infilled), FALSE, TRUE)) %>%
  mutate(airtemp_avg_ave = ifelse(is.na(airtemp_avg_ave), aT_infilled,
                                  airtemp_avg_ave)) %>%
  #one more that did not get filled
  mutate(airtemp_avg_ave = ifelse(sensornode == 'sn_14' & ymd =='2020-05-15',
                                  -0.512236, airtemp_avg_ave))%>%
  mutate(
    year = lubridate::year(ymd),
    yday = lubridate::yday(ymd)
  ) 


write.csv(all_met_infilled_seq, "data_deriv/met_data_infilled.csv", row.names = FALSE)

rm(all_met_infilled_seq) #free up memory
###END infilling temps 



