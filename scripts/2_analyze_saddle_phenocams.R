#### setup#####################################################################
library(tidyverse)

### modified phenex code to allow missing data, add in the -1 per original spec
# and allow user-defined weighting
source("functions/FitDoubleLogBeck_modified.R")

#speed up by not regenerating all the plots
make_plots = TRUE
# read data####################################################################
#all_phen comes from script 1a_prep_data.R
all_phen = read.csv("data/phen_clim_all.csv")


#### infill snow days###########################################################
#for each node, calculate the mean value and sd of max_filtered and spline_filtered
# at snowmelt doy
baselines = all_phen %>%
  rowwise() %>%
  filter(doy == snowmelt_doy_infilled) %>%
  ungroup() %>%
  group_by(sensornode) %>%
  summarize(
    baseline_mn_spline = mean(spline_filtered, na.rm = TRUE), baseline_sd_spline =
      sd(spline_filtered, na.rm = TRUE),
    baseline_mx_filtered = mean(max_filtered, na.rm = TRUE), baseline_sd_mx =
      sd(max_filtered, na.rm = TRUE)
  )

#if there is only one year where we have gcc at snowmelt, then sd is undefined
#for that node, so use others nodes
#estimate sd if only 1 year use other nodes
baselines$baseline_sd_spline[is.na(baselines$baseline_sd_spline)] =
  mean(baselines$baseline_sd_spline[!is.na(baselines$baseline_sd_spline)])

#
baselines$baseline_sd_mx[is.na(baselines$baseline_sd_mx)] =
  mean(baselines$baseline_sd_mx[!is.na(baselines$baseline_sd_mx)])

all_phen = left_join(all_phen, baselines)

#subset to node-years where there are at least 50 days of gcc recorded between
#doy 180 and doy 280; discarding some partial years before fitting
nrec = all_phen %>%
  filter(!is.na(spline_filtered) & doy > 180 & doy < 280) %>%
  group_by(sensornode, year) %>%
  summarize(ct = dplyr::n()) %>%
  filter(ct > 50) %>%
  select(sensornode, year) %>%
  distinct()

all_phen = all_phen %>% inner_join(., nrec)

if (make_plots){
#here is the real data
(ggplot(all_phen %>% mutate(year = factor(year)), aes(x = doy)) +
  geom_point(aes(y = max_filtered, group = year, color = year)) +
  geom_vline(aes(xintercept = snowmelt_doy_infilled), color = "black") +
  geom_vline(aes(xintercept = snowfall_doy_infilled), color = "black") +
  geom_hline(aes(yintercept = baseline_mx_filtered), color = "blue") +
  facet_wrap(~ sensornode + year) + 
  ggtitle('gcc max trimmed to snowfree period \n inferred snowmelt/snowfall days in black')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
)%>% 
  ggsave(file='plots/gcc_trimmed_with_infilled_snowmelt.jpg', ., width = 8, height = 15)

(ggplot(all_phen %>% mutate(year = factor(year))%>%filter(year==2020), aes(x = doy)) +
    geom_point(aes(y = max_filtered, group = year, color = year)) +
    geom_vline(aes(xintercept = snowmelt_doy_infilled), color = "black") +
    geom_vline(aes(xintercept = snowfall_doy_infilled), color = "black") +
    geom_hline(aes(yintercept = baseline_mx_filtered), color = "blue") +
    facet_wrap(~ sensornode + year) + 
    ggtitle('gcc max trimmed to snowfree period \n inferred snowmelt/snowfall days in black')+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
)%>% 
  ggsave(file='plots/gcc_trimmed_with_infilled_snowmelt_2020.jpg', ., width = 8, height = 5)
}

#infill based on (1) date less than snowmelt doy
# or 2 if not known, then infill just doys before 118 (earliest)
#for fall, infill after either snowfall doy OR doy 300 if snowfall doy not known

df =
  all_phen %>%
  rowwise() %>%
  mutate(spline_filtered_infilled =
           ifelse(((doy < snowmelt_doy_infilled & !is.na(snowmelt_doy_infilled)) | doy < 118),
    rnorm(1, baseline_mn_spline, baseline_sd_spline), spline_filtered
  )) %>%
  mutate(spline_filtered_infilled = ifelse(doy >= snowfall_doy_infilled,
    rnorm(1, baseline_mn_spline, baseline_sd_spline), spline_filtered_infilled
  )) %>%
  mutate(max_filtered_infilled =
           ifelse(((doy < snowmelt_doy_infilled & !is.na(snowmelt_doy_infilled)) | doy < 118),
    rnorm(1, baseline_mx_filtered, baseline_sd_mx), max_filtered
  )) %>%
  mutate(max_filtered_infilled =
           ifelse(doy >= snowfall_doy_infilled,
                  rnorm(1, baseline_mx_filtered, baseline_sd_mx), max_filtered_infilled)) %>%
  mutate(spline_filtered_infilled =
           ifelse(doy > 300, rnorm(1, baseline_mn_spline, baseline_sd_spline),
    spline_filtered_infilled
  )) %>%
  mutate(max_filtered_infilled =
           ifelse(doy > 300, rnorm(1, baseline_mx_filtered, baseline_sd_mx),
    max_filtered_infilled
  )) %>%
  # other option is just infill in Nov, Dec, and Jan
  mutate(spline_filtered_infilled_winter = ifelse(doy < 60,
    rnorm(1, baseline_mn_spline, baseline_sd_spline), spline_filtered
  )) %>%
  mutate(spline_filtered_infilled_winter = ifelse(doy > 330,
    rnorm(1, baseline_mn_spline, baseline_sd_spline), spline_filtered_infilled_winter
  )) %>%
  mutate(max_filtered_infilled_winter = ifelse(doy < 60,
    rnorm(1, baseline_mx_filtered, baseline_sd_mx), max_filtered
  )) %>%
  mutate(max_filtered_infilled_winter = ifelse(doy > 330,
    rnorm(1, baseline_mx_filtered, baseline_sd_mx), max_filtered_infilled_winter
  )) %>%
  filter(sensornode != 18) #this one flat and weird and only one yar then taken out.

####fit double logistic curve###################################################

fit_beck = function(df, colname, weighting) {
  out = list()
  for (node in unique(df$sensornode)) {
    df1 = df %>% filter(sensornode == node)
    for (yr in unique(df1$year)) {
      df2 <- df1 %>% filter(year == yr)
      if (sum(!is.na(df2[[colname]][df2$doy > 150 & df2$doy < 250])) > 70) {
        wts = as.numeric(!is.na(df2$spline_filtered))
        # wts[wts==0]=0.05
        wts[wts == 0] <- 0.15
        myfit = FitDoubleLogBeck.4(
          df2[[colname]],
          t = df2$doy, weighting = weighting, wts = wts
        )
        out[[paste(node, yr, sep = "_")]] = myfit
      }
    }
  }
  return(out)
}

#compare raw vs infilled any day that is snow dovered (_infilled)
# vs infille only Dec, Jan, Feb (_infilled_winter)
#comparing iteratively reweighted (TRUE)
#vs just one round of weights (FALSE)
# by trial and error I figured out weights of 0.15 for winter seem to
# show most fidelity to the data
max_raw_rwt = fit_beck(df, 'max_filtered', TRUE)
max_raw = fit_beck(df, 'max_filtered', FALSE)
max_infilled_rwt = fit_beck(df, 'max_filtered_infilled', TRUE)
max_infilled = fit_beck(df, 'max_filtered_infilled', FALSE)
max_infilled_winter_rwt = fit_beck(df, 'max_filtered_infilled_winter', TRUE)
max_infilled_winter = fit_beck(df, 'max_filtered_infilled_winter', FALSE)


#define function to get preds
get_preds = function(fitlist, name) {
  preds = lapply(fitlist, "[[", "predicted")
  preds_max = lapply(preds, function(x) {
    df = data.frame(x)
    df$doy = 1:nrow(df)
    names(df)[1] <- name
    return(df)
  }) %>%
    data.table::rbindlist(., idcol = "what") %>%
    separate(what, into = c("sensornode", "year"), sep = "_") %>%
    mutate(
      sensornode = as.numeric(sensornode),
      year = as.numeric(year)
    )
}

# focus on the max for now as most people don't use spline
with_preds = df %>%
  left_join(., get_preds(max_raw, "max_pred")) %>%
  left_join(., get_preds(max_raw_rwt, "max_rwt_pred")) %>%
  left_join(., get_preds(max_infilled, "max_infilled_pred")) %>%
  left_join(., get_preds(max_infilled_rwt, "max_infilled_rwt_pred")) %>%
  left_join(., get_preds(max_infilled_winter, "max_infilled_winter_pred")) %>%
  left_join(., get_preds(max_infilled_winter_rwt, "max_infilled_winter_rwt_pred"))

if (make_plots){
#no infilling, weighted vs not, all terrible
(ggplot(with_preds %>% mutate(year = factor(year)), aes(x = doy)) +
  geom_vline(aes(xintercept = snowmelt_doy_infilled), color = "black") +
  geom_vline(aes(xintercept = snowfall_doy_infilled), color = "black") +
  geom_hline(aes(yintercept = baseline_mx_filtered), color = "grey") +
  geom_point(aes(y = max_filtered, group = year), color = "black", size = 0.5) +
  geom_line(aes(y = max_pred), color = "green") +
  geom_line(aes(y = max_rwt_pred), color = "blue") +
  facet_wrap(~ sensornode + year) +
    ylim(0.25,0.5)+
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
  ggtitle("without infilling snow, weighted (blue), unweighted(green)")) %>%
    ggsave(file='plots/beckfit_max_raw.jpg', ., width = 5, height = 8)
}
#calc rmse of the fitteds to put on the fig
#because the max_filtered has already been trimmed to the sf period
# this will just show the rmse over the vegetation period
# the rmses to the vegetated
# period are better without the iterative reweighting
# so keep those
rmses=with_preds %>%
  filter(!is.na(max_filtered)) %>%
  group_by(year, sensornode) %>%
  summarise(rmse=Metrics::rmse(max_filtered,  max_infilled_pred))

#better rmse off te
rmses_rwt=with_preds %>%
  filter(!is.na(max_filtered)) %>%
  group_by(year, sensornode) %>%
  summarise(rmse=Metrics::rmse(max_filtered,  max_infilled_winter_rwt_pred))

if (make_plots){
#this is the one I think we should use
(ggplot(with_preds %>% mutate(year = factor(year)), aes(x = doy)) +
    geom_vline(aes(xintercept = snowmelt_doy_infilled), color = "black") +
    geom_vline(aes(xintercept = snowfall_doy_infilled), color = "black") +
    geom_hline(aes(yintercept = baseline_mx_filtered), color = "grey") +
    geom_point(aes(y = max_filtered, group = year), color = "black", size = 0.5) +
    #geom_line(aes(y = max_pred), color = "green") +
    geom_line(aes(y = max_infilled_pred), color = "blue") +
    facet_wrap(~ sensornode + year) +
    ylim(0.25,0.5)+
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
    ggtitle("infilled when confident snow-covered, weighted (blue)")+
    geom_label(data=rmses, aes(x=300, y=0.5, label = round(rmse,3)), size=2)) %>%
  ggsave(file='plots/beckfit_max_infilled.jpg', ., width = 8, height = 15)


#or just all years overtop eachother
  
#sce start here make sure the baselines seem right
(ggplot(with_preds %>%
          filter(year >2017) %>%
          mutate(year = factor(year)), aes(x = doy)) +
    #geom_vline(aes(xintercept = snowmelt_doy_infilled), color = "black") +
    #geom_vline(aes(xintercept = snowfall_doy_infilled), color = "black") +
    #geom_hline(aes(yintercept = baseline_mx_filtered), color = "grey") +
    #geom_point(aes(y = max_filtered, group = year), color = "black", size = 0.5) +
    geom_line(aes(y = max_infilled_pred, color = year)) +
    #geom_line(aes(y = max_infilled_rwt_pred), color = "blue") +
    facet_wrap(~ sensornode) +
    ylim(0.25,0.5)+
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
    ggtitle("phenology time series all years"))%>%
    #geom_label(data=rmses, aes(x=300, y=0.5, label = round(rmse,3)), size=2)) %>%
  ggsave(file='plots/beckfit_max_infilled_allyrs_overplotted.jpg', ., width = 8, height = 15)
  
  #sce start here make sure the baselines seem right
  (ggplot(with_preds %>%
            filter(year >2017) %>%
            mutate(year = factor(year)), aes(x = doy)) +
      geom_line(aes(y = max_infilled_rwt_pred, color = year)) +
      #geom_line(aes(y = max_infilled_rwt_pred), color = "blue") +
      facet_wrap(~ sensornode) +
      ylim(0.25,0.5)+
      theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
      ggtitle("phenology time series all years"))%>%
    #geom_label(data=rmses, aes(x=300, y=0.5, label = round(rmse,3)), size=2)) %>%
    ggsave(file='plots/beckfit_max_infilled_allyrs_overplotted_rwt.jpg', ., width = 8, height = 15)
  
(ggplot(with_preds %>% mutate(year = factor(year)), aes(x = doy)) +
    geom_vline(aes(xintercept = snowmelt_doy_infilled), color = "black") +
    geom_vline(aes(xintercept = snowfall_doy_infilled), color = "black") +
    geom_hline(aes(yintercept = baseline_mx_filtered), color = "grey") +
    geom_point(aes(y = max_filtered, group = year), color = "black", size = 0.5) +
    geom_line(aes(y = max_infilled_winter_pred), color = "green") +
    geom_line(aes(y = max_infilled_winter_rwt_pred), color = "blue") +
    facet_wrap(~ sensornode + year) +
    ylim(0.25,0.5)+
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
    ggtitle("infilled just dec, jan, feb, weighted (blue), unweighted(green)")) %>%
  ggsave(file='plots/beckfit_max_infilled_winter.jpg', ., width = 5, height = 8)

}
#########phenometrics###########################################################

#what is what. 
#pop is peak of season position (fitted)
#rau is max rate of decline doy

sumry = lapply(max_infilled, function(y) phenopix::PhenoDeriv(y$predicted)) %>%
  bind_rows(.) %>% data.frame()
sumry$what = names(max_infilled)
sumry = sumry %>%
  separate(what, into = c('sensornode', 'year'), sep = '_')

#same but using the white method rather than derive method for phenometrics
sumry_white = lapply(max_infilled, function(y)
  phenopix::PhenoTrs(y$predicted, approach = "White", trs = 0.5)) %>%
  bind_rows(.) %>%
  data.frame()
sumry_white$what = names(max_infilled)
sumry_white = sumry_white %>%
  separate(what, into = c("sensornode", "year"), sep = "_")

if (make_plots){
#it's pretty much the same if you do trs at 50% of greenness or the derivs
cor (sumry$sos, sumry_white$sos, use="complete.obs")
cor (sumry$eos, sumry_white$eos, use="complete.obs")
cor (sumry$los, sumry_white$los, use="complete.obs")
}
#sn SOS was really not captured at all by the phenocams
# so remove all the start/peak things from node 9 in 2017
# before further analysis
sumry=sumry%>%mutate(
  sos=ifelse(sensornode==9&year==2017, NA, sos),
  los=ifelse(sensornode==9&year==2017, NA, los),
  rsp=ifelse(sensornode==9&year==2017, NA, rsp),
  msp=ifelse(sensornode==9&year==2017, NA, msp),
  pop=ifelse(sensornode==9&year==2017, NA, pop),
  peak=ifelse(sensornode==9&year==2017, NA, peak),
  mgs=ifelse(sensornode==9&year==2017, NA, mgs))

#save these - that's probably
# all this script should do going forward to
# piecemeal it

write.csv (sumry, "data_deriv/phenometrics.csv",
           row.names = FALSE)

#plot with the phenometrics over it
#or just all years overtop eachother

cols <- c("snowmelt" = "black",
  "spring" = "green", "peak" = "blue",
  "fall"="brown", "snowfall" = "black")

if (make_plots){
(ggplot(with_preds %>% mutate(year = factor(year)), aes(x = doy)) +
    geom_line(aes(y = max_infilled_pred), color='grey', size=2) +
    geom_point(aes(y = max_filtered, group = year), color = "black", size = 0.5) +
    geom_vline(aes(xintercept = snowmelt_doy_infilled, color = "snowmelt")) +
    geom_vline(aes(xintercept = snowfall_doy_infilled, color = "snowfall")) +
    geom_vline(data=sumry, 
               aes(xintercept = sos, color = "spring")) +
    geom_vline(data=sumry, 
               aes(xintercept = pop, color = "peak")) +
    geom_vline(data=sumry, 
               aes(xintercept = eos, color = "fall")) +
    facet_grid(year ~ sensornode) +
    ylim(0.34,0.475)+
    scale_colour_manual(name="estimate dates",values=cols)+ 
                       #guide = guide_legend())+
    theme(axis.text.x=element_text(angle = 90, hjust = 1)) +
    geom_label(data=rmses, aes(x=300, y=0.46, label = round(rmse,3)), size=2)+ 
    ylab('Green Chromatic Coordinate (filtered)')+
    ggtitle("phenology time series all years"))%>%
  #geom_label(data=rmses, aes(x=300, y=0.5, label = round(rmse,3)), size=2)) %>%
  ggsave(file='plots/beckfit_max_infilled_allyrs_overplotted_with_phenometrics.jpg', ., width = 25, height = 8)


#sce space for time graph here
#how snowmelt constrains greenup over time
time_pop_peak=ggplot(data=sumry%>%filter(!is.na(pop)&!is.na(peak)),
                aes(x=pop, y=peak, group=sensornode))+
  #geom_smooth(method='lm', aes(color=sensornode), se=FALSE)+
  geom_smooth(method='lm', se=FALSE, color='darkgrey')+
  geom_point(aes(shape=year))+
  geom_abline(slope=1, intercept=20, color='black', size=1, linetype='dashed')+
  labs(x='timing_of_peak', y='magnitude_of_peak')+
  #note some oddity, 21 isn't showing upin lines and is on the odd side
  # geom_line(data=sumry%>%filter(sensornode==21), aes(x=snowmelt_doy_infilled, y=sos,),
  #           color="black")+
  #ggtitle('snowmelt strongly constrains \n interannual variation in greenup')+
  #ggtitle('start of season')+ 
  theme(plot.title = element_text(size = 10, face = "bold"))+
  guides(color = FALSE, shape=FALSE)
}
