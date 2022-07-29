# NWT_phenocam
NWT ridge phenocam analysis code

## Overview: 
The basic workflow for running is to go to the scripts directory and run the
scripts in order. The scripts have not all been tested from a fresh repository
so some subdirectories may need to be made manually/may be missing some library
statements at the top. It is a work in progress


## Processing details
- `0_download_data.R` will download 3 data packages from EDI.
  - The phenocam data v2 (knb-lter-nwt.192.2)  
Elwood, K., J. Smith, and Niwot Ridge LTER. 2021. Time-lapse camera (phenocam)
imagery of Sensor Network plots from 2017 to ongoing ver 2.
Environmental Data Initiative.`

  - The saddle sensor node met data v3 (knb-lter-nwt.210.3)  
Morse, J. and Niwot Ridge LTER. 2021. Saddle catchment
sensor network data, 2017- ongoing. ver 3. Environmental Data Initiative.
https://doi.org/10.6073/pasta/c1a7a58e355112c362d35092071fa1f0 (Accessed 2021-10-22).

  - 5cm multispectral UAV imagery (knb-lter-nwt.67.1)
Wigmore, O. and Niwot Ridge LTER. 2021. 5cm multispectral imagery from UAV
campaign at Niwot Ridge, 2017 ver 1. Environmental Data Initiative.
https://doi.org/10.6073/pasta/a4f57c82ad274aa2640e0a79649290ca (Accessed 2022-07-06).

- `1a_prep_data.R`  
Does some initial preprocessing of the phenocam and met data,
removing flagged values from the met data and averaging up to daily values. It
also estimates snowmelt/snowfall doy where not recorded by the cameras as the
first doy after doy 100 where the 5d running mean diurnal temperature range at
5cm exceeds the mean 5cm dtr at recorded snowmelt at that node. The snowfall
occurred more synchronously than snowmelt, so was done by inspection in some 
cases. Where there are not node-specific met data, the snowmelt dates are
not infilled (I did not have sufficient confidence that infilling the met data
would give an accurate estimate of snowmelt doy). 
Also pads out missing data to 365 days per year. Generates files
data/all_met_no_freezefill.csv
data/phen_clim_all.csv'
Note I had experimented with infilling moisture values when soil temperatures
dropped below freezing, but this approach yielded some peculiar results, and 
after some reading that indicates the soil moisture may be roughly accurate
even a few degrees below 0, I left in those values. They may not be perfectly
accurate, although for the phenology data, moisture values below freezing are 
rarely used since they don't tend to occur in the weeks prior to greenup.

- `1b_infill_met.R Infills the daily met data. Met data are infilled sequentially
starting with air temperature, so that infilled values (e.g. air T at a node
can be used to predict missing soilT at that node, etc). The algorithm uses
MICE for infilling the air/soil temperatures, and mtdsi for infilling soil
moisture, based on experimentation. Soil moisture shows a high degree of 
temporal autocorrelation, so mtdsi gave much more sensible results for this
in experimental runs. In general, the algorithms try to infill based on the
relationships in a narrow window (e.g. one month) around the datagaps, but
the windows are extended in some cases where insufficient data to generate
the relationships exist within that gap. This generates the file
data/all_met_infilled_seq_6_no_freeze_fill.rds, which is used in subsequent
analysis. It also generates 5 intermediary file
all_met_infilled_seq_1...5_no_freeze_fill.rds. These are for sanity making
in case you need to go back and monkey with the script again or want to run 
one step at a time and go do other things with your computer in between, 
because each one takes ages.

- `2_analyze_saddle_phenocams.R` generates phenometrics from the phenocam data
using a modified double-log beck workflow (sourced from file 
`functions/FitDoubleLogBeck_modified.r`. SCE took the basic code for the 
double-log beck from the phenopix package and modified it to
  - 1. use the original Beck equation (a -1 is missing from the phenopix
implementation - though this does not affect the results much)
  - 2. as well as allowing missing data in the time series. 
  - 3. Weight the data with 0.15 weights given to any greenness values that 
were snow-covered at the time of measurement. These values are infilled
based on the sensor-specific GCC values at snowmelt (random draws
from a normal distribution based on the empirical mean +sd for that node).
  - Generates files data_deriv/phenometrics
 
- `3_phen_surv.R` runs survival models on the derived phenometrics,
using AIC to select the best model for each phenophase and then doing leave-one-out
cross-validation, sequentially dropping nodes and years. Generates some figures
at the end.

- `4_plot_locations.R` makes some nice plots of the locations of the study


- plots used in the manuscript are rendered to /ms_plots; others to /plots
