# NWT_phenocam
NWT ridge phenocam analysis code

## Overview: 
The basic workflow for running is to go to the scripts directory and run the
scripts in order. Before running the scripts, be sure to make a /plots directory -- this houses some interim figures not used in the final analysis without it the code will have some errors trying to write to non-existing locations. 


## Processing details
- `0_download_data.R` will download 3 data packages from EDI.
  - The phenocam data v2 (knb-lter-nwt.192.3)  
Elwood, K., J. Smith, S. Elmendorf, and Niwot Ridge LTER. 2022.
Time-lapse camera (phenocam) imagery of sensor network plots,
2017 - ongoing. ver 3. Environmental Data Initiative.
https://doi.org/10.6073/pasta/285918fbf5cc4bd2ed2c1241db9a1b2d (Accessed 2022-05-23).

  - The saddle sensor node met data v3 (knb-lter-nwt.210.4)  
Morse, J. and Niwot Ridge LTER. 2022. Climate data for saddle catchment sensor
network, 2017 - ongoing. ver 4. Environmental Data Initiative.
https://doi.org/10.6073/pasta/598894834ea3bae61d7550c30da06565
(Accessed 2022-05-23).

  - 5cm multispectral UAV imagery (knb-lter-nwt.67.1)
Wigmore, O. and Niwot Ridge LTER. 2021. 5cm multispectral imagery from UAV
campaign at Niwot Ridge, 2017 ver 1. Environmental Data Initiative.
https://doi.org/10.6073/pasta/a4f57c82ad274aa2640e0a79649290ca (Accessed 2022-07-06). Select files from this need to be unzipped as well outside of the code.

- `1a_prep_data.R`  
Does some initial preprocessing of the phenocam and met data,
removing flagged values from the met data and averaging up to daily values. It
also estimates snowmelt/snowfall doy where not recorded by the cameras as the
first doy after doy 100 where the 5d running mean diurnal temperature range at
5cm exceeds the mean 5cm dtr at recorded snowmelt at that node. Note that snowfall
occurred more synchronously than snowmelt across the nodes, so snowfall dates were assigned by inspection in some  cases. Where node-specific met data were absent due to sensor failures, the snowmelt dates are not infilled, as the ability of the infilling to exactly replicate temperatures around the freezing point was unknown. This script also pads out missing data to 365 days per year. Generates files:

  - data/all_met_no_freezefill.csv
  - data/phen_clim_all.csv'

Note I had experimented with infilling moisture values when soil temperatures
dropped below freezing, but this approach yielded some peculiar results. Some
sources indicate that the soil moisture may be roughly accurate
even a few degrees below 0C, so moisture data where soil temps were below
freezing are left as they were recorded by the instruments. For the phenology
analyses here, moisture values while soil is frozen are largely irrelevant but
for applications where precise overwinter soil temperatures are needed an
alternate treatment of overwinter soil moisture might be desirable.

- `1b_infill_met.R` Infills the daily met data. Met data are infilled sequentially
starting with air temperature, so that infilled values (e.g. air T at a node
can be used to predict missing soilT at that node, etc). The algorithm uses
MICE for infilling the air/soil temperatures, and mtdsi for infilling soil
moisture. Soil moisture shows a high degree of 
temporal autocorrelation, so mtdsi gave much more sensible results for soil moisture than using MICE, which generated unrealistic discontinuities. In general, the algorithms try to infill based on the
relationships in a narrow window (e.g. one month) around the datagaps, but
the windows are extended in some cases where insufficient data to generate
the relationships exist within that gap. This generates the file
data/met_data_infilled.csv, which is used in subsequent
analysis. It also generates 5 intermediary file
all_met_infilled_seq_1...5_no_freeze_fill.rds. These are not used but can be
helpful if you need to rerun this slow-running script in chunks, because each
step takes ages (for future improvement this code could certainly be optimized
to run more efficiently in parallel).


- `2_analyze_saddle_phenocams.R` generates phenometrics from the phenocam data
using a modified double-log beck workflow (sourced from file 
`functions/FitDoubleLogBeck_modified.r`. The basic code for the 
double-log Beck function was derived from the phenopix package and modified as 
follows:
  - 1. use the original Beck equation (a -1 is missing from the phenopix
implementation - though this does not affect the results much)
  - 2. as well as allowing missing data in the time series. 
  - 3. Weight the data with 0.15 weights given to any greenness values that 
were snow-covered at the time of measurement. These values are infilled
based on the sensor-specific GCC values at snowmelt (random draws
from a normal distribution based on the empirical mean +sd for that node).


Generates files data_deriv/phenometrics.csv
 
- `3_phen_surv.R` runs survival models on the derived phenometrics,
using AIC to select the best model for each phenophase and then doing leave-one-out
cross-validation, sequentially dropping nodes and years. Generates some figures
at the end.

- `4_plot_locations.R` makes some nice plots of the locations of the study


- NOTE: plots used in the manuscript are rendered to /ms_plots; others to /plots
