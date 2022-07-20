# Code to plot sensor node locations on UAV imagery
# SCE 6 July 22
# with inspiration from Elisa Van Cleemput's example code

# setup -------------------------------------------------------------------

library(readxl)
library (raster)
library (rgdal)

data_dir = 'data/UAV_clipped_multispec.ow.data/'

# read data ---------------------------------------------------------------


RGBNIR_stack_0814 <- stack(paste0(data_dir, '20170814_MultiB_RGBNIR.tif'))


# Read plot catchment boundary
Niwot_SaddleCatchment <- readOGR(paste0(data_dir, 'UAV_clipped_multispec.ow.data/ClipBoundary_MinRGBNIROverlapExtent.shp'))


#EPSG:4269
Niwot_plots_df <- read_excel(paste0(data_dir, 'NWT_phenocam/nwt_sdl_locations.xlsx'))%>%
  filter(grepl("Saddle Sensor Array", ALT_SITECODE))%>%
  # site 18 no longer used
  filter(!grepl("18", ALT_SITECODE)) %>%
  mutate (node = str_sub(ALT_SITECODE, start = 26))

xy <- Niwot_plots_df[,c("LONGITUDE","LATITUDE")]
Niwot_plots_spdf <- SpatialPointsDataFrame(coords = xy, data = Niwot_plots_df,
                                           proj4string = CRS("+init=epsg:4269"))

Niwot_plots_spdf_SaddleCatchment_UTM <- spTransform(Niwot_plots_spdf, CRS(proj4string(RGBNIR_stack_0814)))


tiff(paste0("ms_plots/UAV_RGB_0814.tiff"),
     width = 5, height = 6, units = "in", res = 300)


plot(Niwot_SaddleCatchment, col="gray40", lwd=4)
plotRGB(RGBNIR_stack_0814, r = 3, g = 2, b = 1, add = TRUE)
plot(Niwot_SaddleCatchment, line="gray40", lwd=4, add = TRUE)
raster::plot(Niwot_plots_spdf_SaddleCatchment_UTM, col=alpha('white', 0.5), pch=16, cex=4, add=T)
text(coordinates(Niwot_plots_spdf_SaddleCatchment_UTM), 
     labels = parse(text = as.character(Niwot_plots_spdf_SaddleCatchment_UTM$node)), cex = 1)
scalebar(200, below = 'meters', divs =4, type = 'bar')
dev.off()