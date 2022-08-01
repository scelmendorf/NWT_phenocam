# script to download source files from EDI
# sce 22 Oct 2021
# resultant files in /data directory 
# updated 5/23/2022 with new data

# Setup ########################################################################

#where data should go
path2data = 'data/'

#note metajam could do this more quickly if it would install more happily on linux
#function to download the EML file from EDI
getEML<-function(packageid){
  require(magrittr)
  myurl<-paste0("https://portal.edirepository.org/nis/metadataviewer?packageid=",
                packageid,
                "&contentType=application/xml")
  myeml<-xml2::read_xml(paste0("https://portal.edirepository.org/nis/metadataviewer?packageid=",
                               packageid,
                               "&contentType=application/xml"))%>%EML::read_eml()
}

#download phenocam####
#copy citation from here:
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nwt.192.3

# Elwood, K., J. Smith, S. Elmendorf, and Niwot Ridge LTER. 2022.
# Time-lapse camera (phenocam) imagery of sensor network plots,
# 2017 - ongoing. ver 3. Environmental Data Initiative.
# https://doi.org/10.6073/pasta/285918fbf5cc4bd2ed2c1241db9a1b2d (Accessed 2022-05-23).

phen_eml = getEML('knb-lter-nwt.192.3')
for (i in (1:length(phen_eml$dataset$dataTable))){
  download.file(url = phen_eml$dataset$dataTable[[i]]$physical$distribution$online$url$url,
                destfile = paste0(path2data, phen_eml$dataset$dataTable[[i]]$physical$objectName)
  )
}

#download met####

#copy citation from here:
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nwt.210.4

#Morse, J. and Niwot Ridge LTER. 2022. Climate data for saddle catchment sensor
#network, 2017 - ongoing. ver 4. Environmental Data Initiative.
#https://doi.org/10.6073/pasta/598894834ea3bae61d7550c30da06565
#(Accessed 2022-05-23).


met_eml = getEML('knb-lter-nwt.210.4')
for (i in (1:length(met_eml$dataset$dataTable))){
  download.file(url = met_eml$dataset$dataTable[[i]]$physical$distribution$online$url$url,
                destfile = paste0(path2data, met_eml$dataset$dataTable[[i]]$physical$objectName)
  )
}

#download imagery ####

#Wigmore, O. and Niwot Ridge LTER. 2021. 5cm multispectral imagery from UAV
#campaign at Niwot Ridge, 2017 ver 1. Environmental Data Initiative.
#https://doi.org/10.6073/pasta/a4f57c82ad274aa2640e0a79649290ca
#(Accessed 2022-07-06).


#increase timeout for downloading large files
options(timeout = max(300, getOption("timeout")))
download.file(url = uav_eml$dataset$otherEntity$physical$distribution$online$url$url,
              destfile = paste0(path2data, uav_eml$dataset$otherEntity$physical$objectName),
              mode = 'wget'
)

#note unix does not like to unzip these in R
#, so these must be unzipped manually on linux before continuing
#UAV_clipped_multispec.ow.data/20170814_MultiB_RGBNIR.tif,
#UAV_clipped_multispec.ow.data/20170814_MultiB_RGBNIR.tfw'
#UAV_clipped_multispec.ow.data/20170814_MultiB_RGBNIR.tif,
#UAV_clipped_multispec.ow.data/20170814_MultiB_RGBNIR.tif.aux.xml,
#UAV_clipped_multispec.ow.data/ClipBoundary_MinRGBNIROverlapExtent.CPG,
#UAV_clipped_multispec.ow.data/ClipBoundary_MinRGBNIROverlapExtent.dbf,
#UAV_clipped_multispec.ow.data/ClipBoundary_MinRGBNIROverlapExtent.prj,
#UAV_clipped_multispec.ow.data/ClipBoundary_MinRGBNIROverlapExtent.sbn,
#UAV_clipped_multispec.ow.data/ClipBoundary_MinRGBNIROverlapExtent.sbx,
#UAV_clipped_multispec.ow.data/ClipBoundary_MinRGBNIROverlapExtent.shp,
#UAV_clipped_multispec.ow.data/ClipBoundary_MinRGBNIROverlapExtent.shp.xml,
#UAV_clipped_multispec.ow.data/ClipBoundary_MinRGBNIROverlapExtent.shx