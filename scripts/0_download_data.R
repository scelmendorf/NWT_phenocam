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

