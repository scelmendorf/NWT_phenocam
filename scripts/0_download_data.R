# script to download source files from EDI
# sce 22 Oct 2021
# resultant files in /data directory 

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
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nwt.192.2

#Elwood, K., J. Smith, and Niwot Ridge LTER. 2021. Time-lapse camera (phenocam)
#imagery of Sensor Network plots from 2017 to ongoing ver 2.
#Environmental Data Initiative.
#https://doi.org/10.6073/pasta/89e8189093392325ee139eccc6b2ff85 (Accessed 2021-10-22).

phen_eml = getEML('knb-lter-nwt.192.2')
for (i in (1:length(phen_eml$dataset$dataTable))){
  download.file(url = phen_eml$dataset$dataTable[[i]]$physical$distribution$online$url$url,
                destfile = paste0(path2data, phen_eml$dataset$dataTable[[i]]$physical$objectName)
  )
}

#download met####

#copy citation from here:
#https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-nwt.210.3

#Morse, J. and Niwot Ridge LTER. 2021. Saddle catchment
#sensor network data, 2017- ongoing. ver 3. Environmental Data Initiative.
#https://doi.org/10.6073/pasta/c1a7a58e355112c362d35092071fa1f0 (Accessed 2021-10-22).


met_eml = getEML('knb-lter-nwt.210.3')
for (i in (1:length(met_eml$dataset$dataTable))){
  download.file(url = met_eml$dataset$dataTable[[i]]$physical$distribution$online$url$url,
                destfile = paste0(path2data, met_eml$dataset$dataTable[[i]]$physical$objectName)
  )
}

