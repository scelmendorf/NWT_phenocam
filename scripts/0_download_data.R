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

phen_eml = getEML('knb-lter-nwt.192.2')
for (i in (1:length(phen_eml$dataset$dataTable))){
  download.file(url = phen_eml$dataset$dataTable[[i]]$physical$distribution$online$url$url,
                destfile = paste0(path2phendata, phen_eml$dataset$dataTable[[i]]$physical$objectName)
  )
}



