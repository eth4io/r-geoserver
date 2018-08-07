# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Copyright 2016-2018 University of Melbourne
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# Written by: Dr. Yiqun Chen    yiqun.c@unimelb.edu.au
# DevLogs:
# v2.8 2018-07-15
# (1) update loadGeoJSON2SP and loadGeoJSON2SPWithAuth methods so it can work with Mac/Linux and windows

# v2.7 2017-08-19
# (1) add utils.loadGeoJSON2SPWithAuth and utils.loadGeoJSON2DFWithAuth methods to access wfs protected auth info
#
# v2.6 2017-08-15
# (1) add utils.project2WGS84 method
#
# v2.4 2017-07-19
# (1) a new function "utils.loadGeoJSON2DF" for handling null geometry geojson urls
# (2) "utils.createWorkspace" and "utils.getWorkspace" functions are created to manipulate workspaces
# (3) a temporay fix for the glitch when calling "utils.createFeatureType" function. In some scenarios, the rest api raises a 500 error while the datalayer can still be published successfully in geoserver. just skip setting procflag = FALSE to make it bit more robust and tolerable. .
#
# v2.3 2017-06-22
# (1) a new function "utils.publishSP2GeoServerWithMultiStyles" to publish one geolayer with multiple styles. This function is particularly useful when developers need to render various attribute of a single output geolayer.
# (2) a new function "utils.project2UTM" to automatically reproject(transform) a sp object into a proper UTM crs.
#
# v2.2 2017-06-16
# (1) a new function "utils.publishSP2GeoServerWithStyle" is created to combine geolayer publish and styling in one place.
# (2) add missing titles in the charts of outputs
#
# v2.1 2017-06-09 
# (1) fix bug (null is returned when warning occurs in "utils.loadGeoJSON2SP" method). A valid sp object will be returned and warning messages will print out.
# (2) add timeout (36000 seconds) for "utils.addShp2DataStore" and "utils.loadGeoJSON2SP" methods
#
# v1.3 2017-05-29 
# (1) add getLayerNameFromWFSUrl and getBbox functions to prepare output values 
#
# v1.2 2017-05-19 
# (1) add df2jsonlist function to convert dataframe to json list(useful in output wrapping)
#
# v1.1 2017-04-19 
# (1) add wms style creation methods
# (2) add utm zone code calculation methods
#
# v1.0 2017-04-06
# (1) add Geoserver interaction methods
# (2) add GeoJSON->sp parsing methods
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(maptools) 
library(rgdal)
library(rgeos)
library(geojsonio)
library(RCurl)
library(uuid)
library(jsonlite)
library(XML)

# this variable contains all credentials for accessing a geoserver instance to publish data layers and create styles
globalGSCredentials = list()

globalGSCredentials["gsRESTURL"] = "http://115.146.85.153:8080/geoserver" #change it to your own geoserver url
globalGSCredentials["gsRESTUSER"] = "admin" #change it to your own geoserver username (admin role required) 
globalGSCredentials["gsRESTPW"] = "nk39qc12" #change it to your own password 
globalGSCredentials["gsWORKSPACENAME"] = "a2z_ws_unique_student_id" #change it to your own unique workspace name, a new workspace will be automatically created if it does not exist  
globalGSCredentials["gsDATASTORESNAME"] = "a2z_ds_unique_student_id" #change it to your own unique datastore name, a new datastore will be automatically created if it does not exist 
globalGSCredentials["tempDirPath"] = sprintf("%s/%s",getwd(),"tempdata")
globalGSCredentials["wfsUrlTemplate"] = "%s/wfs?service=wfs&version=1.0.0&request=GetFeature&typeName=%s:%s&outputFormat=json" #DON'T MODIFY THIS LINE


WMSStyleCreateUrl = "http://apps.csdila.ie.unimelb.edu.au/gs-styling-service/styling/wms/create" #DON'T MODIFY THIS LINE

proj4string_epsg4326 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" #DON'T MODIFY THIS LINE
proj4string_utm_template = "+proj=utm +zone=%i %s+ellps=WGS84 +datum=WGS84 +units=m +no_defs" #DON'T MODIFY THIS LINE

#' debug print
#'
#' @param arg1 A string to be printed
#'
#' @return debug print string started with "==="
#' @export
#'
#' @examples
#' utils.debugprint("all done")
utils.debugprint <- function(arg1){
  
  print(sprintf("=== %s",arg1 ))
  
}


#' load geojson from url into a sp object
#'
#' @param url A geojson url string
#'
#' @return A sp object
#' @export
#'
#' @examples
utils.loadGeoJSON2SP <- function(url){
    
  # create a unique temp file name for geojson
  tmpFilePath = sprintf("%s\\%s.geojson", globalGSCredentials$tempDirPath, UUIDgenerate(FALSE))
  
  # if tempDirPath not existed, create
  if(dir.exists(globalGSCredentials$tempDirPath)==FALSE){
    
    dir.create(globalGSCredentials$tempDirPath, showWarnings=FALSE, recursive=TRUE)
    
    utils.debugprint(sprintf("%s created",globalGSCredentials$tempDirPath))
  }
  
  spobj <- tryCatch(
    {
      # load data from url
      geojson = getURL(url, timeout=36000)
      
      # save data as a local copy
      write(geojson, tmpFilePath)

      # load it as sp object
      #readOGR(tmpFilePath, "OGRGeoJSON")
      readOGR(tmpFilePath, disambiguateFIDs=TRUE)
    },
    error=function(cond) {
      return(NULL)
    },
    finally={
      # remove local temp file
      file.remove(tmpFilePath)
    }
  )    
  return(spobj)
  
}

#' load geojson from url (username/password protected) into a sp object
#'
#' @param url A geojson url string
#' @param username username required to access url
#' @param password password required to access url
#'
#' @return A sp object
#' @export
#'
#' @examples
utils.loadGeoJSON2SPWithAuth <- function(url, username, password){
  
  # create a unique temp file name for geojson
  tmpFilePath = sprintf("%s\\%s.geojson", globalGSCredentials$tempDirPath, UUIDgenerate(FALSE))
  
  # if tempDirPath not existed, create
  if(dir.exists(globalGSCredentials$tempDirPath)==FALSE){
    
    dir.create(globalGSCredentials$tempDirPath, showWarnings=FALSE, recursive=TRUE)
    
    utils.debugprint(sprintf("%s created",globalGSCredentials$tempDirPath))
  }
  
  spobj <- tryCatch(
    {
      # load data from url
      geojson = getURL(url, timeout=36000, userpwd=paste(username, password, sep=":"), httpauth = 1L)
      
      # save data as a local copy
      write(geojson, tmpFilePath)
      
      # load it as sp object
      #readOGR(tmpFilePath, "OGRGeoJSON")
      readOGR(tmpFilePath, disambiguateFIDs=TRUE)
    },
    error=function(cond) {
      return(NULL)
    },
    finally={
      # remove local temp file
      file.remove(tmpFilePath)
    }
  )    
  return(spobj)
  
}

#' load geojson from url into a data.frame object. this is particularly useful to handle geojson wfsurl with null geometry
#'
#' @param url A geojson url string
#'
#' @return A data.frame object
#' @export
#'
#' @examples
utils.loadGeoJSON2DF <- function(url){
  
  dfobj <- tryCatch(
    {
      # load data from url
      geojson = getURL(url, timeout=36000)
      
      tmp = fromJSON(geojson)
      # only return the properties data frame
      tmp$features$properties
    },
    error=function(cond) {
      return(NULL)
    },
    finally={

    }
  )    
  return(dfobj)
  
}

#' load geojson from url (username/password protected) into a data.frame object. this is particularly useful to handle geojson wfsurl with null geometry
#'
#' @param url A geojson url string
#' @param username username required to access url
#' @param password password required to access url
#' 
#' @return A data.frame object
#' @export
#'
#' @examples
utils.loadGeoJSON2DFWithAuth <- function(url, username, password){
  
  dfobj <- tryCatch(
    {
      # load data from url
      geojson = getURL(url, timeout=36000, userpwd=paste(username, password, sep=":"), httpauth = 1L)
      
      tmp = fromJSON(geojson)
      # only return the properties data frame
      tmp$features$properties
    },
    error=function(cond) {
      return(NULL)
    },
    finally={
      
    }
  )    
  return(dfobj)
  
}


#' publish a sp object to geoserver using RESTful API via curl
#'
#' @param spobj A sp object. *IMPORTANT: Each column of spobj MUST NOT have attribute name exceeding length (int)10, failed to do so may result in NO Style being generated on GeoServer*
#'
#' @return A wfs url string of successfully published data layer
#' @export
#'
#' @examples
utils.publishSP2GeoServer <- function(spobj){
  
  procflag = TRUE
  publishedWfsUrl = NULL #if error occurs, return publishedWfsUrl as NULL
  
  # check if sp is null, return null if null
  if(is.null(spobj)){
    return(NULL)
  }

  # save spobj as shp file
  tmpFileName = UUIDgenerate(FALSE)
  tmpFilePath = sprintf("%s\\%s", globalGSCredentials$tempDirPath, tmpFileName)
  dir.create(tmpFilePath, showWarnings=FALSE, recursive=TRUE)
  writeOGR(spobj, dsn=tmpFilePath, layer = tmpFileName,  driver="ESRI Shapefile", check_exists=TRUE, overwrite_layer=TRUE)
  
  # zip shp file
  # ref: for windows, install Rtools from https://cran.r-project.org/bin/windows/Rtools/ and make sure it is on the system Path
  # flag -j to store just the name of a saved file (junk the path), and do not store directory names. By default, zip will store the full path
  zip(zipfile = tmpFilePath, files = dir(tmpFilePath, full.names = TRUE), flags="-j")  

  # upload zip for geoserver
  out = utils.addShp2DataStore(sprintf("%s.zip", tmpFilePath))
  if(nchar(out)>0){
    procflag = FALSE
  }
  
  if(procflag){
    # publish it as new featuretype
    out = utils.createFeatureType(tmpFileName)
    if(nchar(out)>0){
      # it rarely happens (current testing indicates when a new ws and new ds is created at the same time, server might raise a 500 error), it seems that the layer is still published as usual. need more investigation, just skip seeting procflag = FALSE for now
      #procflag = FALSE 
    }
  }
 
  # return wfs url for the uploaded datalayer
  #wfsUrlTemplate: %s/wfs?service=wfs&version=1.0.0&request=GetFeature&typeName=%s:%s&outputFormat=json
  if(procflag){
    publishedWfsUrl = sprintf(globalGSCredentials$wfsUrlTemplate, globalGSCredentials$gsRESTURL,globalGSCredentials$gsWORKSPACENAME, tmpFileName)
  }
  
  # remove tmp zip file
  file.remove(sprintf("%s.zip", tmpFilePath))
  # remove tmp shp file folder
  unlink(tmpFilePath, recursive=TRUE)
  
  return(publishedWfsUrl)
}

#
#' publish a sp object to geoserver with wms stlying information
#'
#' @param spobj a sp object. *IMPORTANT: Each column of spobj MUST NOT have attribute name exceeding length (int)10, failed to do so may result in NO Style being generated on GeoServer* 
#' @param displayname a readable display name for data layer
#' @param attrname the attribute name that the style will be created for
#' @param palettename check colors defined under each palette name at http://colorbrewer2.org 
#'  acceptable color PaletteNames are: YlOrRd: 9 colors PRGn: 11 colors
#'  PuOr:11 colors RdGy: 11 colors Spectral: 11 colors Grays: 9 colors PuBuGn: 9 colors RdPu: 9 colors BuPu: 9 colors YlOrBr: 9 colors Greens: 9 colors
#'  BuGn: 9 colors Accents: 8 colors GnBu: 9 colors PuRd: 9 colors Purples: 9 colors RdYlGn: 11 colors Paired: 12 colors Blues: 9 colors RdBu: 11 colors 
#'  Oranges: 9 colors RdYlBu: 11 colors PuBu: 9 colors OrRd: 9 colors Set3: 12 colors Set2: 8 colors Set1: 9 colors Reds: 9 colors PiYG: 11 colors 
#'  Dark2: 8 colors YlGn: 9 colors BrBG: 11 colors YlGnBu: 9 colors
#'  Pastel2: 8 colors Pastel1: 9 colors
#' 
#' @param colorreverseorder true or false(default), whether to reverse the color order defined in a palette
#' @param geomtype acceptable geomtypes are: point, multipoint, polygon, multipolygon, linestring, multilinestring, geometry
#' @param colornum the number of colors to be classified. default 5.  in range of 2-15
#' @param classifier acceptable classifier function names are: Jenks(default), EqualInterval, Quantile, StandardDeviation 
#'
#' @return a list contains single element to be included in geolayers
#' @export
#'
#' @examples
utils.publishSP2GeoServerWithStyle <- function(spobj, 
                                        attrname, 
                                        palettename="Reds", 
                                        colorreverseorder=FALSE, 
                                        geomtype="Geometry", 
                                        colornum=5, 
                                        classifier="Jenks")
  {
  
  # TODO: this is the place to add any additional column (i.e. normalised values) to spobj 
  
  # make sure spobj is projected in EPSG:4326 (WGS84) before publishing it. 
  spobj = spTransform(spobj,CRS(proj4string_epsg4326))
  
  # publish sp to geoserver
  outputWfsUrl = utils.publishSP2GeoServer(spobj)
  if(is.null(outputWfsUrl)){
    return(NULL)
  }
  
  # calculte wms style 
  wmsStyleResult = fromJSON(utils.createWMSStyle(wfsurl=outputWfsUrl,
                                                 attrname=attrname,
                                                 palettename=palettename,
                                                 colorreverseorder=colorreverseorder,
                                                 geomtype=geomtype,
                                                 colornum=colornum,
                                                 classifier=classifier
                                                 ))
  wmsStyleparams = list()
  if(wmsStyleResult$status==0){
    wmsStyleparams = wmsStyleResult$data
  }
  
  
  # create content of the 1st element

  geolayerinfo = list(
    layername = utils.getLayerNameFromWFSUrl(wfsurl=outputWfsUrl),
    bbox = utils.getBbox(spobj),
    wfsurl = outputWfsUrl,
    wmsstyle = wmsStyleparams
  )
  
  return(list(geolayerinfo))
  
}


#' publish a sp object to geoserver with multiple wms stlying information
#'
#' @param spobj a sp object. *IMPORTANT: Each column of spobj MUST NOT have attribute name exceeding length (int)10, failed to do so may result in NO Style being generated on GeoServer*
#' @param displayname_vec a vector of readable display names for data layer
#' @param attrname_vec a vector of attribute names that the style will be created for
#' @param palettename_vec a vector of palette names. check colors defined under each palette name at http://colorbrewer2.org 
#'  acceptable color PaletteNames are: YlOrRd: 9 colors PRGn: 11 colors
#'  PuOr:11 colors RdGy: 11 colors Spectral: 11 colors Grays: 9 colors PuBuGn: 9 colors RdPu: 9 colors BuPu: 9 colors YlOrBr: 9 colors Greens: 9 colors
#'  BuGn: 9 colors Accents: 8 colors GnBu: 9 colors PuRd: 9 colors Purples: 9 colors RdYlGn: 11 colors Paired: 12 colors Blues: 9 colors RdBu: 11 colors 
#'  Oranges: 9 colors RdYlBu: 11 colors PuBu: 9 colors OrRd: 9 colors Set3: 12 colors Set2: 8 colors Set1: 9 colors Reds: 9 colors PiYG: 11 colors 
#'  Dark2: 8 colors YlGn: 9 colors BrBG: 11 colors YlGnBu: 9 colors
#'  Pastel2: 8 colors Pastel1: 9 colors
#' @param colorreverseorder_vec a vector of color reverse order. true or false(default), whether to reverse the color order defined in a palette
#' @param geomtype single value (not vector) required. acceptable geomtypes are: point, multipoint, polygon, multipolygon, linestring, multilinestring, geometry
#' @param colornum_vec a vector of the number of colors to be classified. default 5.  in range of 2-15
#' @param classifier_vec a vector of classification function name. acceptable classifier function names are: Jenks(default), EqualInterval, Quantile, StandardDeviation 
#'
#' @return a list elements to be included in geolayers
#' @export
#'
#' @examples
utils.publishSP2GeoServerWithMultiStyles <- function(spobj, 
                                                    attrname_vec=c(""), 
                                                    palettename_vec=c("Reds"), 
                                                    colorreverseorder_vec=c(FALSE), 
                                                    geomtype = "Geometry", 
                                                    colornum_vec=c(5), 
                                                    classifier_vec=c("Jenks"))
  {
  # check data
  if(length(attrname_vec)==0){
    utils.debugprint("length of attrname_vec should not be 0")
    return(NULL)
  }
  # get the number of styles
  n_style = length(attrname_vec)
  
  # check the length of other parameter vectors, fill any missing elements by duplicating the value of first element
  n_vec = length(palettename_vec)
  if(n_vec < n_style){
    for (i in 1:(n_style - n_vec)){
      palettename_vec = c(palettename_vec, palettename_vec[1])
    }
  }
  
  n_vec = length(colorreverseorder_vec)
  if(n_vec < n_style){
    for (i in 1:(n_style - n_vec)){
      colorreverseorder_vec = c(colorreverseorder_vec, colorreverseorder_vec[1])
    }
  }
  
  n_vec = length(colornum_vec)
  if(n_vec < n_style){
    for (i in 1:(n_style - n_vec)){
      colornum_vec = c(colornum_vec, colornum_vec[1])
    }
  }
  
  n_vec = length(classifier_vec)
  if(n_vec < n_style){
    for (i in 1:(n_style - n_vec)){
      classifier_vec = c(classifier_vec, classifier_vec[1])
    }
  }
  
  # TODO: this is the place to add any additional column (i.e. normalised values) to spobj 
  
  # make sure spobj is projected in EPSG:4326 (WGS84) before publishing it. 
  spobj = spTransform(spobj,CRS(proj4string_epsg4326))

  
  # publish sp to geoserver
  outputWfsUrl = utils.publishSP2GeoServer(spobj)
  if(is.null(outputWfsUrl)){
    utils.debugprint("error occurs in utils.publishSP2GeoServer method, try it again")
    return(NULL)
  }
  
  # loop all styles to be created
  geolayers_elements = list()
  
  for(i in 1: n_style){
    
    attrname = attrname_vec[i]
    palettename=palettename_vec[i]
    colorreverseorder=colorreverseorder_vec[i]
    colornum=colornum_vec[i]
    classifier=classifier_vec[i]
    
    # calculte wms style 
    wmsStyleResult = fromJSON(utils.createWMSStyle(wfsurl=outputWfsUrl,
                                                   attrname=attrname,
                                                   palettename=palettename,
                                                   colorreverseorder=colorreverseorder,
                                                   geomtype=geomtype,
                                                   colornum=colornum,
                                                   classifier=classifier
    ))
    
    wmsStyleparams = list()
    if(wmsStyleResult$status==0){
      wmsStyleparams = wmsStyleResult$data
    }
    
    
    geolayers_element = list(
      layername=utils.getLayerNameFromWFSUrl(wfsurl=outputWfsUrl),
      bbox=utils.getBbox(spobj),
      wfsurl=outputWfsUrl,
      wmsstyle=wmsStyleparams
    )
    
    geolayers_elements = append(geolayers_elements, list(geolayers_element))
    
  }

  return(geolayers_elements)
  
}

#' create a featuretype in geoserver datastore, which makes the uploaded shpfile 'published'
#'
#' @param filename 
#'
#' @return empty string if success or error message
#' @export
#'
#' @examples
utils.createFeatureType <- function(filename){
  
  # publish uploaded datalayer
  h <- basicTextGatherer()
  url <- sprintf('%s/rest/workspaces/%s/datastores/%s/featuretypes.xml'
                 ,globalGSCredentials$gsRESTURL
                 ,globalGSCredentials$gsWORKSPACENAME
                 ,globalGSCredentials$gsDATASTORESNAME
                 )
  
  body <- sprintf('<featureType><enabled>true</enabled><metadata /><keywords /><metadataLinks /><attributes /><name>%s</name><title>%s</title><srs>EPSG:4326</srs><projectionPolicy>FORCE_DECLARED</projectionPolicy></featureType>'
                  ,filename
                  ,filename)
  
  # add a new featuretype by sending a POST request
  curlPerform(url = url
              ,httpheader=c(Accept="text/xml", 'Content-Type'="text/xml")
              ,username = globalGSCredentials$gsRESTUSER
              ,password = globalGSCredentials$gsRESTPW
              ,httpauth=AUTH_BASIC
              ,post=1
              ,postfields=body
              ,writefunction = h$update
              ,verbose = TRUE)
  
  #utils.debugprint(sprintf("utils.createFeatureType output: %s",h$value()))
  
  return(h$value())
}

#' create a new workspace
#'
#' @param wsname workspace name
#'
#' @return empty string if success or error message
#' @export
#'
#' @examples
utils.createWorkspace <- function(wsname){
  
  wsContentXML = xmlInternalTreeParse(utils.getWorkspace())
  
  if(is.null(wsContentXML)) return("fail to load workspaces")
    
  isExisting = (length(getNodeSet(wsContentXML, sprintf("/workspaces/workspace[name='%s']",wsname))) > 0)

  # if workspace already exists, do nothing 
  if(isExisting) return("")
  
  # otherwise, create a new workspace
  h <- basicTextGatherer()
  
  url <- sprintf('%s/rest/workspaces' ,globalGSCredentials$gsRESTURL)
  
  body <- sprintf('<workspace><name>%s</name></workspace>', wsname)
  
  curlPerform(url = url
              ,httpheader=c(Accept="text/xml", 'Content-Type'="text/xml")
              ,username = globalGSCredentials$gsRESTUSER
              ,password = globalGSCredentials$gsRESTPW
              ,httpauth=AUTH_BASIC
              ,post=1
              ,postfields=body
              ,writefunction = h$update
              ,verbose = TRUE)
  utils.debugprint(sprintf("workspace: %s created",wsname))
  return(h$value())
}


#' get all workspaces
#'
#' @param
#'
#' @return empty string if success or error message
#' @export
#'
#' @examples
utils.getWorkspace <- function(){
  
  h <- basicTextGatherer()
  
  url <- sprintf('%s/rest/workspaces' ,globalGSCredentials$gsRESTURL)
  
  curlPerform(url = url
              ,httpheader=c(Accept="text/xml", 'Content-Type'="text/xml")
              ,username = globalGSCredentials$gsRESTUSER
              ,password = globalGSCredentials$gsRESTPW
              ,httpauth=AUTH_BASIC
              ,writefunction = h$update
              ,verbose = TRUE)
  
  return(h$value())
}

#' upload a shpfile (zip) to geoserver 
#' 
#' @param filepath 
#'
#' @return empty string if success or error message
#' @export
#'
#' @examples
utils.addShp2DataStore <- function(filepath){
  
  #create workspace if it doesn't exist
  utils.createWorkspace(globalGSCredentials$gsWORKSPACENAME)
  
  #ref: https://github.com/omegahat/RCurl/issues/18
  
  h <- basicTextGatherer()
  url <- sprintf('%s/rest/workspaces/%s/datastores/%s/file.shp'
                 ,globalGSCredentials$gsRESTURL
                 ,globalGSCredentials$gsWORKSPACENAME
                 ,globalGSCredentials$gsDATASTORESNAME
                )
  
  content.type <- guessMIMEType(filepath, "application/zip")
  
  # upload shpfile by sending a PUT request
  res <- ftpUpload(what=filepath
                   ,to=url
                   ,httpheader = c('Content-Type'=content.type[[1]])
                   ,customrequest='PUT'
                   ,upload=TRUE
                   ,httpheader=c(Accept="*/*",'Content-Type'="application/zip")
                   ,username = globalGSCredentials$gsRESTUSER
                   ,password = globalGSCredentials$gsRESTPW
                   ,timeout = 36000
                   ,httpauth=AUTH_BASIC
                   ,verbose = TRUE
                   ,writefunction = h$update)
  
  #utils.debugprint(sprintf("utils.addShp2DataStore output: %s",h$value()))
  return(h$value())
}



#' returns the function that is called as the READFUNCTION callback.
#'
#' @param file 
#'
#' @return
#' @export
#'
#' @examples
utils.uploadFunctionHandler <- function(file){
  # this method is learnt from RCurl\R\upload.R in https://github.com/omegahat/RCurl
      function(size) {
        readBin(file, raw(), size)
      }
}


#
#' create wms style parameters for a given wfs datalayer
#'
#' @param wfsurl wfs url for the data layer
#' @param attrname the attribute name that the style will be created for
#' @param palettename check colors defined under each palette name at http://colorbrewer2.org 
#'  acceptable color PaletteNames are: YlOrRd: 9 colors PRGn: 11 colors
#'  PuOr:11 colors RdGy: 11 colors Spectral: 11 colors Grays: 9 colors PuBuGn: 9 colors RdPu: 9 colors BuPu: 9 colors YlOrBr: 9 colors Greens: 9 colors
#'  BuGn: 9 colors Accents: 8 colors GnBu: 9 colors PuRd: 9 colors Purples: 9 colors RdYlGn: 11 colors Paired: 12 colors Blues: 9 colors RdBu: 11 colors 
#'  Oranges: 9 colors RdYlBu: 11 colors PuBu: 9 colors OrRd: 9 colors Set3: 12 colors Set2: 8 colors Set1: 9 colors Reds: 9 colors PiYG: 11 colors 
#'  Dark2: 8 colors YlGn: 9 colors BrBG: 11 colors YlGnBu: 9 colors
#'  Pastel2: 8 colors Pastel1: 9 colors
#' 
#' @param colorreverseorder true or false(default), whether to reverse the color order defined in a palette
#' @param geomtype acceptable geomtypes are: point, multipoint, polygon, multipolygon, linestring, multilinestring, geometry
#' @param colornum the number of colors to be classified. default 5.  in range of 2-15
#' @param classifier acceptable classifier function names are: Jenks(default), EqualInterval, Quantile, StandardDeviation 
#'
#' @return
#' @export
#'
#' @examples
utils.createWMSStyle <- function(wfsurl, attrname, palettename="Reds", colorreverseorder=FALSE, geomtype="Geometry", colornum=5, classifier="Jenks"){
  
  wfsurl = sprintf("%s%s%s", wfsurl, "&propertyName=", attrname)
  
  createStyleUrl = sprintf("%s%s%s",WMSStyleCreateUrl, "?attrname=", attrname)
  createStyleUrl = sprintf("%s%s%s",createStyleUrl, "&palettename=", palettename)
  createStyleUrl = sprintf("%s%s%s",createStyleUrl, "&colorreverseorder=", colorreverseorder)
  createStyleUrl = sprintf("%s%s%s",createStyleUrl, "&geomtype=", geomtype)
  createStyleUrl = sprintf("%s%s%s",createStyleUrl, "&colornum=", colornum)
  createStyleUrl = sprintf("%s%s%s",createStyleUrl, "&classifier=", classifier)
  createStyleUrl = sprintf("%s%s%s",createStyleUrl, "&url=", URLencode(wfsurl,reserved=TRUE))
  createStyleUrl = sprintf("%s%s%s",createStyleUrl, "&gsresturl=", URLencode(globalGSCredentials$gsRESTURL,reserved=TRUE))
  createStyleUrl = sprintf("%s%s%s",createStyleUrl, "&gsusername=", globalGSCredentials$gsRESTUSER)
  createStyleUrl = sprintf("%s%s%s",createStyleUrl, "&gspassword=", globalGSCredentials$gsRESTPW)
  
  return (getURL(createStyleUrl))
}


#' find utm zone for given longitude
#'
#' @param long Longitude
#'
#' @return UTM zone code
#' @export
#'
#' @examples
utils.long2UTM <- function(long){

 (floor((long + 180)/6) %% 60) + 1

}


#' convert data.frame to json list structure
#'
#' @param df 
#'
#' @return list
#' @export
#'
#' @examples
utils.df2jsonlist <- function(df){
  
  result = list()
  for (i in 1:nrow(df)){
    rowlist = list()
    for (attrname in colnames(df)){
      
      if(is.factor(df[i, attrname])){
        rowlist[attrname] <- as.character(df[i, attrname])
      }
      else{
        rowlist[attrname] <- df[i, attrname]
      }
    }
    result[[i]] <- rowlist
  }
  return (result)
  
}


#' get layname from wfs url 
#'
#' @return 
#' @export
#'
#' @examples
utils.getLayerNameFromWFSUrl <-function(wfsurl){
  
  result = ""
  components <- strsplit(wfsurl, "&")[[1]]
  
  for (i in 1:length(components)){
    if(startsWith(components[i], "typeName")){
      
        result=strsplit(components[i], "=")[[1]][2]
        break
      }
  }
  
  return(result)
  
}


#' get boundingbox of a sp object
#'
#' @param spobj 
#'
#' @return list(minX, minY, maxX, maxY)
#' @export
#'
#' @examples
utils.getBbox <-function (spobj){
  return(list(spobj@bbox[1,1],spobj@bbox[2,1],spobj@bbox[1,2],spobj@bbox[2,2]))
  }


#' reproject(transform) a sp object to a proper utm crs
#'
#' @param spobj 
#'
#' @return a reprojected sp object
#' @export
#'
#' @examples
utils.project2UTM <-function (spobj){
  
  # get bbox of sp
  bbox = utils.getBbox(spobj)
  
  # find the centre x and y
  centreX = (bbox[[1]]+bbox[[3]]) / 2.0
  centreY = (bbox[[2]]+bbox[[4]]) / 2.0
  
  # determie the utm zone number by centre y
  utmzone = utils.long2UTM(centreX)
  
  # construct a valid proj4string_utm
  proj4string_utm = ""
  if(centreY>0) # in the Northern Hemisphere
  {
    proj4string_utm=sprintf(proj4string_utm_template,utmzone,"")
  }else
  {
    proj4string_utm=sprintf(proj4string_utm_template,utmzone,"+south ")
  }
  
  return(spTransform(spobj,CRS(proj4string_utm)))
}

#' reproject(transform) a sp object to a WGS84 (EPSG:4326)
#'
#' @param spobj 
#'
#' @return a reprojected sp object
#' @export
#'
#' @examples
utils.project2WGS84 <-function (spobj){
  
  return(spTransform(spobj,CRS(proj4string_epsg4326)))
}
