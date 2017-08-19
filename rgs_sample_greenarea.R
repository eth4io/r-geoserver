# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Copyright 2016-2017 University of Melbourne
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
# 
# sample datasets wfs url:
# (1)greenarea: http://144.6.224.184:8080/geoserver/UADI/wfs?request=GetFeature&version=1.0.0&typeName=UADI:mcc_base_prop_use&outputFormat=json&cql_filter=usecode=%27LO%27%20OR%20usecode=%27LR%27
# 
# (2)pop (Melbourne Meshblock): http://144.6.224.184:8080/geoserver/wfs?service=wfs&version=1.0.0&request=GetFeature&typeName=UADI:mb_2011_aust_pop&outputFormat=json&cql_filter=INTERSECTS(geom,%20buffer(collectGeometries(queryCollection(%27UADI:lga_2011_aust%27,%27geom%27,%27lga_code11=%27%2724600%27%27%27)),-0.00001))
#
# DevLogs:
#
# v1.0 2017-07-18
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

library(maptools) 
library(rgdal)
library(rgeos)
library(jsonlite)
library(doParallel)

# using 4 cores for parallel computing
registerDoParallel(cores=4) 

# change working directory to your own dir path where the r-geoserver.zip is unzipped to
setwd("C:\\Users\\chen1\\SourceCodeRepos\\r-geoserver")

# load utils methods for use
source("rgs_utils.R")

# calcuate green area index for city of melbourne using meshblock population
execIndicatorGreenArea <- function(){
  
  # the follow two lines are for testing
  greenarea_wfsurl = "http://144.6.224.184:8080/geoserver/UADI/wfs?request=GetFeature&version=1.0.0&typeName=UADI:mcc_base_prop_use&outputFormat=json&cql_filter=usecode=%27LO%27%20OR%20usecode=%27LR%27"
  pop_wfsurl = "http://144.6.224.184:8080/geoserver/wfs?service=wfs&version=1.0.0&request=GetFeature&typeName=UADI:mb_2011_aust_pop&outputFormat=json&cql_filter=INTERSECTS(geom,%20buffer(collectGeometries(queryCollection(%27UADI:lga_2011_aust%27,%27geom%27,%27lga_code11=%27%2724600%27%27%27)),-0.00001))"
  
  # load spatial object direct from geojson
  sp_greenarea = utils.loadGeoJSON2SP(greenarea_wfsurl)
  # check if data layer can be successfully loaded
  if(is.null(sp_greenarea)){
    utils.debugprint("fail to load data layer for greenarea")
    return(FALSE)
  }
  
  sp_pop = utils.loadGeoJSON2SP(pop_wfsurl)
  # check if data layer can be successfully loaded
  if(is.null(sp_pop)){
    utils.debugprint("fail to load data layer for population")
    return(FALSE)
  }
  
  pop_basenum = 100000 # green area index calculation base is for 100,000 persons
  
  # # # # # # # # # # # # # 
  # implement indicator logic here, such as
  # 1. spatial data operations like projection, intersection, union, 
  # 2. statistics generation
  # 3. etc.
  # # # # # # # # # # # # # 
  
  # project(transform) sp into UTM to enable area calculation
  sp_greenarea_prj = utils.project2UTM(sp_greenarea)
  
  # # # # # # # # # # # # # # # # # # # # # # 
  # calculation greenarea for each population polygon 
  # # # # # # # # # # # # # # # # # # # # # # 
  
  # project(transform) sp into UTM to enable area calculation, since the orginal population polygin in WGS84 coordinate reference system
  sp_pop_prj = utils.project2UTM(sp_pop)
  
  # add two more attributes for sp_pop_prj, one is for the actual size of greenarea, the other is for the greenarea index
  sp_pop_prj@data[,"gaarea"] = 0.0
  sp_pop_prj@data[,"idxval"] = 0.0
  
  # ==== main loop starts here ====
  # for each population polygon, find all green areas it intersects and calcuate the size of intersected area
  # two methods are implemented:
    
  # ===================================
  # method 1 : using parallel computing
  # ===================================
  
  # bulid foreach result as a matrix containing 3 columns storing values for gaarea, idxval and mb_code11 respectively
  result = foreach(i=1:nrow(sp_pop_prj), .combine = rbind, .export=c("SpatialPolygons","over","gIntersection","gArea")) %dopar% {

    # get the geometry polgyon of population, return 0 for gaarea and idxval if geometry is NULL
    if(is.null(sp_pop_prj@polygons[i])){
      out = c(0, 0)
    }else{

      geom_pop = SpatialPolygons(sp_pop_prj@polygons[i], proj4string=sp_pop_prj@proj4string)

      # accumulate the total size of intersected greenarea for the current population geometry
      intersectedGreenArea = 0.0

      # this 'over' method is much faster to find all intersected green area polygons of current pop polygon
      # temporarily save all intersected greenarea into a sub spatialdataframe
      intersectedGADF = sp_greenarea_prj[!is.na(over(sp_greenarea_prj,sp_pop_prj[i,]))[,1],]

      # if intersected with one or more greenarea polygon, calculate and accumulate the intersected area for each population meshblock
      if(nrow(intersectedGADF)>0){

        for(j in nrow(intersectedGADF):1){

          geom_greenarea = SpatialPolygons(intersectedGADF@polygons[j], proj4string=intersectedGADF@proj4string)

          # do the actual intersction process
          intsectedGeom = gIntersection(geom_pop, geom_greenarea)
          # accumulate the size of intersected greenarea
          intersectedGreenArea = intersectedGreenArea + gArea(intsectedGeom)

        }
      }

      # check population attribute, make sure it is valid
      population = sp_pop_prj@data[i,"persons"]

      if(is.null(population)||is.na(population)) population=0

      # for those polygons with 0 population, assign idxval = 0
      idx_val = 0
      if(population>0){
        idx_val = intersectedGreenArea / (population / (pop_basenum * 1.0))
      }

      out = c(intersectedGreenArea, idx_val)
    }
  }

  # assign calculated values back to sp_pop_prj@data. use as.numberic() to assure the values are numeric
  sp_pop_prj@data[,"gaarea"] = as.numeric(result[,1])
  sp_pop_prj@data[,"idxval"] = as.numeric(result[,2])
  
  # ===================================
  # method 2: using normal for loop
  # ===================================
  
  # this process takes long time to accomplish. 
  # in RStudio, use Ctrl+Shift+C to uncomment/comment it for testing
  
  
  # for(i in nrow(sp_pop_prj):1){
  # 
  #   utils.debugprint(sprintf("processing [%i/%i]", i, nrow(sp_pop_prj)))
  # 
  #   # get the geometry polgyon of population, skip if it is NULL
  #   if(is.null(sp_pop_prj@polygons[i])){
  #     next
  #   }
  #   geom_pop = SpatialPolygons(sp_pop_prj@polygons[i], proj4string=sp_pop_prj@proj4string)
  # 
  #   # accumulate the total size of intersected greenarea for the current population geometry
  #   intersectedGreenArea = 0.0
  # 
  #   # this 'over' method is much faster to find all intersected green area polygons of current pop polygon
  #   # temporarily save all intersected greenarea into a sub spatialdataframe
  #   intersectedGADF = sp_greenarea_prj[!is.na(over(sp_greenarea_prj,sp_pop_prj[i,]))[,1],]
  # 
  #   # if intersected with one or more greenarea polygon, calculate and accumulate the intersected area for each population meshblock
  #   if(nrow(intersectedGADF)>0){
  # 
  #     for(j in nrow(intersectedGADF):1){
  # 
  #       geom_greenarea = SpatialPolygons(intersectedGADF@polygons[j], proj4string=intersectedGADF@proj4string)
  # 
  #       # do the actual intersction process
  #       intsectedGeom = gIntersection(geom_pop, geom_greenarea)
  #       # accumulate the size of intersected greenarea
  #       intersectedGreenArea = intersectedGreenArea + gArea(intsectedGeom)
  # 
  #     }
  #   }
  # 
  #   # check population attribute, make sure it is valid
  #   population = sp_pop_prj@data[i,"persons"]
  # 
  #   if(is.null(population)||is.na(population)) population=0
  # 
  #   # for those polygons with 0 population, assign idxval = 0
  #   if(population>0){
  #     sp_pop_prj@data[i,"idxval"] = intersectedGreenArea / (population / (pop_basenum * 1.0))
  #   }
  #   # assgin intersectedGreenArea to gaarea attribute
  #   sp_pop_prj@data[i,"gaarea"] = intersectedGreenArea
  # 
  # }
  
  
  # ==== main loop ends here ====
  
  # this example shows how to publish a geolayer by creating multiple wms styles on various attributes of the same data layer. 
  # the data layer will be only published one time, with various wms styles generated for selected attributes 
  publishedinfo = utils.publishSP2GeoServerWithMultiStyles(spobj=sp_pop_prj, 
                                                               attrname_vec=c("gaarea","idxval"),
                                                               palettename_vec=c("Reds","Blues"), 
                                                               colorreverseorder_vec=c(FALSE,FALSE), 
                                                               geomtype = "Geometry", 
                                                               colornum_vec=c(6,8), 
                                                               classifier_vec=c("Jenks","Jenks")
                                                               )
  
  if(is.null(publishedinfo) || length(publishedinfo)==0){
    utils.debugprint("fail to save data to geoserver")
    return(FALSE)
  }
  
  
  # print the outputs in json format
  utils.debugprint(sprintf("outputs: %s", toJSON(publishedinfo, auto_unbox=TRUE)))
  

  return(TRUE)
}

execIndicatorGreenArea()

accessWFSWithAuth <- function(){
  
  # this wfs is protected by username/password
  wfsurl = "http://115.146.92.213:8080/geoserver/wfs?service=wfs&version=1.0.0&request=GetFeature&typeName=topp:tasmania_cities&outputFormat=json"
  username = "wfs_reader"
  password = "wfs_password"
  
  sp = utils.loadGeoJSON2SPWithAuth(wfsurl, username, password)
  #sp_failed  = utils.loadGeoJSON2SP(wfsurl)
  
  # load dataframe
  wfsurl2 = "http://115.146.92.213:8080/geoserver/wfs?service=wfs&version=1.0.0&request=GetFeature&propertyName=CITY_NAME,ADMIN_NAME&typeName=topp:tasmania_cities&outputFormat=json"

  df = utils.loadGeoJSON2DFWithAuth(wfsurl2, username, password)
  #df_failed  = utils.loadGeoJSON2DF(wfsurl2)
  
}


