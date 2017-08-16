# README #

This repo contains a utility  creating UADI indicator tools using R.

"rgs_sample_greenarea.R" contains a sample source code of creating greenarea indicator using R and GeoServer.

"rgs_utils.R" contains all utility functions to interact GeoServer with R as simple as possible. 


All codes are tested in R 3.3 environment on windows and linux.


### How do I get set up? ###

A couple of softwares and tools are required to create UADI indicators usring R:

(1) Install R (version 3.3 or above) from the official website

https://cran.r-project.org/


(2) RStudio is recommended for developing R scripts. The "RStudio Desktop Open Source License" version just works great.

https://www.rstudio.com/products/rstudio/download/


(3) Install Rtools (for Windows developers only) from here 

https://cran.r-project.org/bin/windows/Rtools/ 

since zip.exe in the Rtools packages will be used in the "uadi_utils.R", please make sure "YOUR_Rtools_DIRECTORY/bin" is included in the system "Path" environment variable.

(4) When testing your indicator script, RScript.exe which locates in "YOUR_R_DIRECTORY\bin" will be used, make sure this directory is included in the system "Path" environment variable.

(5) Download the source code from gitbub

(6) install R packages on your dev environment. In RStudio, go to "Tools"->"Install Packages", make sure the "Install dependencies" option is checked, then install the following packages one by one:

'sp', maptools', 'rgdal', 'rgeos', 'jsonlite, 'geojsonio', 'RCurl', 'uuid'

you might need to close and reopen RStudio to use the installed packages



### Who do I talk to? ###

If you have any questions or suggestions, please contact:

Dr Yiqun Chen

Centre for Disaster Management & Public Safety

Centre for Spatial Data Infrastructures & Land Administration

The University of Melbourne

E: yiqun.c@unimelb.edu.au