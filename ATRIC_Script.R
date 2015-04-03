
########################################################################
########################################################################
##                                                                    ##
##                               ATRIC                                ##
##                                                                    ##
##         Automated Accumulation Threshold selection                 ##
##                               and                                  ##
##                   RIparian Corridor delineation                    ##
##                                                                    ##
########################################################################
########################################################################

## Required R packages:
library(raster)
library(rgdal)
library(maptools)
library(rgeos)
library(spgrass6)
library(plyr)
library(tcltk)
## For more information about the packages and how to install them,
# please visit: http://www.r-project.org/

## Initiating GRASS GIS from the installed location
## Linux
location <- initGRASS("/usr/lib/grass-6.4.3", override=TRUE)

## Mac OS X
#location <- initGRASS("/Applications/GRASS-6.4.app/Contents/MacOS",home=tempdir(),override=TRUE)

## MS Windows 
# Note for MS Windows users:
# Windows directories come with backslashes "\" that are not identified by R.
# Please change them to regular slashes "/" or to double backslashes "\\".
# Here, for example, "C:\Program Files (x86)\GRASS GIS 6.4.3" has been changed to
# "C:/Program Files (x86)/GRASS GIS 6.4.3"
# location <- initGRASS("C:/Program Files (x86)/GRASS GIS 6.4.3", override=TRUE)

## For details on spgrass6 and GRASS GIS initiation instructions please visit
# "http://cran.r-project.org/web/packages/spgrass6/index.html"
# and
# "http://grasswiki.osgeo.org/wiki/GRASS_Help#First_Day_Documentation"

## Set the path to the working directory containing the mapped stream network (MSN),
# digital elevation model (DEM), mask and stream sampling points (SSP).
# For MS Windows OS, see above note on directory separators.
# Linux
path <- "/media/storage/projects/Uni_Landau/Data/"
# Mac OS X
# path <- "/Users/<username>/Desktop/Stream Derivation/Thueringen"
# MS Windows 
# path <- "C:/<username>/Stream_Extraction/Thueringen"

setwd(path)

## Load the required digital elevation model (DEM), the first input for automated accumulation threshold (AT)
# selection into R. Note that the DEM raster file (.tif) is named as DEM.
DEM<-raster("DEM.tif")

## The DEM has the spatial reference information. We will create projection information file for the
# GRASS GIS temporary directory
currmapset <- execGRASS("g.gisenv", parameters=list(get="MAPSET"), intern=TRUE)
execGRASS("g.mapset", parameters=list(mapset="PERMANENT"))
execGRASS("g.proj", flags=c("c"), parameters=list(georef="DEM.tif"))
execGRASS("g.mapset", parameters=list(mapset=currmapset))

## Load and store the DEM into the temporary directory for the GRASS GIS usage.
execGRASS("r.in.gdal", flags=c("o", "overwrite"), parameters=list(input=paste(path,"DEM.tif", sep="/"), output="DEM"))

## Load the required mapped stream network (MSN) shapefile, the second input for automated AT
# selection into R and then store it into the temporary directory for GRASS usage. 
# Note that the MSN file (.shp) is named as SSN (matching the earlier version of the paper, will be renamed soon).
# For R usage the spatial reference system ('NA' at the beginning ) of the MSN will be set
# to the spatial reference system of the DEM. The region of analysis will be set to the 
# spatial extent of the DEM in GRASS
MSN<-readShapeLines("SSN")
execGRASS("v.in.ogr", flags=c("o", "overwrite"), parameters=list(dsn=path, output="MSN", layer="SSN", type="line"))
proj4string(MSN) <- proj4string(DEM)
execGRASS("g.region", parameters=list(rast="DEM")) 

## A DEM, like any raster file comes with a rectengular spatial extent which doesn't necessarily
# reflect the actual spatial extent of the MSN (line or polygon). To avoid artefacts during the
# stream extraction e.g. in coastal areas, the DEM and the extracted streams will be masked with
# the spatial extent of the MSN. For this purpose, a polygon shapefile of the geographic region 
# of interest that contains the MSN will be loaded into R (third input for automated AT selection).
# Alternatively, one could create a 'Convex Hull' that contains the MSN (for example if a shapefile
# of the study region was unavailable) and import it into R and use it to mask the DEM and the
# accumulation map derived stream network. Note that the provided mask file (.shp) is named as Mask
Mask <- readShapePoly("Mask")
proj4string(Mask) <- proj4string(DEM)
execGRASS("v.in.ogr", flags=c("o", "overwrite"), parameters=list(dsn=path, output="Mask", layer="Mask"))

## Alternative: Create and import a covex hull
# execGRASS("v.hull", flags="overwrite", parameters=list(input="MSN", output="Mask"))
# Mask <- readVECT6("Mask")

## Create an Accumulation and a Drainage map in GRASS and import it to R.
# To avoid mismatch between the DEM and Accumulation map, the accumulation map will be 
# generated from the same DEM that will be used for stream extraction
execGRASS("r.watershed", flags=c("overwrite"), parameters=list(elevation="DEM", accumulation="Accumulation", drainage="Drainage"))
Accumulation<-raster(readRAST6("Accumulation"))

## The MSN will be converted to a raster MSN 
execGRASS("v.to.rast", flags="overwrite", parameters=list(input="MSN", output="MSN_raster", use="val", rows=as.integer(nrow(Accumulation))))
MSN_raster <- raster(readRAST6("MSN_raster"))

## To check whether the spatial extent of the DEM, known stream and mask
# match, one can plot them on top of each other
plot(DEM)
plot(Mask, col=rgb(red=0, green=50, blue=0, alpha=50, max=255), add=TRUE)
plot(MSN, col="blue", add=TRUE)


########################################################################
##         Automated accumulation threshold selection method          ##
########################################################################

### Steps:

## At the beginning of this method, the nodes in the MSN will be identified.
# A node is a point on a stream network where a stream starts or ends,
# or meets or crisscrosses with another stream
# (see http://www.ing.unitn.it/~grass/docs/tutorial_64_en/htdocs/esercitazione/network_analysis/
# for details). For the purpose, the polyline shapefile of the stream  network will be
# converted to a network shapefile and the nodes will be identified with
# the 'Node' operation by the 'v.net' module in GRASS GIS
execGRASS("v.net", flags=c("c", "overwrite"), parameters=list(input="MSN", output="MSN_Network", operation="nodes"))

## The created nodes will be extracted using the 'v.extract' module in GRASS GIS.
# The node information are stored in the second layer, the number of nodes created
# can be extracted with the 'vInfo' function. The extracted nodes will be saved as
# a new shapefile to store all the corresponding information in a database and then
# the new node shapefile will be loaded into R. Note that the node information will
# be lost if a direct import is performed using 'readVECT6'
execGRASS("v.extract", flags="overwrite", parameters=list(list=paste("1", vInfo("MSN_Network")[[1]], sep="-"), input="MSN_Network", output="MSN_Nodes",  type="point", layer=as.integer(2)))
execGRASS("v.out.ogr", parameters=list(input="MSN_Nodes", type="point", dsn=path, olayer="MSN_Nodes"))
MSN_Nodes<-readShapePoints("MSN_Nodes")

## The next target is to obtain the DEM stream sources. For this purpose, the number
# of stream segments connected to each of the nodes will be counted and stored.
# The number of stream segments connected to each of the nodes will be identified 
# within the 'v.net' module using 'nreport' operation.
segment_info <- execGRASS("v.net", flags=c("c", "overwrite"), parameters=list(input="MSN_Network", operation="nreport"), intern = TRUE)
number_segments <- sapply(strsplit(segment_info, split=" ", fixed=TRUE), function(x) length(strsplit(x[[2]], split=",", fixed=TRUE)[[1]]))

## The nodes that have only one stream segment connected are either sources or outlets
# (not known at this stage), therefore they will be called 'mapped single segment nodes'.
# The single segment nodes will be extracted from the nodes shapefile and stored
MSN_Nodes@data <- data.frame(MSN_Nodes@data, number_segments)
MSN_single_segment_nodes <- subset(MSN_Nodes, MSN_Nodes@data$number_segments==1)

## Nodes in a stream network are connected by parts of the network (connections).
# At this stage, connections between the identified single segment nodes will be
# extracted. This allows for grouping the single segment nodes by connections and
# thus every connection will have one (or multiple) source and one outlet.
# This will be done by 'v.net.components' module in GRASS GIS where the 'weakly'
# connected streams will be identified
# (see http://grass.osgeo.org/grass64/manuals/v.net.components.html for details).
# As an outcome, all the streams that connect a (or multiple) source and an outlet
# will be given the same component ids. The output will be imported to R
execGRASS("v.net.components", flags="overwrite", parameters=list(input="MSN_Network", output="MSN_connectivity", method="weak"))
MSN_connectivity <- readVECT6("MSN_connectivity")

## In the next step, previously identified single segment nodes will be evaluated
# on the basis of the connections they belong to. Therefore, with the GRASS GIS
# module 'v.distance' the ids of the stream they belong to will be extracted.
# The extracted ids will be stored in the single segment nodes database 
MSN_single_segment_nodes@data$lcat <- rep(0, nrow(MSN_single_segment_nodes))
writeVECT6(MSN_single_segment_nodes, "MSN_single_segment_nodes", v.in.ogr_flags=c("o", "overwrite"))
execGRASS("v.distance", flags="overwrite", parameters=list(from="MSN_single_segment_nodes", to="MSN_connectivity", upload="cat", column="lcat"))
MSN_single_segment_nodes_stream <- readVECT6("MSN_single_segment_nodes")

## For Windows OS users the following need to be done instead of L 171
# execGRASS("v.out.ogr", parameters=list(input="MSN_single_segment_nodes", type="point", dsn=path, olayer="MSN_single_segment_nodes"))
# MSN_single_segment_nodes_stream <- readShapePoints("MSN_single_segment_nodes")

MSN_single_segment_nodes_stream@data <- merge(MSN_single_segment_nodes_stream@data, MSN_connectivity@data, by.x="lcat", by.y="cat")

## A DEM-derived stream newtork (DSN) usually shows lateral displacements (d) from
# the MSN. This implies that the cell with highest accumulations within the d from
# the mapped single segment nodes indicate the corresponding DEM single segment nodes
# (please see the paper for details). Here, with the purpose of automatically selecting
# the d between MSN and DSN, we will identify the DEM single segment nodes for varying ds
# and will group them by connection as previously shown. Thereafter the accumulation values
# of the cells at the grouped DEM single segment nodes will be identified and compared for
# every d.The accumulation values at the sources are always lower than at the outlets or 
# vice versa. Thus we will identify the source in each connection for each d will be identified
# (in case of multiple sources the source with the lowest accumulation will be selected as
# this corresponds to accumulation threshold (AT)) and then for each d the mean of the
# sources accumulations will be selected as the AT for that d. Then we will extract a stream
# network for each d using the resulted AT from DEM. As the number of DEM extracted stream cells (nDS)
# is an inverse function of AT and thus an inverse function of d (please see the paper for details)
# the d that derives an AT for extracting optimal nDS will be selected. For this purpose, nDS
# extracted by different ds will be compared to the number of mapped stream cells (nSS).
# The selected d will be the one for which the extracted AT will result in the nDS that is 
# approx. +5% of the nSS (please see the paper for the rationale). The whole algorithm is 
# implemented on a 'Repeat loop' in R with 'ifelse' arguments

## The initial d will be set two-fold the resolution of the DEM
initd <- 2*res(DEM)[[1]]
i <- initd

# Empty lists will be created where the nDS and corresponding ATs will be stored for different ds
nDS <- list()
AT_mean <- list()
AT_percentile <- list()

# The first stream network will be derived using the initial d and then the nDS will be computed
# for it and stored. Then in the repeat loop the d will first be consecutively increased by 50m 
# and for every increase the new nDS will be compared to the nSS. This process will continue until
# the nDS decreases. This might happen at the very first increase of 50m in the d or at the 
# initial d. Consequently, at the next steps the d will be decreased by 10m with the target of
# increasing the nDS. This process will continue until this yields the approx. +5% nMS. 
# The corresponding AT values will be stored for further optimization. Note that this procedure
# could be expanded to the nearest 5 or 2 but we regard this level of precision as sufficient.
# In the mean time, if the nDS is obtained as the same of the nSS (highly unlikely),
# the training will also be stopped 
repeat{
  Accumulation_values <- extract(Accumulation, MSN_single_segment_nodes_stream, method="simple", buffer=i, small=TRUE)
  Accumulation_values <- sapply(Accumulation_values, FUN="max", simplify=TRUE)
  MSN_single_segment_nodes_stream@data$Accumulation_values <- Accumulation_values
  Sources_Accumulation_Values <- aggregate(MSN_single_segment_nodes_stream@data, by=list(MSN_single_segment_nodes_stream@data$comp), FUN="min")
  AT_mean[[i]] <- mean(na.omit(Sources_Accumulation_Values$Accumulation_values))
  AT_percentile[[i]] <- quantile(na.omit(Sources_Accumulation_Values$Accumulation_values), c(0.05, 0.5, 0.95))
  execGRASS("r.stream.extract", flags="overwrite", parameters=list(elevation="DEM", accumulation="Accumulation",  threshold=round(AT_mean[[i]], digit=0), mexp=0.01, stream_rast="Stream_Comparison"))
  Stream_Comparison <- raster(readRAST6("Stream_Comparison"))
  Stream_Comparison <- mask(Stream_Comparison, Mask)
  nDS[[i]] <- length(Which(Stream_Comparison>0, cells=TRUE))
  if(nDS[[i]]<length(Which(MSN_raster==1, cells=TRUE))){
    repeat{
      i <- i-initd/5
      Accumulation_values <- extract(Accumulation, MSN_single_segment_nodes_stream, method="simple", buffer=i, small=TRUE)
      Accumulation_values <- sapply(Accumulation_values, FUN="max", simplify=TRUE)
      MSN_single_segment_nodes_stream@data$Accumulation_values <- Accumulation_values
      Sources_Accumulation_Values <- aggregate(MSN_single_segment_nodes_stream@data, by=list(MSN_single_segment_nodes_stream@data$comp), FUN="min")
      AT_mean[[i]] <- mean(na.omit(Sources_Accumulation_Values$Accumulation_values))
      AT_percentile[[i]] <- quantile(na.omit(Sources_Accumulation_Values$Accumulation_values), c(0.05, 0.5, 0.95))
      execGRASS("r.stream.extract", flags="overwrite", parameters=list(elevation="DEM", accumulation="Accumulation",  threshold=round(AT_mean[[i]], digit=0), mexp=0.01, stream_rast="Stream_Comparison"))
      Stream_Comparison <- raster(readRAST6("Stream_Comparison"))
      Stream_Comparison <- mask(Stream_Comparison, Mask)
      nDS[[i]] <- length(Which(Stream_Comparison>0, cells=TRUE))
      if(nDS[[i]]>=length(Which(MSN_raster==1, cells=TRUE)) & nDS[[i]]-length(Which(MSN_raster==1, cells=TRUE))<=0.05*length(Which(MSN_raster==1, cells=TRUE))){
        repeat{
          i <- i-initd/5
          Accumulation_values <- extract(Accumulation, MSN_single_segment_nodes_stream, method="simple", buffer=i, small=TRUE)
          Accumulation_values <- sapply(Accumulation_values, FUN="max", simplify=TRUE)
          MSN_single_segment_nodes_stream@data$Accumulation_values <- Accumulation_values
          Sources_Accumulation_Values <- aggregate(MSN_single_segment_nodes_stream@data, by=list(MSN_single_segment_nodes_stream@data$comp), FUN="min")
          AT_mean[[i]] <- mean(na.omit(Sources_Accumulation_Values$Accumulation_values))
          AT_percentile[[i]] <- quantile(na.omit(Sources_Accumulation_Values$Accumulation_values), c(0.05, 0.5, 0.95))
          execGRASS("r.stream.extract", flags="overwrite", parameters=list(elevation="DEM", accumulation="Accumulation",  threshold=round(AT_mean[[i]], digit=0), mexp=0.01, stream_rast="Stream_Comparison"))
          Stream_Comparison <- raster(readRAST6("Stream_Comparison"))
          Stream_Comparison <- mask(Stream_Comparison, Mask)
          nDS[[i]] <- length(Which(Stream_Comparison>0, cells=TRUE))
          if (nDS[[i]]>nDS[[i+initd/5]]){
          d <- i+initd/5
          print(d)
          print(AT_mean[[i+initd/5]])
          print(AT_percentile[[i+initd/5]])
          break()
        }
      }
        break()
    }
      else{
        if(nDS[[i]]>=length(Which(MSN_raster==1, cells=TRUE))){
        repeat{
          i <- i-initd/5
          Accumulation_values <- extract(Accumulation, MSN_single_segment_nodes_stream, method="simple", buffer=i, small=TRUE)
          Accumulation_values <- sapply(Accumulation_values, FUN="max", simplify=TRUE)
          MSN_single_segment_nodes_stream@data$Accumulation_values <- Accumulation_values
          Sources_Accumulation_Values <- aggregate(MSN_single_segment_nodes_stream@data, by=list(MSN_single_segment_nodes_stream@data$comp), FUN="min")
          AT_mean[[i]] <- mean(na.omit(Sources_Accumulation_Values$Accumulation_values))
          AT_percentile[[i]] <- quantile(na.omit(Sources_Accumulation_Values$Accumulation_values), c(0.05, 0.5, 0.95))
          execGRASS("r.stream.extract", flags="overwrite", parameters=list(elevation="DEM", accumulation="Accumulation",  threshold=round(AT_mean[[i]], digit=0), mexp=0.01, stream_rast="Stream_Comparison"))
          Stream_Comparison <- raster(readRAST6("Stream_Comparison"))
          Stream_Comparison <- mask(Stream_Comparison, Mask)
          nDS[[i]] <- length(Which(Stream_Comparison>0, cells=TRUE))
          if (nDS[[i]]>nDS[[i+initd/5]]){
            d <- i+initd/5
            print(d)
            print(AT_mean[[i+initd/5]])
            print(AT_percentile[[i+initd/5]])
            break()
          }
        }
        break()
      }
    }     
  }
    break()
}
  i <- i+initd
}

## The raster MSN will be buffered using the selected d. This enables evaluation of the match
# between the DSN and the MSN
execGRASS("r.buffer", flags="overwrite", parameters=list(input="MSN_raster", output="MSN_buffer",  distances=d))
MSN_buffer <- raster(readRAST6("MSN_buffer"))

## In this step, the number of overlapped DEM stream cells with the mapped stream cells will
# be optimized. Different stream networks will be extracted from the DEM using the accumulation
# map and different ATs. The first AT value will be the one obtained for the selected d.
# The other ATs will be chosen as consecutive increases or decreases to this value. For every AT a 
# raster stream network will be created and then will be overlaid with the buffered MSN.
# The percentage of the overlapped (nDS(Ov)) (with the buffer zone of the MSN) cells of the DSN
# (nDS(Ov)/nDS) will be computed (please see the paper for details regarding rationale).
# Similar to before, the whole algorithm is implemented on a 'Repeat loop' in R with 'ifelse' arguments

## Set the previously obtained AT (ATd) value as the premilinary AT value
ATd <- round(AT_mean[[i+initd/5]], digit=0)
i <- ATd

## An empty list will be created where the nDS(Ov)/nDS (perc_nDS_Ov to avoid syntax complexity)
# will be stored for different thresholds
perc_nDS_Ov <- list()

## The first stream network will be extracted outside of the loop and then the nDS(Ov)/nDS will
# be computed for it and stored
execGRASS("r.stream.extract", flags="overwrite", parameters=list(elevation="DEM", accumulation="Accumulation",  threshold=i, mexp=0.01, stream_rast="Stream_Comparison"))
Stream_Comparison <- raster(readRAST6("Stream_Comparison"))
Stream_Comparison <- mask(Stream_Comparison, Mask)
Comparison_Stream <- MSN_buffer+Stream_Comparison
perc_nDS_Ov[[i]] <- length(Which(Comparison_Stream>=0, cells=TRUE))/length(Which(Stream_Comparison>=0, cells=TRUE))*100

## Then in the repeat loop the AT will first be consecutively increased by 1000 unit (as 1000<ATd<10000 unit,
# see the paper for details) and for every increase the nDS(Ov)/nDS will be compared to the preceding
# (see paper for discussion). This process will continue until the derived nDS satisfies Eq. 2, i.e.
# derived max. +5% nMS (please see the paper for rationale).This might happen at the very first increase of
# 1000 units in the AT. Consequently, at the next steps the AT will be decreased by 1000 units with the target
# of increasing nDS(Ov)/nDS.

repeat{
  i <- i+10^(nchar(ATd)-1)
  execGRASS("r.stream.extract", flags="overwrite", parameters=list(elevation="DEM", accumulation="Accumulation",  threshold=i, mexp=0.01, stream_rast="Stream_Comparison"))
  Stream_Comparison <- raster(readRAST6("Stream_Comparison"))
  Stream_Comparison <- mask(Stream_Comparison, Mask)
  Comparison_Stream <- MSN_buffer+Stream_Comparison
  perc_nDS_Ov[[i]] <- length(Which(Comparison_Stream>=0, cells=TRUE))/length(Which(Stream_Comparison>=0, cells=TRUE))*100
  if(perc_nDS_Ov[[i]]<perc_nDS_Ov[[i-10^(nchar(ATd)-1)]] | length(Which(Stream_Comparison>0, cells=TRUE))<length(Which(MSN_raster==1, cells=TRUE))){
    repeat{
      i <- i-10^(nchar(ATd)-1)
      execGRASS("r.stream.extract", flags="overwrite", parameters=list(elevation="DEM", accumulation="Accumulation",  threshold=i, mexp=0.01, stream_rast="Stream_Comparison"))
      Stream_Comparison <- raster(readRAST6("Stream_Comparison"))
      Stream_Comparison <- mask(Stream_Comparison, Mask)
      Comparison_Stream <- MSN_buffer+Stream_Comparison
      perc_nDS_Ov[[i]] <- length(Which(Comparison_Stream>=0, cells=TRUE))/length(Which(Stream_Comparison>=0, cells=TRUE))*100
      ## As soon as the percentage decreases compared to the previous step and the nDS becomes lower than the nMS,
      # the upper boundary value of the AT is reached. To determine the optimal AT more precisely, the 1000 unit
      # will now be reduced to 100 unit (as 1000<AT<10000 unit, see the paper for the rationale). Thus in the next
      # part of the loop the AT will be increased by 100 units from the current stage and again the nDS(Ov)/nDS
      # between the current and preceding stage will be compared. This process will continue until this yields
      # a negative value (percentage decreases with a higher AT) or a decrease in the nDS compared to the nMS.
      # This is regarded as the optimal AT to start a stream network on DEM. Note that this procedure could be
      # expanded to the nearest 10 but we regard this level of precision as sufficient
      
      if(perc_nDS_Ov[[i]]<perc_nDS_Ov[[i+10^(nchar(ATd)-1)]] & length(Which(Stream_Comparison>0, cells=TRUE))>length(Which(MSN_raster==1, cells=TRUE))){
        repeat{
          i <- i+10^(nchar(ATd)-2)
          execGRASS("r.stream.extract", flags="overwrite", parameters=list(elevation="DEM", accumulation="Accumulation",  threshold=i, mexp=0.01, stream_rast="Stream_Comparison"))
          Stream_Comparison <- raster(readRAST6("Stream_Comparison"))
          Stream_Comparison <- mask(Stream_Comparison, Mask)
          Comparison_Stream <- MSN_buffer+Stream_Comparison
          perc_nDS_Ov[[i]] <- length(Which(Comparison_Stream>=0, cells=TRUE))/length(Which(Stream_Comparison>=0, cells=TRUE))*100
          if(perc_nDS_Ov[[i]]<perc_nDS_Ov[[i-10^(nchar(ATd)-2)]] | length(Which(Stream_Comparison>0, cells=TRUE))<length(Which(MSN_raster==1, cells=TRUE))){
            AT <- i-10^(nchar(ATd)-2)
            print(AT)
            break()
          }
        }
        break()
      }
    }
    break()
  }
}
##(length(Which(Stream_Comparison>0, cells=TRUE))-length(Which(MSN_raster==1, cells=TRUE)))/length(Which(MSN_raster==1, cells=TRUE))*100
## Finally, the stream network will be extracted from the DEM using the obtained AT value as before
execGRASS("r.stream.extract", flags="overwrite", parameters=list(elevation="DEM", accumulation="Accumulation", threshold=round(AT, digits=0), mexp=0.01, stream_vect="Stream_Comparable_known"))
execGRASS("v.overlay", flags=c("t","overwrite"), parameters=list(ainput="Stream_Comparable_known", atype="line", binput="Mask", btype="area", output="DSN", operator="and"))

# The output can be stored, plotted and compared. For better visual comparison, we suggest to display
# the DSN and MSN with GUI GIS softwares. The files that would not be used in the next algorithm will
# be removed from the workspace
execGRASS("v.out.ogr", parameters=list(input="DSN", type="line", dsn=path, olayer="DSN"))
DSN <- readShapeLines("DSN")
rm(Accumulation)
rm(MSN_raster)
rm(MSN_buffer)
rm(Stream_Comparison)
rm(Comparison_Stream)
rm(MSN_Nodes)
rm(MSN_single_segment_nodes)
rm(MSN_connectivity)
plot(DSN)
plot(MSN, col="blue", add=TRUE)
legend(3699000, 5699000, c("MSN","DSN"), cex=0.8, col=c("blue","black"), lty=1)



########################################################################
#               Upstream riparian corridor delineation                 #
#                  for given stream sampling points                    #
########################################################################


### Steps:

## Load the stream sampling points (SSP) from governmental biomonitoring
# shapefile into R
SSP<-readShapePoints("SSP")
proj4string(SSP) <- proj4string(DEM)

## The biomonitoring data were collected at the MSN and therefore ideally
# the SSP should exactly fall on the MSN when overlaid. But while compared
# many SSP showed substantial lateral displacements from MSN (see paper).
# Obviously, this an artefact of the different human resources and equipments
# used by the responsible authorities for collecting and digitizing the 
# biomonitoring sites and mapped streams. But the name of the streams
# where the samples were collected were input in the database of the
# biomoitoring, which allows to evaluate whether a particular sample site
# corresponds to the nearest stream if it is not exactly on that stream.
# Therefore each of the SSP will be snapped to the stream that matches by name
# where the biomonitoring data was sampled. The stream names in the SSP and
# the MSN are stored under "GEWAESSERN" and "HAUPTNAME" ids, respectively.
# First, it will be checked whether or not all the stream names of the
# SSP matches to the stream names of the MSN

SSP@data$GEWAESSERN <- as.character(SSP@data$GEWAESSERN)
MSN@data$RS_NAME <- as.character(MSN@data$RS_NAME)
no_SSP<- NULL
nr_match_MSN <-NULL
for (i in 1:nrow(SSP)){
  no_SSP<- c(no_SSP, i)
  nr_match_MSN <- c(nr_match_MSN, length(which(MSN@data$RS_NAME==SSP[i,]@data$GEWAESSERN)))
}
match_MSN_nrs <- data.frame(no_SSP, nr_match_MSN)
non_matched_SSP_perc <- nrow(subset(match_MSN_nrs, nr_match_MSN==0))/nrow(match_MSN_nrs)*100
print(non_matched_SSP_perc)

## As observed, 26% of the SSP has not matched by the stream names.
# While exploring the non-matched SSP, it was observed that the most of the
# mismatches occurred due to writing the stream names differently (spelling)
# in the two different databases, different abbreviations and versions
# of the names and the presence of the special characters in the name. These
# problems can be solved by hand and then almost all the SSP stream names will
# be matched to the ones of MSN. To be simplistic, in this algorithm only the
# matched SSP will be used for the further analyses and therefore be snapped to
# the MSN. The function 'snapPointsToLines' and the 'spRbind' functions in the
# 'maptools' package of R will be applied for the purpose
matched_SSP <- subset(match_MSN_nrs, nr_match_MSN!=0)
Snapped_SSP_MSN <- snapPointsToLines(SSP[1,], MSN[which(MSN@data$RS_NAME==SSP[1,]@data$GEWAESSERN),])
for (i in matched_SSP$no_SSP[-1]){
  Snapped_SSP_MSN <- spRbind(Snapped_SSP_MSN, snapPointsToLines(SSP[i,], MSN[which(MSN@data$RS_NAME==SSP[i,]@data$GEWAESSERN),]))
}

## Now all the SSP are exactly lying on the streams where they were sampled.
# But this does not allow for any hydrological analyis on DEM since the MSN
# from the government authority does not represent the stream that could be
# derived from the DEM. Mostly, the DEM derived streams are slightly shifted
# from the equivalent MSN. For example, while attempted to identify the upstream UC
# for a particular SSP using 'r.water.outlet', the results were a few scattered cell on the
# DEM. But when the same SSP is snapped to the DSN that is comparable to the MSN, it resulted
# in a reliable UC. Please consult the paper and the graphics for further details. Therefore
# the DSN that is comparable to the MSN can be applied to snap the SSP and furthermore
# for upstream riparian corridor (URC) delineation. The length of the URC is set to 10km upstream
# from the given SSP and the width is 100m in each side parallel to the stream 

## First the SSP on the MSN will be snapped to the comparable DSN. 
# The function 'snapPointsToLines' in the 'maptools' package of R will be applied for the purpose.
# The GRASS GIS function 'v.distance' can be used to reduce execution time
Snapped_SSP_DSN<- snapPointsToLines(Snapped_SSP_MSN, DSN)

## Alternative with 'v.distance' function in GRASS
#Snapped_SSP_MSN@data$x <- rep(0, nrow(Snapped_SSP_MSN)) 
#Snapped_SSP_MSN@data$y <- rep(0, nrow(Snapped_SSP_MSN)) 
#writeVECT6(Snapped_SSP_MSN, "Snapped_SSP_MSN", v.in.ogr_flags=c("o", "overwrite"))
#execGRASS("v.distance", parameters=list(from="Snapped_SSP_MSN", to="DSN", upload=c("cat", "to_x", "to_y"), column=c("near_stream_cat", "x", "y")))
#Snapped_SSP_DSN<- readVECT6("Snapped_SSP_MSN")
#Snapped_SSP_DSN@coords[,1] <- Snapped_SSP_DSN@data$x
#Snapped_SSP_DSN@coords[,2] <- Snapped_SSP_DSN@data$y

## Then, two empty objects will be created to store the computed
# upstream UC and riparian corridor polygons
UC <- NULL
URC <- NULL

## A progressbar will be created to check the status of the process
ProgressBar <- tkProgressBar(title = "progress bar", min = 0,max = nrow(Snapped_SSP_DSN), width = 300)

## A 'for' loop will be applied for the purpose.
## Note for users: URC delineation time depends on the number of SSP
# and computer configuration. To save memory space we will remove the
# processed files at each step
for(i in 1:nrow(Snapped_SSP_DSN)) {
  ## First the UCs for each of the SSP will be extracted using the 'r.water.outlet' 
  # module in GRASS GIS. For the purpose, the snapped SSP will be exported to GRASS,
  # the raster UCs will be extracted, converted to polygons and imported to R 
  writeVECT6(Snapped_SSP_DSN[i,], "Snapped_SSP_DSN", v.in.ogr_flags=c("o", "overwrite"))
  x <- as.character(Snapped_SSP_DSN[i,]@coords[,1][[1]])
  y <- as.character(Snapped_SSP_DSN[i,]@coords[,2][[1]])
  execGRASS("r.water.outlet",flags="overwrite", drainage="Drainage", basin="UC_rast", easting=x, northing=y)
  execGRASS("r.to.vect",flags=c("overwrite"), parameters=list(input="UC_rast", output="UC_Poly", feature="area"))
  UCpoly <- readVECT6("UC_Poly")
  proj4string(UCpoly)<-proj4string(DEM)
  UC <- c(UC, UCpoly)
  rm(UCpoly)
  UC_rast <- readRAST6("UC_rast")
  proj4string(UC_rast)<-proj4string(DEM)
  UC_rast <- raster(UC_rast)
  
  ## The snapped SSP will be connected to the equivalent DSN as new nodes 
  execGRASS("v.net", flags="overwrite", parameters=list(input="DSN", points="Snapped_SSP_DSN", output="Stream_network", operation="connect", thresh=0.5))
  
  ## The iso-stream sections of 10km will be identified, extracted and imported to R
  execGRASS("v.net.iso", flags="overwrite", parameters=list(input="Stream_network", output="Zone_stream_10000m", ccats="1", costs=as.integer(10000)))
  execGRASS("v.extract", flags="overwrite", parameters=list(list="1", input="Zone_stream_10000m", output="Stream_within_10000m",  type="line"))
  Stream_within_10000m <- readVECT6("Stream_within_10000m")
  proj4string(Stream_within_10000m)<-proj4string(DEM)
  
  ## Buffer polygons of 150m will be computed for each of the iso-stream sections.
  # The buffer polygons will be masked with previously extracted raster UC. 
  # The functions 'gBuffer' and 'mask' of the 'rgeos' and 'raster' packages, respectively,
  # of R will be used for the purposes. The masked UC is the desired URC in
  # raster format. The raster URC will be exported to GRASS and converted into polygons
  Buffer10000m_100m <- gBuffer(Stream_within_10000m, byid=FALSE, id=NULL, width=100, quadsegs=5, capStyle="ROUND", joinStyle="ROUND")
  URC_rast <- mask(UC_rast, Buffer10000m_100m)
  rm(UC_rast)
  writeRaster(URC_rast, file="URC_rast", format="GTiff", overwrite=TRUE)
  rm(URC_rast)
  execGRASS("r.in.gdal",flags=c("o","overwrite"), parameters=list(input=paste(path,"URC_rast.tif", sep="/"), output="URC_rast"))
  execGRASS("r.to.vect",flags=c("overwrite"), parameters=list(input="URC_rast", output="URC_Poly", feature="area"))
  URCpoly <- readVECT6("URC_Poly")
  proj4string(URCpoly)<-proj4string(DEM)
  URC <- c(URC, URCpoly)
  rm(URC_poly)
  
  ## The progressbar will report the current status of the whole process
  Sys.sleep(0.1)
  setTkProgressBar(ProgressBar, i, label=paste(round(i/nrow(Snapped_SSP_DSN)*100, 0), "% done"))
  
}

save(UC, file="UC.Rdata")
save(URC, file="URC.Rdata")
# The computed upstream UC and riparian corridor polygons are stored in the 'UC' and 'URC' lists.
# They can be applied for further analyses in freshwater research, e.g. landuse matrices derivation etc.

########################################################################
########################################################################
#                              (end)                                   #
########################################################################
########################################################################