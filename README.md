## ATRIC ##

Automated Accumulation Threshold selection and RIparian Corridor delineation (ATRIC) is a combination of two algorithms that allows for automated accumulation threshold selection for stream network extraction, and watershed and riparian corridor delineation for given points on streams, from digital elevation models. It is also useful for quantification of stressors such as land use and climate change components, pollutants, insecticide and pesticide runoff within the watersheds and upstream riparian zones of stream sampling points, where stream communities, e.g. macroinvertebrates are sampled. ATRIC is the first open source tool of its kind, the algorithms are developed by co-interfacing R and GRASS GIS.

For details on ATRIC, please read the published paper: http://www.sciencedirect.com/science/article/pii/S1364815214003077 (copies available upon request).


## Note for users ##

You need to have the latest version of R and GRASS GIS version 6.4 installed in your computer. Details on how to initiate GRASS GIS from R interface for different operating systems are available in the commented R script "ATRIC_script.R". Note that a new stable version of GRASS GIS (7.0) is now available where parameter names for some modules have been changed. Consequently, the script might not work if you install or update your GRASS GIS 6.4 with GRASS GIS 7.0. I will release a script for GRASS GIS 7.0 soon.

## Reproduce results of the paper ##

If you would like to reproduce our results in the paper or just would like to get familiar with the functionalities of ATRIC, download the "ATRIC_Data.zip" and run the script.  


