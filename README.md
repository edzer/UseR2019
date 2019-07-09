## Spatial and Spatiotemporal Data Analysis in R

UseR! Workshop, Jul 9, 2019, Edzer Pebesma, Roger Bivand, Angela Li (helper)

This tutorial dives into some of the modern spatial and
spatiotemporal analysis packages available in R. It will show how
support, the spatial size of the area to which a data value refers,
plays a role in spatial analysis, and how this is handled in R. It
will show how package stars complements package sf for handling
spatial time series, raster data, raster time series, and more
complex multidimensional data such as dynamic origin-destination
matrices. It will also show how stars handles out-of- memory
datasets, with an example that uses Sentinel-2 satellite time
series. This will be connected to analysing the data with packages
that assume spatial processes as their modelling framework,
including gstat, spdep, and R-INLA. Familiarity with package sf
and the tidyverse will be helpful for taking this tutorial.

## Workshop program

* 14:00-15:30: Introduction, spatial, spatio-temporal data, data cubes + exercises (Edzer); [materials](https://edzer.github.io/UseR2019/part1.html), [solutions for exercises](https://github.com/edzer/UseR2019/blob/master/solutions_1.R)
* 15:30-16:00: Break and time for questions
* 16:00-17:30: Spatial modelling, spatial weights, spatial regression + exercises (Roger);
[materials](https://edzer.github.io/UseR2019/part2.html), [R script](https://raw.githubusercontent.com/edzer/UseR2019/master/part2.R)

## Packages used in this workshop

#### Part 1

* `sf`
* `stars`
* `gstat`
* `units`
* `tidyverse`
* `xts`
* `viridis`
* `abind`

#### Part 2 (new packages only)

* `spatstat` 
* `spdep` 
* `tmap` 
* `spatialreg` 
* `igraph` 
* `hglm` 
* `metafor` 
* `sp` 
* `spData` 
* `RANN`

If you are having trouble installing these packages, you may need to update their geospatial library dependencies (GDAL, GEOS, PROJ.4, or UDUNITS). Please see [this guide from Data Carpentry](https://datacarpentry.org/geospatial-workshop/setup.html) for more information on how to install these dependencies.

## Reference materials

* `sf` [cheatsheet](https://github.com/rstudio/cheatsheets/blob/master/sf.pdf) and [reference website](https://r-spatial.github.io/sf/index.html)
* `stars` [reference website](https://r-spatial.github.io/stars/)
* UseR! 2017 workshop "Spatial Data in R: New Directions", mostly discussing `sf`: [slides](https://edzer.github.io/UseR2017/)
* Geostat summer school workshop "New R packages for spatial and spatiotemporal vector and raster data": [video](https://www.youtube.com/watch?v=yhpkx_xO-LE&list=PLXUoTpMa_9s3T-K7m8LO3Mf29g9E4EJLs)
* Spatial Data Science with R: [book](https://r-spatial.org/book)
* rstudio::conf presentations _Tidy spatial data analysis_ (2018; [slides](https://edzer.github.io/rstudio_conf/index.html#1), [video](https://www.rstudio.com/resources/videos/tidy-spatial-data-analysis/)), _Spatial data science in the Tidyverse_ (2019; [slides](https://edzer.github.io/rstudio_conf/2019/#1), [video](https://resources.rstudio.com/rstudio-conf-2019/spatial-data-science-in-the-tidyverse))
* [Bergen ML slides](https://github.com/bergen-ml/2019-02-19-bivand)
* [Paris workshop slides](https://github.com/rsbivand/sew19)
