## Spatial and Spatiotemporal Data Analysis in R

UseR! Workshop, Jul 9, 2019, 

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

Workshop program:

* 14:00-15:15: Introduction, spatial, spatio-temporal data, data cubes (Edzer); [slides](https://edzer.github.io/UseR2019/part1.html)
* 15:15-16:00: Exercises (and break) with help from Angela, Edzer, Roger (break: 15:30-16:00)
* 16:00-16:15: Recap of the exercises
* 16:15-17:30: Spatial modelling, spatial weights, disease modelling (Roger)
[slides](https://edzer.github.io/UseR2019/part2.html)

Potentially useful materials for preparing:

* UseR! 2017 workshop "Spatial Data in R: New Directions", mostly discussing `sf`: [slides](https://edzer.github.io/UseR2017/)
* Geostat summer school workshop "New R packages for spatial and spatiotemporal vector and raster data": [video](https://www.youtube.com/watch?v=yhpkx_xO-LE&list=PLXUoTpMa_9s3T-K7m8LO3Mf29g9E4EJLs)
* Spatial Data Science with R: [book](https://r-spatial.org/book)
* rstudio::conf presentations _Tidy spatial data analysis_ (2018; [slides](https://edzer.github.io/rstudio_conf/index.html#1), [video](https://www.rstudio.com/resources/videos/tidy-spatial-data-analysis/)), _Spatial data science in the Tidyverse_ (2019; [slides](https://edzer.github.io/rstudio_conf/2019/#1), [video](https://resources.rstudio.com/rstudio-conf-2019/spatial-data-science-in-the-tidyverse))
* [Bergen ML slides](https://github.com/bergen-ml/2019-02-19-bivand)
* [Paris workshop slides](https://github.com/rsbivand/sew19)

Datasets that will be used:


