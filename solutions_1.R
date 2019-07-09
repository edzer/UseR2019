#Part A:

#1. start R, install package `sf` if not already present, and load it

library(sf)

#2. load the `nc` sample dataset as done above

nc = read_sf(system.file("gpkg/nc.gpkg", package="sf")) # read as sf-tibble

#3. use `st_intersects` to discover which counties in `nc` intersect with county `Rowan`; experiment with the order of arguments to `st_intersects`

rowan = nc[nc$NAME == "Rowan",]
st_intersects(nc, rowan)
st_intersects(rowan, nc)

#4. try to understand the object returned by `st_intersects`

i = st_intersects(rowan, nc)
class(i)
str(i)

#5. look up which methods are available for objects of class `sgbp` (`methods(class = "sgbp")`), and try some of them

methods(class = "sgbp")
t(i)
as.matrix(i)
as.data.frame(i)
!i

#6. does `Rowan` intersect with itself? 

which(nc$NAME == "Rowan") %in% i[[1]]

#7. Find all counties that touch `Rowan` by using `st_touches`; does `Rowan` touch itself?

i2 = st_touches(rowan, nc)
which(nc$NAME == "Rowan") %in% i2[[1]]

#8. How can we find what `intersects` and `touches` mean, more formally?

# The help page of ?st_intersects points to https://en.wikipedia.org/wiki/DE-9IM

## B:
#1. Start R, install package `stars` (if not already present), and load it

library(stars)

#2. load the Landsat image test set, as above

tif = system.file("tif/L7_ETMs.tif", package = "stars")
x = read_stars(tif)

#3. Create an RGB composite plot using argument `rgb = c(3,2,1)` for RGB, and `rgb = c(4,3,2)` for false color.

plot(x, rgb = c(3, 2, 1))
plot(x, rgb = c(4, 3, 2))

#4. Use `x6 = split(x, "band")` to create a 2-dimensional raster with 6 attributes

(x6 = split(x, "band"))

#5. Plot `x6`. What has changed?

plot(x6)
# only a single band is plotted!

#6. Compute the mean of all six bands by adding all six attributes and dividing by 6, and assigning the resulting matrix as a new attribute in `x6`

x6$mean = (x6[[1]] + x6[[2]] + x6[[3]] + x6[[4]] + x6[[5]] + x6[[6]])/6

#7. As an alternative, compute the mean of all six bands by applying function `mean` over the `x` and `y` dimensions (essentially reducing "band"), by using `st_apply`; compare the results of the two approaches

xm = st_apply(x, c("x", "y"), mean)
all.equal(xm[[1]], x6$mean)

## C:
# 1. rasterize the `nc` dataset above, using `st_rasterize`. 

# Which attribute is retained? How could you do this for another attribute?

st_rasterize(nc)

# Which attribute is retained? 

# AREA (the first)

# How could you do this for another attribute?

# name the argument:
st_rasterize(nc["SID74"])

# 2. Modify the rasterization such that you get 30 rows and 60 columns (hint: look at the exercises help of `?st_rasterize`)

st_rasterize(nc["SID74"], st_as_stars(st_bbox(nc), nx = 30, ny = 60, values = NA_real_))

# 3. Look up what `st_crop` does, and run its examples.
# 4. Using the precipitation grid data shown above, find the time series with for each county in North Carolina the _maximum_ precipitation, and plot this. Note that
#   - `nc` has a different datum from `prec`, use `st_transform` for datum transformation
#   - the time values should be hourly, but may come out differently

a = aggregate(prec, st_transform(nc, st_crs(prec_slice)), max)
plot(a, max.plot = 23)

# 5. read in the file `system.file("nc/tos_O1_2001-2002.nc", package = "stars")

r = read_stars(system.file("nc/tos_O1_2001-2002.nc", package = "stars"))

# 6. what is the bounding box of this data set? What is the CRS, and the temporal reference system?

st_bbox(r)
st_crs(r)

# time reference: PCICt 360-day calendar

# 7. obtain the temperature time series of this dataset at point `POINT(200 10)`, using `pt = st_sfc(st_point(c(200, 10)))`

pt = st_point(c(200, 10))
r[st_sfc(pt)]

# or 

aggregate(r, st_sfc(pt), identity)
