

## ---- echo=TRUE----------------------------------------------------------
suppressPackageStartupMessages(library(spatstat))
intenfun <- function(x, y) 200 * x
set.seed(1)
(ihpp <- rpoispp(intenfun, lmax = 200))


## ---- echo=TRUE, out.width='90%', fig.align='center'---------------------
plot(density(ihpp), axes=TRUE)
points(ihpp, col="green2", pch=19)


## ---- echo=TRUE, out.width='90%', fig.align='center'---------------------
opar <- par(mfrow=c(1,2))
plot(envelope(ihpp, Kest, verbose=FALSE), main="Homogeneous")
plot(envelope(ihpp, Kinhom, verbose=FALSE), main="Inhomogeneous")
par(opar)


## ---- echo=TRUE, out.width='90%', fig.align='center'---------------------
library(sf)
#st_as_sf(ihpp) %>% dplyr::filter(label == "point") -> sf_ihpp
st_as_sf(ihpp) ->.; .[.$label == "point",] -> sf_ihpp
crds <- st_coordinates(sf_ihpp)
sf_ihpp$x <- crds[,1]
sf_ihpp$y <- 100 + 50 * sf_ihpp$x + 20 * rnorm(nrow(sf_ihpp))
plot(sf_ihpp[,"y"], pch=19)


## ---- echo=TRUE, out.width='90%', fig.align='center'---------------------
suppressPackageStartupMessages(library(gstat))
vg0 <- variogram(y ~ 1, sf_ihpp)
vg1 <- variogram(y ~ x, sf_ihpp)
opar <- par(mfrow=c(1,2))
plot(gamma ~ dist, vg0, ylim=c(0,550), main="Trend ignored")
s <- seq(0, 0.5, 0.01)
p0 <- predict(loess(gamma ~ dist, vg0, weights=np), newdata=data.frame(dist=s), se=TRUE)
lines(s, p0$fit)
lines(s, p0$fit-2*p0$se.fit, lty=2)
lines(s, p0$fit+2*p0$se.fit, lty=2)
plot(gamma ~ dist, vg1, ylim=c(0,550), main="Trend included")
p1 <- predict(loess(gamma ~ dist, vg1, weights=np), newdata=data.frame(dist=s), se=TRUE)
lines(s, p1$fit)
lines(s, p1$fit-2*p1$se.fit, lty=2)
lines(s, p1$fit+2*p1$se.fit, lty=2)
par(opar)


## ---- echo=TRUE, results='hide', message=FALSE---------------------------
suppressPackageStartupMessages(library(spdep))
nb_tri <- tri2nb(crds)


## ---- echo=TRUE----------------------------------------------------------
(nb_soi <- graph2nb(soi.graph(nb_tri, crds), sym=TRUE))


## ---- echo=TRUE, out.width='90%', fig.align='center'---------------------
plot(nb_soi, crds)


## ---- echo=TRUE----------------------------------------------------------
comps <- n.comp.nb(nb_soi)
sf_ihpp$comps <- comps$comp.id
comps$nc


## ---- echo=TRUE, warning=FALSE, message=FALSE----------------------------
lwB <- nb2listw(nb_soi, style="B")
out <- broom::tidy(moran.test(sf_ihpp$y, listw=lwB, randomisation=FALSE, alternative="two.sided"))[1:5]
names(out)[1:3] <- c("Moran's I", "Expectation", "Variance"); out


## ---- echo=TRUE----------------------------------------------------------
lm_obj <- lm(y ~ x, data=sf_ihpp)
out <- broom::tidy(lm.morantest(lm_obj, listw=lwB, alternative="two.sided"))[1:5]
names(out)[1:3] <- c("Moran's I", "Expectation", "Variance"); out


## ---- echo=TRUE----------------------------------------------------------
lmor0 <- localmoran(sf_ihpp$y, listw=lwB, alternative="two.sided")
lmor1 <- as.data.frame(localmoran.sad(lm_obj, nb=nb_soi, style="B", alternative="two.sided"))
sf_ihpp$z_value <- lmor0[,4]
sf_ihpp$z_lmor1_N <- lmor1[,2]
sf_ihpp$z_lmor1_SAD <- lmor1[,4]


## ---- echo=TRUE, warning=FALSE, out.width='90%', fig.align='center', width=7, height=4----
suppressPackageStartupMessages(library(tmap))
tm_shape(sf_ihpp) + tm_symbols(col=c("z_value", "z_lmor1_N", "z_lmor1_SAD"), midpoint=0) + tm_facets(free.scales=FALSE, nrow=1) + tm_layout(panel.labels=c("No trend", "Trend, normal", "Trend, SAD"))


## ---- echo=TRUE----------------------------------------------------------
library(sf)
nc <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
st_crs(nc) <- "+proj=longlat +datum=NAD27"
row.names(nc) <- as.character(nc$FIPSNO)
head(nc)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
library(spdep)
gal_file <- system.file("weights/ncCR85.gal", package="spData")[1]
ncCR85 <- read.gal(gal_file, region.id=nc$FIPSNO)
ncCR85


## ---- echo=TRUE, warning=TRUE, out.width='90%', fig.align='center', width=7, height=4----
plot(st_geometry(nc), border="grey")
plot(ncCR85, st_centroid(st_geometry(nc), of_largest_polygon), add=TRUE, col="blue")


## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
nc$rand <- rnorm(nrow(nc))
lw <- nb2listw(ncCR85, style="B")
moran.test(nc$rand, listw=lw, alternative="two.sided")


## ---- echo=TRUE----------------------------------------------------------
nc$LM <- as.numeric(interaction(nc$L_id, nc$M_id))
alpha <- 1
beta <- 0.5
sigma <- 2
nc$trend <- alpha + beta*nc$LM + sigma*nc$rand
moran.test(nc$trend, listw=lw, alternative="two.sided")


## ---- echo=TRUE----------------------------------------------------------
lm.morantest(lm(trend ~ LM, nc), listw=lw, alternative="two.sided")


## ---- echo=TRUE----------------------------------------------------------
aggLM <- aggregate(nc[,"LM"], list(nc$LM), head, n=1)
(aggnb <- poly2nb(aggLM))


## ---- echo=TRUE----------------------------------------------------------
set.seed(1)
LMrand <- rnorm(nrow(aggLM))


## ---- echo=TRUE----------------------------------------------------------
moran.test(LMrand, nb2listw(aggnb, style="B"))


## ---- echo=TRUE----------------------------------------------------------
nc$LMrand <- LMrand[match(nc$LM, aggLM$LM)]
plot(nc[,"LMrand"])


## ---- echo=TRUE----------------------------------------------------------
moran.test(nc$LMrand, listw=lw, alternative="two.sided")


## ---- echo=TRUE----------------------------------------------------------
library(sf)
# if(!require("spData", quietly=TRUE)) install.packages("spData")
# if(!require("spDataLarge", quietly=TRUE)) install.packages("spDataLarge",
#   repos = "https://nowosad.github.io/drat/", type = "source")
data(pol_pres15, package="spDataLarge")
head(pol_pres15[, c(1, 4, 6)])


## ---- echo=TRUE, eval=FALSE----------------------------------------------
## ?spDataLarge::pol_pres15


## ---- echo=TRUE----------------------------------------------------------
suppressPackageStartupMessages(library(spdep))
system.time(nb_q <- poly2nb(pol_pres15, queen=TRUE))


## ---- echo=TRUE----------------------------------------------------------
nb_q


## ---- echo=TRUE----------------------------------------------------------
opar <- par(mar=c(0,0,0,0)+0.5)
plot(st_geometry(pol_pres15), border="grey", lwd=0.5)
coords <- st_centroid(st_geometry(pol_pres15),
  of_largest_polygon=TRUE)
plot(nb_q, coords=st_coordinates(coords), add=TRUE, points=FALSE, lwd=0.5)
par(opar)


## ---- echo=TRUE----------------------------------------------------------
system.time({
  fB1 <- st_intersects(st_as_sfc(lapply(st_geometry(pol_pres15), function(x) {
    st_as_sfc(st_bbox(x))[[1]]
  })))
  fB1a <- lapply(seq_along(fB1), function(i) fB1[[i]][fB1[[i]] > i])
  fB1a <- fB1a[-length(fB1a)]
  nb_sf_q1 <- poly2nb(pol_pres15, queen=TRUE, foundInBox=fB1a)
})


## ---- echo=TRUE----------------------------------------------------------
all.equal(nb_q, nb_sf_q1, check.attributes=FALSE)


## ---- echo=TRUE----------------------------------------------------------
n.comp.nb(nb_q)$nc


## ---- echo=TRUE----------------------------------------------------------
lw <- nb2listw(nb_q, style="B")


## ---- echo=TRUE, results='hide', message=FALSE---------------------------
library(Matrix)
W <- as(lw, "CsparseMatrix")


## ---- echo=TRUE----------------------------------------------------------
isTRUE(all.equal(W, t(W)))


## ---- echo=TRUE----------------------------------------------------------
image(W)


## ---- echo=TRUE----------------------------------------------------------
WW <- W %*% W
image(WW)


## ---- echo=TRUE----------------------------------------------------------
W3 <- WW %*% W
image(W3)


## ---- echo=TRUE----------------------------------------------------------
W4 <- W3 %*% W
image(W4)


## ---- echo=TRUE, message=FALSE-------------------------------------------
library(igraph)
g1 <- graph.adjacency(W, mode="undirected")
class(g1)


## ---- echo=TRUE----------------------------------------------------------
B1 <- get.adjacency(g1)
mat2listw(B1)$neighbours


## ---- echo=TRUE----------------------------------------------------------
c1 <- clusters(g1)
c1$no


## ---- echo=TRUE----------------------------------------------------------
is.connected(g1)


## ---- echo=TRUE----------------------------------------------------------
dg1 <- diameter(g1)
dg1


## ---- echo=TRUE----------------------------------------------------------
sp_mat <- shortest.paths(g1)
max(sp_mat)
str(sp_mat)


## ---- echo=TRUE----------------------------------------------------------
gra <- R2BayesX::nb2gra(ncCR85)
str(gra)


## ---- echo=TRUE----------------------------------------------------------
tf <- tempfile()
nb2INLA(tf, ncCR85)
file.size(tf)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
coords <- st_centroid(st_geometry(nc), of_largest_polygon)
(knn_5_nb <- knn2nb(knearneigh(st_coordinates(coords), k=5)))


## ---- echo=TRUE----------------------------------------------------------
klw <- nb2listw(knn_5_nb, style="B")
kW <- as(klw, "CsparseMatrix")
isTRUE(all.equal(kW, t(kW)))


## ---- echo=TRUE----------------------------------------------------------
image(kW)


## ---- echo=TRUE, messages=FALSE------------------------------------------
library(igraph)
g1 <- graph.adjacency(kW, mode="directed")
B1 <- get.adjacency(g1)
mat2listw(B1)$neighbours


## ---- echo=TRUE----------------------------------------------------------
diameter(g1)


## ---- echo=TRUE----------------------------------------------------------
nc_sps <- shortest.paths(g1)
mr <- which.max(apply(nc_sps, 2, max))
nc$sps1 <- nc_sps[,mr]
plot(nc[,"sps1"], breaks=0:21)

## ---- echo=TRUE----------------------------------------------------------
plot(nc$sps1, c(st_distance(coords[mr], coords))/1000, xlab="shortest path count", ylab="km distance")


## ---- echo=TRUE----------------------------------------------------------
ncCR85a <- ncCR85
attr(ncCR85a, "region.id") <- as.character(nc$CRESS_ID)
nc_gra <- R2BayesX::nb2gra(ncCR85a)
nc_tf <- tempfile()
nb2INLA(nc_tf, ncCR85)
nc_lw <- nb2listw(ncCR85, style="B")
nc_W <- as(nc_lw, "CsparseMatrix")
nc_mat <- listw2mat(nc_lw)


## ---- echo=TRUE----------------------------------------------------------
nc$ft.SID74 <- sqrt(1000)*(sqrt(nc$SID74/nc$BIR74) + sqrt((nc$SID74+1)/nc$BIR74))
nc$ft.NWBIR74 <- sqrt(1000)*(sqrt(nc$NWBIR74/nc$BIR74) + sqrt((nc$NWBIR74+1)/nc$BIR74))
tm_shape(nc) + tm_fill(c("ft.SID74", "ft.NWBIR74"))


## ---- echo=TRUE----------------------------------------------------------
plot(ft.SID74 ~ ft.NWBIR74, nc)


## ---- echo=TRUE----------------------------------------------------------
moran.test(nc$ft.SID74, nc_lw, alternative="two.sided", randomisation=FALSE)


## ---- echo=TRUE----------------------------------------------------------
lm.morantest(lm(ft.SID74 ~ 1, weights=BIR74, data=nc), nc_lw, alternative="two.sided")


## ---- echo=TRUE----------------------------------------------------------
lm.morantest(lm(ft.SID74 ~ ft.NWBIR74, data=nc), nc_lw, alternative="two.sided")


## ---- echo=TRUE----------------------------------------------------------
lm.morantest(lm(ft.SID74 ~ ft.NWBIR74, weights=BIR74, data=nc), nc_lw, alternative="two.sided")


## ---- echo=TRUE, results='hide', message=FALSE, warning=FALSE------------
library(spatialreg)
m1 <- spautolm(ft.SID74 ~ ft.NWBIR74, weights=BIR74, data=nc, listw=nc_lw, family="SAR")


## ---- echo=TRUE----------------------------------------------------------
summary(m1)


## ---- echo=TRUE----------------------------------------------------------
nc$SAR_ssre <- 2*as.vector((m1$lambda * nc_W) %*% m1$Y - (m1$lambda * nc_W) %*% (m1$X %*% m1$fit$coefficients))
tm_shape(nc) + tm_fill(c("ft.SID74", "SAR_ssre"), midpoint=c(NA, 0))


## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
m1a <- errorsarlm(ft.SID74 ~ ft.NWBIR74, weights=BIR74, data=nc, listw=nc_lw)
summary(m1a, Hausman=TRUE)


## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
m1b <- errorsarlm(ft.SID74 ~ ft.NWBIR74, weights=BIR74, data=nc, listw=nc_lw, Durbin=TRUE)
summary(m1b, Hausman=TRUE)


## ---- echo=TRUE, message=FALSE, warning=FALSE----------------------------
summary(impacts(m1b))


## ---- echo=TRUE, results='hide', message=FALSE, warning=FALSE------------
library(hglm)
E <- nc$BIR74 * sum(nc$SID74)/sum(nc$BIR74)
HGLM_iid <- hglm(fixed=SID74 ~ ft.NWBIR74, random= ~ 1|CRESS_ID, offset=log(E), weights=BIR74,
                 data=nc, family=poisson(link=log))


## ---- echo=TRUE----------------------------------------------------------
ranef_iid <- unname(summary(HGLM_iid, print.ranef=TRUE)$RandCoefMat)
metafor::forest(ranef_iid[,1], ranef_iid[,2], subset=order(ranef_iid[,1], decreasing=TRUE), 
        slab=NA, annotate=FALSE, lty=c("solid","blank"), pch=19, psize=2, cex.lab=1, cex.axis=1)


## ---- echo=TRUE, results='hide', message=FALSE, cache=TRUE, warning=FALSE----
HGLM_sar <- hglm(fixed=SID74 ~ ft.NWBIR74, random= ~ 1|CRESS_ID, offset=log(E), weights=BIR74, 
                 data=nc, family=poisson(link=log), rand.family=SAR(D=nc_W))
ranef_sar <- unname(summary(HGLM_sar, print.ranef=TRUE)$RandCoefMat)
metafor::forest(ranef_sar[,1], ranef_sar[,2], subset=order(ranef_sar[,1], decreasing=TRUE), 
        slab=NA, annotate=FALSE, lty=c("solid","blank"), pch=19, psize=2, cex.lab=1, cex.axis=1)


## ---- echo=TRUE----------------------------------------------------------
nc$HGLM_re <- ranef_iid[,1]
nc$HGLM_ss_SAR <- ranef_sar[,1]
tm_shape(nc) + tm_fill(c("HGLM_re", "HGLM_ss_SAR"), midpoint=c(0), title="Poisson HGLM RE") +
  tm_facets(free.scales=FALSE) + tm_layout(panel.labels=c("IID", "SAR SSRE"))


## ---- echo=TRUE----------------------------------------------------------
m1c <- spautolm(ft.SID74 ~ ft.NWBIR74, weights=BIR74, data=nc, listw=nc_lw, family="CAR")
summary(m1c)


## ---- echo=TRUE----------------------------------------------------------
nc$CAR_ssre <- as.vector((m1c$lambda * nc_W) %*% m1c$Y - 
                           (m1c$lambda * nc_W) %*% (m1c$X %*% m1c$fit$coefficients))
tm_shape(nc) + tm_fill(c("SAR_ssre", "CAR_ssre"), midpoint=c(0), title="Gauss ML RE") +
  tm_facets(free.scales=FALSE) + tm_layout(panel.labels=c("SAR SSRE", "CAR SSRE"))


## ---- echo=TRUE, results='hide', message=FALSE, cache=TRUE, warning=FALSE----
HGLM_car <- hglm(fixed=SID74 ~ ft.NWBIR74, random= ~ 1|CRESS_ID, offset=log(E), weights=BIR74, 
                 data=nc, family=poisson(link=log), rand.family=CAR(D=nc_W))
ranef_car <- unname(summary(HGLM_car, print.ranef=TRUE)$RandCoefMat)
metafor::forest(ranef_car[,1], ranef_car[,2], subset=order(ranef_car[,1], decreasing=TRUE), 
        slab=NA, annotate=FALSE, lty=c("solid","blank"), pch=19, psize=2, cex.lab=1, cex.axis=1)


## ---- echo=TRUE----------------------------------------------------------
nc$HGLM_ss_CAR <- ranef_car[,1]
tm_shape(nc) + tm_fill(c("HGLM_ss_CAR", "HGLM_ss_SAR"), midpoint=c(0), title="Poisson HGLM RE") +
  tm_facets(free.scales=FALSE) + tm_layout(panel.labels=c("CAR SSRE", "SAR SSRE"))


## ---- echo=TRUE, results='hide', message=FALSE, cache=TRUE, warning=FALSE----
library(mgcv)
names(aggnb) <- as.character(aggLM$Group.1)
nc$LM <- as.factor(nc$LM)
GAM_mrf <- gam(SID74 ~ s(ft.NWBIR74) + s(LM, bs="mrf", xt=list(nb=aggnb)), offset=log(E), weights=BIR74, data=nc, family=poisson(link=log))
summary(GAM_mrf)


## ---- echo=TRUE----------------------------------------------------------
plot(GAM_mrf)


## ---- echo=TRUE----------------------------------------------------------
GAM_mrf_re <- predict(GAM_mrf, type="terms", se=TRUE)
metafor::forest(GAM_mrf_re$fit[,2], GAM_mrf_re$se.fit[,2], subset=order(GAM_mrf_re$fit[,1], decreasing=TRUE), 
        slab=NA, annotate=FALSE, lty=c("solid","blank"), pch=19, psize=2, cex.lab=1, cex.axis=1)


## ---- echo=TRUE----------------------------------------------------------
nc$GAM_mrf_re <- GAM_mrf_re$fit[,2]
tm_shape(nc) + tm_fill(c("GAM_mrf_re"), midpoint=c(0), title="Poisson GAM MRF RE")


## ---- echo=TRUE----------------------------------------------------------
NY8 <- st_read(system.file("shapes/NY8_utm18.shp", package="spData"))


## ---- echo=TRUE----------------------------------------------------------
tm_shape(NY8) + tm_fill("Z")


## ---- echo=TRUE----------------------------------------------------------
NY_nb <- poly2nb(NY8)
NY_lw <- nb2listw(NY_nb, style="B")


## ---- echo=TRUE----------------------------------------------------------
mod1 <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=NY8, family="CAR", listw=NY_lw, weights=POP8)
summary(mod1)


## ----sI, echo = TRUE-----------------------------------------------------
sessionInfo()

