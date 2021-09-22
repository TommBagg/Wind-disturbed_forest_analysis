# Title: Stored volume index 22/09/2021

# Script for computing the stored volume, an index assessing the volume between a smooth layer and the original DSM for the assessment of forest affected by storms
# Author details: Tommaso Baggio Contact details: tommaso.baggio@phd.unipd.it; tbaggio93@gmail.com

# Copyright statement: This script is the product of the work of Tommaso Baggio - 2021

# THE SCRIPT IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Necessary input: point cloud (.las format) from LiDAR or photogrammetric survey
#                  shapefile reporting the areas of interest for which the algorithm compute the analysis  

library(raster)
library(rgdal)
library(lidR)


# setwd("/home/tommaso/computer_lavoro/phd/valanghe_schianti/temporal_analysis_disentis/uniform_layer")
setwd("D:/phd/valanghe_schianti/temporal_analysis_disentis/uniform_layer")
DTM <- raster("dtm_2019_0p2.tif"); plot(DTM)
# areas_H0 <- readOGR("/home/tommaso/computer_lavoro/phd/valanghe_schianti/temporal_analysis_disentis/areas", "area_analysis_21781", stringsAsFactors = T); plot(areas_H0, add=TRUE)
areas_H0 <- readOGR("D:/phd/valanghe_schianti/temporal_analysis_disentis/areas", "area_analysis_21781", stringsAsFactors = T); plot(areas_H0, add=TRUE)

res_max <- 2.5       #resolution to compute the maximum values
res_photog <- 0.50     #resolution to create the DSM based photogrammetry without standing trees


#### 1991 ####

las <- readLAS("di91_lieg_dom1_20cm.laz", select = "xyz")  #input photogrammetric dense point cloud

###########################################

#Remove outliers SOR
las <- classify_noise(las, sor(8,2))
las <- filter_poi(las, Classification != LASNOISE)
print(las)
# Calculate the DSM without trees #CSF
las <- classify_ground(las, csf(sloop_smooth = TRUE, class_threshold = 3.0, cloth_resolution = 0.5))
# writeLAS(las, "DSM_no_trees_classification.las")
#plot(las, color = "Classification")
las_no_trees_CSF_canopy <- filter_poi(las, Classification == 2)
DSM_no_trees <- grid_metrics(las_no_trees_CSF_canopy, res = res_photog, ~max(Z)); plot(DSM_no_trees)
DSM_no_trees <- focal(DSM_no_trees, w=matrix(1,51,51), fun=mean, NAonly=TRUE, na.rm=TRUE, pad=TRUE); plot(DSM_no_trees)
# DSM_no_trees <- grid_terrain(las, res = res_photog, tin())

#compute CHM
DTM_res <- resample(DTM, DSM_no_trees, method='bilinear');plot(DTM_res)
CHM <- DSM_no_trees - DTM_res;
CHM <- CHM * (CHM > 0); plot(CHM)

#### compute H0
factor_max <- res_max/((res(CHM)[1]+res(CHM)[2])/2)
CHM_max <- aggregate(CHM, fact=factor_max, expand=TRUE, fun=max); res(CHM_max); plot(CHM_max)
CHM_max_res <- resample(CHM_max,CHM, method='bilinear')
CHM_max_res <- CHM_max_res * (CHM_max_res > 0); plot(CHM_max_res)
difference <- abs(CHM_max_res - CHM); plot(difference)

values <- extract(difference, areas_H0, fun=sum)
areas_H0$diff <- as.numeric(values)
areas_H0$volume <- areas_H0$diff*((res(CHM)[1]+res(CHM)[2])/2)^2
areas_H0$area <- area(areas_H0)
areas_H0$H0_volume <- areas_H0$volume / areas_H0$area

plot(difference);plot(areas_H0, add=TRUE)

# Outputs 
writeRaster(difference, "difference_1991.tif", overwrite=TRUE)
# writeOGR(areas_H0, dsn = wd, layer = "areas_H0_photogr_las",driver="ESRI Shapefile")
write.table(areas_H0, "volume_difference_H0_1991.txt")


#### 2001 ####

las <- readLAS("di01_lieg_dom1_20cm.laz", select = "xyz")  #input photogrammetric dense point cloud

###########################################

#Remove outliers SOR
las <- classify_noise(las, sor(8,2))
las <- filter_poi(las, Classification != LASNOISE)
print(las)
# Calculate the DSM without trees #CSF
las <- classify_ground(las, csf(sloop_smooth = TRUE, class_threshold = 3.0, cloth_resolution = 0.5))
# writeLAS(las, "DSM_no_trees_classification.las")
# plot(las_no_trees_CSF, color = "Classification")
las_no_trees_CSF_canopy <- filter_poi(las, Classification == 2)
DSM_no_trees <- grid_metrics(las_no_trees_CSF_canopy, res = res_photog, ~max(Z)); plot(DSM_no_trees)
DSM_no_trees <- focal(DSM_no_trees, w=matrix(1,51,51), fun=mean, NAonly=TRUE, na.rm=TRUE, pad=TRUE); plot(DSM_no_trees)
# DSM_no_trees <- grid_terrain(las, res = res_photog, tin())

#compute CHM
DTM_res <- resample(DTM, DSM_no_trees, method='bilinear');plot(DTM_res)
CHM <- DSM_no_trees - DTM_res;
CHM <- CHM * (CHM > 0); plot(CHM)

#### compute H0
factor_max <- res_max/((res(CHM)[1]+res(CHM)[2])/2)
CHM_max <- aggregate(CHM, fact=factor_max, expand=TRUE, fun=max); res(CHM_max); plot(CHM_max)
CHM_max_res <- resample(CHM_max,CHM, method='bilinear')
CHM_max_res <- CHM_max_res * (CHM_max_res > 0); plot(CHM_max_res)
difference <- abs(CHM_max_res - CHM); plot(difference)

values <- extract(difference, areas_H0, fun=sum)
areas_H0$diff <- as.numeric(values)
areas_H0$volume <- areas_H0$diff*((res(CHM)[1]+res(CHM)[2])/2)^2
areas_H0$area <- area(areas_H0)
areas_H0$H0_volume <- areas_H0$volume / areas_H0$area

plot(difference);plot(areas_H0, add=TRUE)

# Outputs 
writeRaster(difference, "difference_2001.tif", overwrite=TRUE)
# writeOGR(areas_H0, dsn = wd, layer = "areas_H0_photogr_las",driver="ESRI Shapefile")
write.table(areas_H0, "volume_difference_H0_2001.txt")



#### 2009 ####

las <- readLAS("di09_lieg_dom1_20cm.laz", select = "xyz")  #input photogrammetric dense point cloud

###########################################

#Remove outliers SOR
las <- classify_noise(las, sor(8,2))
las <- filter_poi(las, Classification != LASNOISE)
print(las)
# Calculate the DSM without trees #CSF
las <- classify_ground(las, csf(sloop_smooth = TRUE, class_threshold = 3.0, cloth_resolution = 0.5))
# writeLAS(las, "DSM_no_trees_classification.las")
# plot(las_no_trees_CSF, color = "Classification")
las_no_trees_CSF_canopy <- filter_poi(las, Classification == 2)
DSM_no_trees <- grid_metrics(las_no_trees_CSF_canopy, res = res_photog, ~max(Z)); plot(DSM_no_trees)
DSM_no_trees <- focal(DSM_no_trees, w=matrix(1,51,51), fun=mean, NAonly=TRUE, na.rm=TRUE, pad=TRUE); plot(DSM_no_trees)
# DSM_no_trees <- grid_terrain(las, res = res_photog, tin())

#compute CHM
DTM_res <- resample(DTM, DSM_no_trees, method='bilinear');plot(DTM_res)
CHM <- DSM_no_trees - DTM_res;
CHM <- CHM * (CHM > 0); plot(CHM)

#### compute H0
factor_max <- res_max/((res(CHM)[1]+res(CHM)[2])/2)
CHM_max <- aggregate(CHM, fact=factor_max, expand=TRUE, fun=max); res(CHM_max); plot(CHM_max)
CHM_max_res <- resample(CHM_max,CHM, method='bilinear')
CHM_max_res <- CHM_max_res * (CHM_max_res > 0); plot(CHM_max_res)
difference <- abs(CHM_max_res - CHM); plot(difference)

values <- extract(difference, areas_H0, fun=sum)
areas_H0$diff <- as.numeric(values)
areas_H0$volume <- areas_H0$diff*((res(CHM)[1]+res(CHM)[2])/2)^2
areas_H0$area <- area(areas_H0)
areas_H0$H0_volume <- areas_H0$volume / areas_H0$area

plot(difference);plot(areas_H0, add=TRUE)

# Outputs 
writeRaster(difference, "difference_2009.tif", overwrite=TRUE)
# writeOGR(areas_H0, dsn = wd, layer = "areas_H0_photogr_las",driver="ESRI Shapefile")
write.table(areas_H0, "volume_difference_H0_2009.txt")


#### 2019 ####

#las <- readLAS("/home/tommaso/computer_lavoro/phd/valanghe_schianti/aree_studio/Disentis/20190618/20190618_dense_cloud.las", select = "xyz") 
las <- readLAS("20190618_dense_cloud_cl_AOIbuffer.las", select = "xyz")  #input photogrammetric dense point cloud
DTM <- raster("DTM_2019_0p2_2056.tif"); plot(DTM)
areas_H0 <- readOGR("D:/phd/valanghe_schianti/temporal_analysis_disentis/areas", "area_analysis_2056", stringsAsFactors = T);plot(DTM); plot(areas_H0, add=TRUE)

###########################################

#Clip LAS with the 10m buffered shapefile 
#areas_H0_buf <- buffer(areas_H0, width = 5); plot(areas_H0_buf)
#las <- clip_roi(las, areas_H0_buf)
#Remove outliers SOR
las <- classify_noise(las, sor(8,2))
las <- filter_poi(las, Classification != LASNOISE)
print(las)
# Calculate the DSM without trees #CSF
las <- classify_ground(las, csf(sloop_smooth = TRUE, class_threshold = 3.0, cloth_resolution = 0.5))
# writeLAS(las, "DSM_no_trees_classification.las")
# plot(las_no_trees_CSF, color = "Classification")
las_no_trees_CSF_canopy <- filter_poi(las, Classification == 2)
DSM_no_trees <- grid_metrics(las_no_trees_CSF_canopy, res = res_photog, ~max(Z)); plot(DSM_no_trees)
DSM_no_trees <- focal(DSM_no_trees, w=matrix(1,51,51), fun=mean, NAonly=TRUE, na.rm=TRUE, pad=TRUE); plot(DSM_no_trees)
# DSM_no_trees <- grid_terrain(las, res = res_photog, tin())

#compute CHM
DTM_res <- resample(DTM, DSM_no_trees, method='bilinear');plot(DTM_res)
CHM <- DSM_no_trees - DTM_res;
CHM <- CHM * (CHM > 0); plot(CHM)

#### compute H0
factor_max <- res_max/((res(CHM)[1]+res(CHM)[2])/2)
CHM_max <- aggregate(CHM, fact=factor_max, expand=TRUE, fun=max); res(CHM_max); plot(CHM_max)
CHM_max_res <- resample(CHM_max,CHM, method='bilinear')
CHM_max_res <- CHM_max_res * (CHM_max_res > 0); plot(CHM_max_res)
difference <- abs(CHM_max_res - CHM); plot(difference)

values <- extract(difference, areas_H0, fun=sum)
areas_H0$diff <- as.numeric(values)
areas_H0$volume <- areas_H0$diff*((res(CHM)[1]+res(CHM)[2])/2)^2
areas_H0$area <- area(areas_H0)
areas_H0$H0_volume <- areas_H0$volume / areas_H0$area

plot(difference);plot(areas_H0, add=TRUE)

# Outputs 
writeRaster(difference, "difference_2019.tif", overwrite=TRUE)
# writeOGR(areas_H0, dsn = wd, layer = "areas_H0_photogr_las",driver="ESRI Shapefile")
write.table(areas_H0, "volume_difference_H0_2019.txt")

