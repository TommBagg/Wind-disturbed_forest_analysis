# Title: Tree identification and crown extraction in windthrow areas 10/12/2021
# Author details: Tommaso Baggio Contact details: tommaso.baggio@phd.unipd.it; tbaggio93@gmail.com

# Copyright statement: This script is the product of the work of Tommaso Baggio - 2021

# THE SCRIPT IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Necessary inputs: point cloud (.las format) from LiDAR or photogrammetric survey
#                  dtm representing the ground level
#                  shapefile reporting the areas of interest for which the algorithm computes the analysis  

## tree top identification -> Popescu & Wynne (2004)
## crown area identification -> Silva 2016

library(raster)
library(rgdal)
library(lidR)
library(rgeos)

#Set the path of the working directory where input data are stored
wd <- ("/home/tommaso/computer_lavoro/phd/valanghe_schianti/temporal_analysis_disentis/crown_analysis"); setwd(wd)

#Input data
las <- readLAS("di91_lieg_dom1.laz", select = "xyz")
dtm <- raster("dtm.tif")
areas_H0 <- readOGR("/home/tommaso/computer_lavoro/phd/valanghe_schianti/temporal_analysis_disentis/areas", "area_analysis_21781", stringsAsFactors = T); plot(areas_H0)

# Parameters for tree extraction and analysis
hmin_tree = 3.75 #minimum height to identify a tree
hmin_neigh = 3.75 #minimum height to identify a tree in the neighborhood procedure
cmin_tree = 2 #minimum crown area to identify a tree
width <- 4 #width of the buffer around the identified crown contours to calculate the corrected height


##############################
# START OF THE ALGORITHM

las <- normalize_height(las,dtm)
las <- classify_noise(las, sor(10,3))
las <- filter_poi(las, Classification != LASNOISE)
chm_las <- grid_canopy(las = las, res = 0.25, p2r(subcircle = 0.25,na.fill = tin()))
chm_las_mask0 <- chm_las > 0
chm_las <-chm_las * chm_las_mask0
plot(chm_las, col = height.colors(50), zlim=c(0,10)); plot(areas_H0, add=TRUE)
ttops <- find_trees(las, lmf(7,hmin = hmin_tree, shape = "circular")); plot(ttops, add = T) #7x7 m fixed window size
#x <- plot(las); add_treetops3d(x, ttops)
crowns = silva2016(chm = chm_las, treetops = ttops, exclusion = 0.3, max_cr_factor = 0.5)()
contour = rasterToPolygons(crowns, dissolve = TRUE)
plot(chm_las, col = height.colors(50)); plot(contour, add = T); plot(areas_H0, add=TRUE)
contour$area <- area(contour)

contour <- subset(contour, area>cmin_tree)

###### difference between the tree top and the mean height of the vegetation on the ground in a neighborhood window
plot(chm_las, col = height.colors(50)); plot(contour, add = T); plot(areas_H0, add=TRUE)
contour$neigh_height <- 0.00

for (i in 1:length(contour)) {
  print(i)
  cont <- contour[1,]
  cont_buf <- buffer(cont, width = width, dissolve=F)
  cont_rest <- contour - cont; plot(cont_rest, add=T, col="green")
  cont_bufcl <- cont_buf - cont_rest; plot(cont_bufcl, add=T, col="white");plot(cont, add=T, col="red")
  n_cont <- as.numeric(nrow(as.data.frame(extract(chm_las, cont))))
  n_buf <- as.numeric(nrow(as.data.frame(extract(chm_las, cont_bufcl))))
  s_cont <- as.numeric(extract(chm_las, cont, fun=sum))
  s_buf <- as.numeric(extract(chm_las, cont_bufcl, fun=sum))
  contour[i,] <- (s_buf-s_cont)/(n_buf-n_cont)
}
###### 

contour$max_chm <- extract(chm_las, contour, fun=max); contour$max <- as.numeric(contour$max_chm); contour$max_chm <- NULL
contour$dif_height <- contour$max - contour$neigh_height

contour_fin <- subset(contour, dif_height > hmin_neigh)
plot(contour_fin, add=TRUE, col="pink")
plot(chm_las, col = height.colors(50)); plot(contour_fin, add = T); plot(areas_H0, add=TRUE)

tree_p <- gCentroid(contour_fin, byid = T)
tree_p <- SpatialPointsDataFrame(tree_p, data= contour_fin@data)
plot(chm_las, col = height.colors(50)); plot(contour_fin, add = TRUE); plot(areas_H0, add=TRUE); plot(tree_p, add = TRUE)

writeOGR(contour_fin, dsn = wd, layer = "tree_crown", driver = "ESRI Shapefile", overwrite_layer=TRUE )
writeOGR(tree_p, dsn = wd, layer = "tree_top", driver = "ESRI Shapefile", overwrite_layer=TRUE )


