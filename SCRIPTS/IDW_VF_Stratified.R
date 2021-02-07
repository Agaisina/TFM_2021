#####################################################################################################
#                                                                                                   #
#                                               CAMPIMETERS                                         #
#                                                                                                   #
#                                        (8) (10) IDW STRATIFIED                                    #
#                                                                                                   #
#                           1. Measuring points and (8) (10) cluster visual fields visualization    #
#                           2. Stratified IDW                                                       #
#                              2.1 Compute stratified IDW CV                                        #
#                              2.2 Display error measurements                                       #
#                           3. IDW Examples                                                         #
#                                                                                                   #          
#####################################################################################################

library(lattice)      # Graphics
library(sp)           # Spatial data
library(fs)           # File manipulations
library(tidyverse)    # Data manipulation
library(tools)        # Files manipulations
library(stringr)      # Character manipulation
library(RColorBrewer) # Color Palettes (Plots)
library(fBasics)      # rowSkewness and rowkurtosis
library(robustbase)   # Functions rowMeans, rowMedians
library(rgl)          # 3D visualization
library(dplyr)        # Function mutate
library(factoextra)   # fviz_nbclust, PCA
library(NbClust)      # Clustering
library(rgdal)        # Read sph
library(glue)         # Equivalent of python 'format' (str)
library(gstat)        # Geostatistics
library(RcmdrMisc)    # Function Hist()
library(spatialEco)   # sp.na.omit()
library(gridExtra)     # Grid arrange
library(raster)
library(rgeos)
library(dismo)

#########################################################################################################
#     0. Load data and define number of visual field clusters                                           #
#########################################################################################################

n_cluster <- 8 # Number of visual field areas:
# Visual field area shapefile:
eye_layout <- glue("Eye_clusters{n_cluster}") # Zonification VF Shapefile 
sn = 1 # Sample num 

# Load data:
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
data_clust_folder <- glue("Data_{n_cluster}Cluster")
data_sample_folder <- glue("5260_Patients_Sample{sn}")
path <- paste0(root_path,"/",data_clust_folder,"/",data_sample_folder)
setwd(path)
file_names <- file_path_sans_ext(list.files(file.path(getwd()), pattern = "*.csv"))

#########################################################################################################
#     1. Measuring points and (8) (10) cluster visual fields visualization                              #
#########################################################################################################

# Dummy campimeter to extract points
path_dummy <- file.path(file.path(getwd()),paste0(file_names[1],".csv"))
dummy <- read.csv2(path_dummy, header = TRUE, sep = ";", skipNul = FALSE)
dummy[is.na(dummy)] <- 0 #Fill NA with 0 Cluster (Center point)
coordinates <- dummy[,1:2]
dummy_sp <- SpatialPointsDataFrame(coordinates, dummy, coords.nrs = numeric(0), 
                                      proj4string = CRS(as.character(NA)), bbox = NULL)
# Visualize measuring points map
x <- dummy$x
y <- dummy$y
map <- xyplot(y~x)
labels <- as.character(row.names(dummy))
myColours <- brewer.pal(6,"Blues") # Blue palette
xyplot(y~x,pch=16, col = myColours[3], main=list(label="Puntos de medida", cex=1))
trellis.focus("panel",1,1)
panel.text(x=map$panel.args[[1]]$x,y=map$panel.args[[1]]$y,labels = labels, pos=3, cex=0.6, col=myColours[6])
trellis.unfocus()
# Visualize internal clusters
dev.off()
path_eye <- "C:/Users/Agata/MASTER/TFM/Eye_Layers/Eye_cluster"
eye_clusters <- readOGR(path_eye, eye_layout) #Eye cluster shapefile 
dummy_sp$Cluster <- as.factor(dummy_sp$Cluster)
levels(dummy_sp$Cluster)
plot_colors <- brewer.pal(n_cluster,"Paired") # Paired palette
plot_colors <- c('dark gray', plot_colors) # Uncomment for 8 and 10 cluster (Center point not include)
plot(eye_clusters,border='gray',lwd = 2,main = "Zonificación del campo visual")
plot(dummy_sp,
     col = (plot_colors)[dummy_sp$Cluster],
     pch = 19,
     add=T)
par(xpd=TRUE) # Allow writting out of plot bounds
legend(6,5.45,
       legend = c(levels(dummy_sp$Cluster)),
       col = plot_colors, # set the color of each legend line
       pch=19,
       bty = "y", # turn off border,
       cex = .9) # adjust legend font size

#########################################################################################################
#                                        2. STRATIFIED IDW                                              #
#########################################################################################################


# 2.1 Compute Stratified IDW CV results 

# Function STRATIFIED CV IDW 
# Only errors
# Kigriging errors with variogram fitted to model residuals

idw.results.strat <- function(file_names,path,n_strat,nmax,power){
        
        # Vectores con todos los valores del conjunto
        x <- c()
        y <- c()
        predicted <- vector()
        observed <- vector()
        residuals <- vector()
        name_camp <- vector()
        strat_name <- vector()
        n_neighbours <- vector()
        # Vectores para guardar medias y desviaciones del resultado de cada campimetria
        # No de todo el conjunto de test
        mean_res <- vector()
        var_res <- vector()
        mspe <- vector()
        rmspe <- vector()
        cor_obs_pred <- vector()
        cor_pred_res <- vector()
        name_list <- vector()
        neighbours_list <- vector()
        power_list <- vector()
        
        
        for (i in file_names){
                filepath <- file.path(file.path(paste0(path,"/",paste0(format(i),".csv"))))
                name <- file_path_sans_ext(basename(filepath))
                camp <- read.csv2(filepath, header = TRUE, sep = ";", skipNul = TRUE)
                coordenadas <- camp[,1:2]
                camp_sp <- SpatialPointsDataFrame(coordenadas, camp, coords.nrs = numeric(0), 
                                                  proj4string = CRS(as.character(NA)), bbox = NULL)
                camp_sp <- sp.na.omit(camp_sp, margin = 1) # Remove center point (does not belong to any cluster)
                
                # Vectores para guardar medias del resultado de cada estrato de cada campimetria 
                # (Ejm: error medio del estrato 1 en la campimetría 1. El vector se llena con los resultados de todos los estratos)
                mean_res.s <- vector()
                var_res.s <- vector()
                mspe.s <- vector()
                rmspe.s <- vector()
                cor_obs_pred.s <- vector()
                cor_pred_res.s <- vector()
                name_list.s <- vector()
                neighbours_list.s <- vector()
                
                for (i in (1:n_strat)){
                        strat <- i
                        data_strat <- camp_sp[camp_sp$Cluster==i,]
                        n_point_strat <- dim(data_strat)[1]
                        coordenadas.x <- data_strat$x
                        coordenadas.y <- data_strat$y
                        idw_camp <- gstat(formula = VNorm ~ 1, #intercept only model
                                          locations=data_strat, 
                                          nmax=nmax,
                                          set=list(idp = power))
                        model.pred.idw <- gstat.cv(idw_camp) #Default LOOCV
                        
                        #Results of IDW
                        pred.idw.s <- model.pred.idw$var1.pred
                        observed.idw.s <- model.pred.idw$observed
                        residuals.idw.s <- model.pred.idw$residual
                        #Means and variances (For each campimetry)
                        mean_res.s <- c(mean_res.s,mean(residuals.idw.s))
                        var_res.s <- c(var_res.s,var(residuals.idw.s))
                        #Error measures (For each campimetry)
                        mspe.s <- c(mspe.s,mean(residuals.idw.s^2))
                        rmspe.s <- c(rmspe.s,sqrt(mean(residuals.idw.s^2)))
                        #correlation observed and predicted, ideally 1
                        cor_obs_pred.s <- c(cor_obs_pred.s, abs((cor(observed.idw.s, pred.idw.s))))
                        #correlation predicted and residual, ideally 0
                        cor_pred_res.s <- c(cor_pred_res.s, abs((cor(pred.idw.s, residuals.idw.s))))
                        
                        #Vector with ALL results 
                        name_camp <- c(name_camp,rep(name,n_point_strat))
                        strat_name <- c(strat_name,rep(strat,n_point_strat))
                        x <- c(x,coordenadas.x)
                        y <- c(y,coordenadas.y)
                        n_neighbours <- c(n_neighbours, rep(nmax,n_point_strat))
                        predicted <- c(predicted, pred.idw.s)
                        observed <- c(observed, observed.idw.s)
                        residuals <- c(residuals, residuals.idw.s)
                        print(length(residuals))
        
                        # Medias de los errores por estrato:
                        # mean_res <- c(mean_res,mean(mean_res.s))
                        # var_res <- c(var_res,mean(var_res.s))
                        # mspe <- c(mspe,mean(mspe.s))
                        # rmspe <- c(rmspe,mean(rmspe.s))
                        # cor_obs_pred <- c(cor_obs_pred,mean(cor_obs_pred.s))
                        # cor_pred_res <- c(cor_pred_res,mean(cor_pred_res.s))
                        # name_list <- c(name_list,name)
                        # neighbours_list <- c(neighbours_list,nmax)
                        # power_list <- c(power_list,power)
                        
                }
                # Media de los errores por campimetría
                mean_res <- c(mean_res,mean(mean_res.s))
                var_res <- c(var_res,mean(var_res.s))
                mspe <- c(mspe,mean(mspe.s))
                rmspe <- c(rmspe,mean(rmspe.s))
                cor_obs_pred <- c(cor_obs_pred,mean(cor_obs_pred.s))
                cor_pred_res <- c(cor_pred_res,mean(cor_pred_res.s))
                name_list <- c(name_list,name)
                neighbours_list <- c(neighbours_list,nmax)
                power_list <- c(power_list,power)
                        
        }
        # Devolver resultados kriging
        all.result <- data.frame(name_camp,x,y,strat_name,n_neighbours,observed,predicted,residuals)
        all.errors <- data.frame(name_list,neighbours_list,mean_res,var_res,mspe,rmspe,cor_obs_pred,cor_pred_res)
        return <- list(all.result,all.errors)
        #return <- all.errors
        return(return)
}


# Read sample 4
ns <- 4
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
clust_folder <- glue("Data_{n_cluster}Cluster")
path_sample <- file.path(paste0(root_path,"/",clust_folder,"/",glue("5260_Patients_Sample{ns}")))
file_names <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))
file_names_test <- file_names # Sample 4
n_strat <- n_cluster

idw_test <- idw.results.strat(file_names_test,path_sample, n_strat, 5,2)
write.csv2(idw_test[[1]], file.path(root_path,paste0("idw_test_results_VF10_strat.csv")))
write.csv2(idw_test[[2]], file.path(root_path,paste0("idw_test_errors_VF10_strat.csv")))


# 2.2 Display error measurements
n_cluster = 8
n_clust_p = 12

path_results <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/RESULTS"
results_folder <- glue("IDW_stratified")
path_results <- file.path(path_results,results_folder)
errors_test <- read.csv2(file.path(path_results,glue("idw_test_errors_VF{n_cluster}_strat.csv")), sep = ";")
results_test <-read.csv2(file.path(path_results,glue("idw_test_results_VF{n_cluster}_strat.csv")), sep = ";")

# Mean of all test campimeters errors by cluster

# Mean of all test campimeters errors 

errors_test_mean <- errors_test %>% group_by(neighbours_list) %>%
        summarize_all(mean,na.rm = TRUE)
#View(errors_test_mean)
#View(errors_test)

# Residual histogram

blues <- brewer.pal(5,"Blues")
breaks_me <- seq(-11,15,0.5)
breaks_mspe <- seq(0,5000,100)
#breaks_mspe_n <- seq(0,8,0.5)
breaks_cor_obs_pred <- seq(0,1,0.05)
breaks_var <- seq(0,5000,150)
#breaks_cor_pred_zscore <- seq(-0.5,1,0.1)
#breaks_percent <- seq(0,15,0.5)

par(mfrow=c(1,3),oma=c(0,0,2,0))
par(mar=c(3,2,2.4,1))
Hist(errors_test$mean_res, main="Errores medios",ylim=c(0,3000),cex.main=0.7,
     ylab="Porcentaje",xlab=NULL, breaks=breaks_me, col=blues[4], scale="percent")
Hist(errors_test$mspe, main="Errores cuadráticos medios",ylim=c(0,3000),cex.main=0.7,
     ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe, col=blues[4], scale="percent")
Hist(errors_test$cor_obs_pred, main="Correlación observado-predicción",ylim=c(0,3000),cex.main=0.7,
     ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_obs_pred, col=blues[4], scale="percent")
#Hist(errors_test$var_res, main="Varianza entre errores",ylim=c(0,3000),cex.main=0.7,
     #ylab="Porcentaje",xlab=NULL, breaks=breaks_var, col=blues[4], scale="percent")
title(glue("Modelo IDW Estratificado (Zonificación en {n_cluster} sectores)"), line = 0.5, outer = TRUE, cex.main=1.3)
dev.off()

# 2.3 Spplot mean residuals

results_test_mean <- results_test %>% group_by(x,y) %>%
        summarize_all(mean,na.rm = TRUE)
View(results_test_mean)

cuts <- c(-10,-5,-2.5,-1,-0.5,0.5,1,2.5,5,10)
cols <- brewer.pal(length(cuts)/2, "OrRd")
cols1 <- rev(cols)
cols2 <- cols
cols <- c(cols1,cols2)

coordinates <- results_test_mean[,1:2]
results_test_mean_sp <- SpatialPointsDataFrame(coordinates, results_test_mean, coords.nrs = numeric(0), 
                                   proj4string = CRS(as.character(NA)), bbox = NULL)
eye_plot(results_test_mean_sp,"residuals","Residuos medios",cols,cuts)



################################################################################################
#                                        3. IDW EXAMPLES                                       #
################################################################################################

n_cluster <- 10 # Number of visual field areas

# Eye circle shapefiles
path_eye <- "C:/Users/Agata/MASTER/TFM/Eye_Layers/Eye_cluster"
eye_cir <- readOGR(path_eye, "Eye_Cir") #Eye circular shapefile 
eye_VF <- readOGR(path_eye, glue("Eye_clusters{n_cluster}")) #Eye circular + VF shapefile 
eye_bs <- readOGR(path_eye, "Eye_Blind") #Eye circular + Blind Spot shapefile 
blindSpot <- readOGR(path_eye, "BlindSpot") #Blind Spot shapefile 
grid <- raster(extent(eye_cir))
res(grid) <- 0.5 # Resolution
grid <- as(grid, "SpatialPixelsDataFrame") #Pixel Data Frame
eye_pixels <- grid
#plot(eye_pixels)

# Eye zones shapefiles
path_eye <- "C:/Users/Agata/MASTER/TFM/Eye_Layers/Eye_cluster"
VF_folder <- glue("VF{n_cluster}_Parts")
path_eye_zones <- file.path(path_eye,VF_folder)
zones <- list.files(path_eye_zones)
zones_names_list <- paste0(glue("VF{n_cluster}_"),1:n_cluster)
n <- 1
for (i in zones_names_list){
        zone_layout <- readOGR(path_eye_zones,i)
        assign(glue("zone_{n}"),zone_layout)
        n <- n+1
}

# Interpolate map

sn = 4
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
data_clust_folder <- glue("Data_{n_cluster}Cluster")
data_sample_folder <- glue("5260_Patients_Sample{sn}")
path_sample <- paste0(root_path,"/",data_clust_folder,"/",data_sample_folder)
file_names_sample <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))
ejm <- #camps_test_c1[1] #"1036210_camp3"
s <- "1036210_camp3"
ejm_idx <- match(s,file_names_sample)
file_names_ejms <- file_names_sample[ejm_idx]


power=2
nmax=5
path <- path_sample
n_strat <- n_cluster
# Plot parameters 
scale.parameter = 1.1
original.bbox <- eye_bs@bbox
xshift = 0  # Shift to right in map units. 
yshift = 0  # Shift to left in map units.
edges = original.bbox
edges[1,] <- (edges[1,] - mean(edges[1,])) * scale.parameter + mean(edges[1,]) + xshift
edges[2,] <- (edges[2,] - mean(edges[2,])) * scale.parameter + mean(edges[2,]) + yshift
for (i in file_names_ejms){
        filepath <- file.path(file.path(paste0(path,"/",paste0(format(i),".csv"))))
        name <- file_path_sans_ext(basename(filepath))
        camp <- read.csv2(filepath, header = TRUE, sep = ";", skipNul = TRUE)
        coordenadas <- camp[,1:2]
        camp_sp <- SpatialPointsDataFrame(coordenadas, camp, coords.nrs = numeric(0), 
                                          proj4string = CRS(as.character(NA)), bbox = NULL)
        puntos_pred <- SpatialPointsDataFrame(coordenadas, camp, coords.nrs = numeric(0), 
                                              proj4string = CRS(as.character(NA)), bbox = NULL)
        camp_sp <- sp.na.omit(camp_sp, margin = 1) # Remove center point (does not belong to any cluster)
        all_coords <- data.frame()
        all_data <- data.frame()
        all_maps <- c()
        for (i in (1:n_strat)){
                strat <- i
                data_strat <- camp_sp[camp_sp$Cluster==i,]
                grid <- raster(extent(eye_cir))
                res(grid) <- 0.5 # Resolution
                grid <- as(grid, "SpatialPixelsDataFrame") #Pixel Data Frame
                eye_pixels <- grid
                idw_map <- idw(formula = VNorm ~ 1, #intercept only model
                               locations=data_strat, eye_pixels,
                                nmax=nmax,idp = power)
                idw_map_strat <- crop(idw_map, get(glue("zone_{strat}"))) # Clip map
                all_coords <- rbind(all_coords, idw_map_strat@coords)
                all_data <- rbind(all_data, idw_map_strat@data)
        }
        all_maps <- cbind(all_coords,all_data)
        rast <- raster()
        extent(rast) <- extent(eye_cir) # this might be unnecessary
        ncol(rast) <- 21 # this is one way of assigning cell size / resolution
        nrow(rast) <- 21
        mapa_ras <- rasterize(all_maps[,1:2], rast, all_maps$var1.pred, fun=mean)
        spplot(mapa_ras, cex = 0.7, do.log = F,
                       main=list(label="Defecto normalizado",cex=1.2),
                       sp.layout = list(blindSpot, fill='darkgray',first = FALSE),
                       xlim = edges[1,], ylim = edges[2,],
                       key.space="right")
}

grid.arrange(spplot(puntos_pred, "VNorm", cex = 0.7, do.log = F,
                    main=list(label="Defecto normalizado",cex=1.2),
                    #col.regions=cols,cuts=cuts,
                    sp.layout = list("sp.polygons",eye_bs,col="dark gray"),
                    xlim = edges[1,], ylim = edges[2,],
                    key.space="right"),
             spplot(mapa_ras, cex = 0.7, do.log = F,
                    main=list(label="Defecto normalizado",cex=1.2),
                    sp.layout = list(blindSpot, fill='darkgray',first = FALSE),
                    xlim = edges[1,], ylim = edges[2,],
                    key.space="right"),ncol=2, widths = c(2,1.66),
             top=glue("IDW estratificado: {n_cluster} Zonas"))

# spplot(mapa_ras, cex = 0.7, do.log = F,
#        main=list(label=glue("Defecto normalizado"),cex=1.2),
#        sp.layout = list(list(eye_VF, col="lightblue",first=FALSE),list(blindSpot, fill='darkgray',first = FALSE)),
#        xlim = edges[1,], ylim = edges[2,],
#        key.space="right")
