#####################################################################################################
#                                                                                                   #
#                                               CAMPIMETERS                                         #
#                                                                                                   #
#                                           "GLOBAL" VARIOGRAMS                                     #
#                                                                                                   #
#                                       0. Read and display tarining dataset                        #
#                                       1. "Global variogram" VNorm ~ 1                             #                                           #
#                                         1.1 Variogram function distance                           #                                                      #
#                                           1.1.1 Distribution of semivariance values               #
#                                         1.2 Variogram of means and medians                        #
#                                       2. Fit variogram models                                     #                      #
#                                         2.1 Means Global Model                                    #                    #
#                                         2.2 Medians Global Model                                  #     
#                                       3. Kriging Cross-Validation  (Validation Data set)          #             
#                                       4. Kriging Cross-Validation  (Test Data set)                #
#                                       5. Kriging examples
#
#####################################################################################################

library(lattice)       # Graphics
library(RColorBrewer)  # Color Palettes for plots
library(tools)         # Files manipulation (file_path_sans_ext function)
library(sp)            # Spatial data
library(gstat)         # Geostatistics
library(fBasics)       # rowSkewness and rowkurtosis
library(robustbase)    # Functions rowMeans, rowMedians
library(ggplot2)       # Plotting
library(resample)      # Function colVar, colStdevs, colMeans
library(glue)          # Format str
library(rgdal)         # Read sph
library(gridExtra)     # Grid arrange
library(raster)
library(rgeos)
library(dismo)
#library(fs)           # File manipulations
#library(tidyverse)    # Data manipulation
#library(stringr)      # Character manipulation
#library(stats)        # Basic statistics
#library(geoR)         # GeoR


#########################################################################################################
#                                            0. READ AND DISPLAY TRAINING DATASET                                   #
#########################################################################################################

# Read data
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Data_IDW"
setwd(root_path)
# Read sample 2 
ns <- 2
path_sample <- file.path(paste0(root_path,"/",glue("5260_Patients_Sample{ns}")))
file_names_train <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))
file_names_train_500 <- sample(file_names_train,500)

# Dummy campimeter
# Coordinates points -> Read a dummy campimeter and extract it
path <- path_sample
setwd(path)   # Fixed working directory
file_names <- list.files(getwd(), pattern = "*.csv")
dummy_camp <- read.csv2(file.path(file.path(getwd(),paste(file_names[1],sep=""))))
x <- dummy_camp$x
y <- dummy_camp$y
#plot(variogram(VNorm ~ 1, dummy_sp,boundaries=c(1,2,3,4,5,6,7,8,9)))

###############################################################################################
# # Visualizar media y desviación de los datos de entrenamiento
# 
# cols <- brewer.pal(10, "Spectral")
# cols2 <- brewer.pal(6, "Spectral")
# cols3 <- rev(cols)
# cuts <- c(40,0,-20,-40,-80,-140,-180,-220,-260,-300)
# cuts_skw <- c(-3,-1.5,-0.5,0.5,1.5,3)
# cuts_kurt <- c(-15,-5,-1,1,5,15)
# cuts_sd <- c(15,30,45,60,75,100)
# 
# # Read data frame that contains all values of the random sample (5260 campimeters) -> Previosly created
# # Calculate statistics by rows
# sn = 2
# data_values_train <- read.csv2(paste0(root_path,"/5260_Values_df",
#                                      glue("/df_5260_Values_campimeters{sn}.csv")))
# data <- data_values_train
# rownames(data) <- data[,1]
# data[,1] <- NULL
# data_train = data[row.names(data) %in% file_names_train,]
# data_train_t <- as.data.frame(t(data_train)); data_train_t <- cbind(y = y, data_train_t); data_train_t <- cbind(x = x, data_train_t)
# coordenadas <- data_train_t[,1:2]
# data_train_t$VNorm_mean<- rowMeans(data_train_t[c(-1,-2)])
# data_train_t$VNorm_median<- rowMedians(as.matrix(data_train_t[c(-1,-2)]))
# data_train_t$VNorm_skew<- rowSkewness(as.matrix(data_train_t[c(-1,-2)]))
# data_train_t$VNorm_kurt<- rowKurtosis(as.matrix(data_train_t[c(-1,-2)]))
# data_train_t$VNorm_sd<- rowStdevs(as.matrix(data_train_t[c(-1,-2)]))  
# data_train_mean_sd <- cbind(coordenadas,data_train_t$VNorm_mean,data_train_t$VNorm_sd)
# colnames(data_train_mean_sd) <- c("x","y","Train_mean","Train_sd")
# data_train_mean_sd <- SpatialPointsDataFrame(coordenadas, data_train_mean_sd, coords.nrs = numeric(0), 
#                                    proj4string = CRS(as.character(NA)), bbox = NULL)
# spplot(data_train_mean_sd, c("Train_mean","Train_sd"),names.attr=c("Media","Desviación Típica"),cuts=c(cuts,cuts_sd),
#        col.regions=c(cols,cols3),key.space="left",
#        main=list(label="Datos de entrenamiento - Defecto visual",cex=1))
# dev.off()

#########################################################################################################
#                                          1. GLOBAL VARIOGRAM                                          #
#########################################################################################################

##############################       1.1  VARIOGRAM FUNCTION DISTANCE        ############################

# Calculate all variograms and put the values together. 
# Each campimetry has a semi-variance value in each lag

file_names_train <- file_names_train_500
variograms_values <- vector()

for(i in file_names_train){
  filepath <- file.path(file.path(getwd()),paste0(format(i),".csv"))
  camp <- read.csv2(filepath, header = TRUE, sep = ";", skipNul = TRUE)
  coordinates = camp[,1:2]
  camp_sp <- SpatialPointsDataFrame(coordinates, camp, coords.nrs = numeric(0), 
                                        proj4string = CRS("+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"), 
                                        bbox = NULL)
  v_i <- variogram(VNorm ~ 1 , camp_sp, boundaries=c(1,2,3,4,5,6,7,8,9)) #Boundaries
  variograms_values <- rbind(variograms_values,as.vector(v_i[,3]))
}
row.names(variograms_values) <- file_names_train

# Plot results and identify outliers
dist_vector <- rep(as.vector(v_i[,2]),length(file_names_train))
plot(dist_vector, variograms_values, xlab="Distance (mm)", ylab="Semivariance", main="Train semivariance values")
#identify(dist_vector, variograms_values, labels=row.names(variograms_values), plot=TRUE) # identify points

# Outliers visualization
#cuts <- c(40,0,-20,-40,-80,-140,-180,-220,-260,-300)
#colors <- brewer.pal(10, "Spectral")
#display.brewer.pal(10,"Spectral")
#outlier_name <- "837066_camp3"
#outlier <- read.csv2(file.path(file.path(getwd()),paste0(format(outlier_name),".csv")))
#View(outlier)
#coords <- outlier[,1:2]
#outlier_sp <- SpatialPointsDataFrame(coords, outlier, coords.nrs = numeric(0), 
                                      #proj4string = CRS("+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"), bbox = NULL)
#spplot(outlier_sp, "VNorm", do.log = F, main=list(label=paste0("VNorm ",outlier_name),cex=1.5),
       #cuts=cuts, cols=colors, key.space="right")


##########################       1.1.1  DISTRIBUTION OF SEMI-VARIANCE VALUES        ##########################

# Boxplot
variograms_values_train <- as.data.frame(variograms_values)
variogram_boxplot <- boxplot(variograms_values_train$V1, variograms_values_train$V2, variograms_values_train$V3, 
                             variograms_values_train$V4, variograms_values_train$V5, variograms_values_train$V6,
                             variograms_values_train$V7, variograms_values_train$V8, variograms_values_train$V9,
                             main = "Boxplots semivariance values",
                             at = c(1,2,3,4,5,6,7,8,9),
                             names = c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9"),
                             las = 2,
                             col = c("light blue","yellow"),
                             border = "dark blue",
                             horizontal = FALSE,
                             notch = TRUE
)

# Calculate upper whiskers
upper_whiskers <- variogram_boxplot$stats[5,]
upper_whiskers

# View table of values
#View(variograms_values_train)

# Display one campimeter
# cuts <- c(-50,10,50,100,150,200,300)
# cols <- brewer.pal(length(cuts), "Spectral")
# cols <- rev(cols)
# 
# camp_name <- "689762_camp2"
# camp_one <- read.csv2(file.path(file.path(getwd()),paste0(format(camp_name),".csv")))
# #View(outlier)
# coordenadas <- outlier[,1:2]
# camp_one_sp <- SpatialPointsDataFrame(coordenadas, camp_one, coords.nrs = numeric(0), 
#                                      proj4string = CRS("+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"), bbox = NULL)
# spplot(camp_one_sp, "VNorm", do.log = F, main=list(label=paste0("VNorm ",camp_name),cex=1.5),
#        cuts=cuts,col.regions=cols, key.space="right")


############################      1.2   VARIOGRAM OF MEANS AND MEDIANS         ##########################

# Calcular Varianza, desviación típica, media y mediana
# de los valores en cada lag de distancia. Un punto para cada lag

variograms_values <- vector()

for(i in file_names_train){
  filepath <- file.path(file.path(getwd()),paste0(format(i),".csv"))
  camp <- read.csv2(filepath, header = TRUE, sep = ";", skipNul = TRUE)
  coordinates(camp) = ~ x + y
  v_i <- variogram(VNorm ~ 1 , camp, boundaries=c(1,2,3,4,5,6,7,8,9))
  variograms_values <- rbind(variograms_values,as.vector(v_i[,3]))
  variance <- colVars(variograms_values)
  std <- colStdevs(variograms_values)
  medias <- colMeans(variograms_values)
  medianas <- colMedians(variograms_values)
  sum <- colSums(variograms_values)
}
myColours <- brewer.pal(6,"Blues") # Paleta Azules
plot(as.vector(v_i[,2]), medias, xlab="Distancia (mm)", ylab = "Media de los valores de semivarianza",
     col=myColours[5], pch=16, main = "Modelo global de medias")
plot(as.vector(v_i[,2]), medianas, xlab="Distancia (mm)", ylab = "Mediana de los valores de semivarianza",
     col=myColours[5], pch=16, main="Modelo global de medianas")


#########################################################################################################
#                                     2. FIT VARIOGRAMS MODELS                                          #
#########################################################################################################

# Variograma prueba (auxiliar)
# Dummy campimeter
dummy <- read.csv2(file.path(file.path(getwd()),paste0(file_names_train[1],".csv")), 
                   header = TRUE, sep = ";", skipNul = TRUE)
coords <- dummy[,1:2]
dummy_sp <- SpatialPointsDataFrame(coords, dummy, coords.nrs = numeric(0), 
                                   proj4string = CRS("+proj=utm +zone=30 +ellps=GRS80 +units=m +no_defs"), 
                                   bbox = NULL)
plot(variogram(VNorm ~ 1, dummy_sp,boundaries=c(1,2,3,4,5,6,7,8,9)))
variogram_prueba <- variogram(VNorm ~ x + y, dummy_sp, boundaries=c(1,2,3,4,5,6,7,8,9))


coordinates(dummy) = ~ x + y 
variogram_prueba <- variogram(VNorm ~ 1, dummy, boundaries=c(1,2,3,4,5,6,7,8,9))
np <- variogram_prueba$np # Número pares puntos
dist <- variogram_prueba$dist # Dist breaks

# Variance
#var.vgm = data.frame(dist=dist_vector,gamma=variance)
#var.vgm$np = np #Number of point pairs (Calculado antes)
#class(var.vgm) = c("gstatVariogram","data.frame")
#vg.exp <- vgm(psill=1,model='Exp', range = 3)
#fit.var.exp <- fit.variogram(var.vgm, model=vg.exp,fit.method=1)

# Standard deviation
#std.vgm = data.frame(dist=dist_vector,gamma=std)
#std.vgm$np = np #Number of point pairs (Calculado antes)
#class(std.vgm) = c("gstatVariogram","data.frame")
#vg.exp <- vgm(psill=1,model='Exp', range = 3)
#fit.std.exp <- fit.variogram(std.vgm, model=vg.exp,fit.method=1)

####################################################################################################
#                                       2.1 MEANS GLOBAL MODEL                                     #
####################################################################################################
# http://r-spatial.github.io/gstat/reference/fit.StVariogram.html
# Default method fit = 7. This criterion is not supported by theory, but by practice

media.vgm = data.frame(dist=dist,gamma=medias)
media.vgm$np = np #Number of point pairs (Calculado antes)
media.vgm$dir.hor = rep(0,length(np))
media.vgm$dir.ver = rep(0,length(np))
class(media.vgm) = c("gstatVariogram","data.frame")

# Exponencial
vg.exp <- vgm(psill=max(media.vgm$gamma)*0.9, model = "Exp", range=max(media.vgm$dist)/2, 
              nugget = mean(media.vgm$gamma)/4)
fit.media.exp <- fit.variogram(media.vgm, vg.exp)
plot(media.vgm,fit.media.exp, main = "Mean semivariogram", sub= "Exponential model")

# Esférico
vg.sph <- vgm(psill=max(media.vgm$gamma)*0.9, model = "Sph", range=max(media.vgm$dist)/2, 
              nugget = mean(media.vgm$gamma)/4)
fit.media.sph <- fit.variogram(media.vgm, vg.sph)
plot(media.vgm,fit.media.sph, col = "green",main = "Mean semivariogram", sub= "Spherical model")

# Gaussiano
vg.gauss <- vgm(psill=max(media.vgm$gamma)*0.9, model = "Gau", range=max(media.vgm$dist)/3, 
                nugget = mean(media.vgm$gamma)/4)
fit.media.gauss <- fit.variogram(media.vgm, vg.gauss)
plot(media.vgm,fit.media.gauss, col = "red",main = "Mean semivariogram", sub= "Gaussian model")

# Lineal
vg.lin <- vgm(psill=max(media.vgm$gamma)*0.9, model = "Lin", range=max(media.vgm$dist)/2, 
              nugget = mean(media.vgm$gamma)/4)
fit.media.lin <- fit.variogram(media.vgm, vg.lin)
plot(media.vgm,fit.media.lin, col = "purple",main = "Mean semivariogram", sub= "Lineal model")

# Circular
vg.cir <- vgm(psill=max(media.vgm$gamma)*0.9, model = "Cir", range=max(media.vgm$dist)/2, 
              nugget = mean(media.vgm$gamma)/4)
fit.media.cir <- fit.variogram(media.vgm, vg.cir)
plot(media.vgm,fit.media.cir, col = "dark grey",main = "Mean semivariogram", sub= "Circular model")


####################################### Plot all fitted models ###########################################

Colores <- brewer.pal(10,"Paired") # Paleta Paired
list_col <- c(Colores[1],Colores[2],Colores[3],Colores[7],Colores[5]) #Todos
models_name <- c("Exponential","Spheric","Gaussian", "Linear", "Circular") #Todos
models_name <- c("Exponencial","Esférico","Gausiano", "Lineal", "Circular") #Todos

acc_par <- list(superpose.line = list(col = list_col, lwd=1.5, lty=1.5))

list_col <- list_col
models_name <- models_name

xyplot(gamma ~ dist, media.vgm, xlim=c(0,8), ylim=c(500,2800), pch=20, col= "black", 
       xlab="Distancia (mm)", ylab="semivarianza",main="Modelo global de medias",
       panel = function(...) {
         panel.xyplot(...)
         
         # First model: exponential
         ret = variogramLine(fit.media.exp, maxdist = 8)
         llines(ret$dist, ret$gamma, col=Colores[1])
         
         # Second model: spheric
         ret = variogramLine(fit.media.sph, maxdist = 8)
         llines(ret$dist, ret$gamma, col=Colores[2])
         
         # Third model: gaussian
         ret = variogramLine(fit.media.gauss, maxdist = 8)
         llines(ret$dist, ret$gamma, col = Colores[3])
         
         # Fourth model: linear
         ret = variogramLine(fit.media.lin, maxdist = 8)
         llines(ret$dist, ret$gamma, col = Colores[7])
         
         # Fith model: Circular
         ret = variogramLine(fit.media.cir, maxdist = 8)
         llines(ret$dist, ret$gamma, col = Colores[5])
         
       },
       par.settings=acc_par,
       auto.key=list(x=0.6,y=0.35, corner=c(0,1),
                     text=models_name,
                     points=FALSE, lines=TRUE, cex=0.8)
)


########################################### ERRORES ######################################################

# Calidad de Ajuste -> SSErr
models_name <- c("Exponential","Spheric","Gaussian", "Linear","Circular")
media_SSErr <- c(attr(fit.media.exp, "SSErr"),attr(fit.media.sph, "SSErr"),attr(fit.media.gauss, "SSErr"),
                 attr(fit.media.lin, "SSErr"),attr(fit.media.cir, "SSErr"))
SSErr_media <- cbind(models_name,media_SSErr)
colnames(SSErr_media) <- c("Model","SSErr")
SSErr_media



####################################################################################################
#                                       2.2 MEDIANS GLOBAL MODEL                                     #
####################################################################################################
# http://r-spatial.github.io/gstat/reference/fit.StVariogram.html
# Default method fit = 7. This criterion is not supported by theory, but by practice

# median.vgm = data.frame(dist=dist,gamma=medianas)
# median.vgm$np = np #Number of point pairs (Calculado antes)
# median.vgm$dir.hor = rep(0,length(np))
# median.vgm$dir.ver = rep(0,length(np))
# class(median.vgm) = c("gstatVariogram","data.frame")
# 
# # Exponencial
# vg.exp <- vgm(psill=max(median.vgm$gamma)*1.2, model = "Exp", range=max(median.vgm$dist)/2, 
#               nugget = mean(median.vgm$gamma)/4)
# fit.median.exp <- fit.variogram(median.vgm, vg.exp)
# plot(median.vgm,fit.median.exp, main = "Median semivariogram", sub= "Exponential model")
# 
# # Esférico
# vg.sph <- vgm(psill=max(median.vgm$gamma)*0.9, model = "Sph", range=max(median.vgm$dist)/2, 
#               nugget = mean(median.vgm$gamma)/4)
# fit.median.sph <- fit.variogram(median.vgm, vg.sph)
# plot(median.vgm,fit.median.sph, col = "green",main = "Median semivariogram", sub= "Spherical model")
# 
# # Gaussiano
# vg.gauss <- vgm(psill=max(median.vgm$gamma)*1.2, model = "Gau", range=max(median.vgm$dist)/3, 
#                 nugget = mean(median.vgm$gamma)/3)
# fit.median.gauss <- fit.variogram(median.vgm, vg.gauss)
# plot(median.vgm,fit.median.gauss, col = "red",main = "Median semivariogram", sub= "Gaussian model")
# 
# 
# # Lineal
# vg.lin <- vgm(psill=max(median.vgm$gamma)*0.9, model = "Lin", range=max(median.vgm$dist)/2, 
#               nugget = mean(median.vgm$gamma)/4)
# fit.median.lin <- fit.variogram(median.vgm, vg.lin)
# plot(median.vgm,fit.median.lin, col = "purple",main = "Median semivariogram", sub= "Lineal model")
# 
# 
# # Circular
# vg.cir <- vgm(psill=max(median.vgm$gamma)*0.9, model = "Cir", range=max(median.vgm$dist)/2, 
#               nugget = mean(median.vgm$gamma)/4)
# fit.median.cir <- fit.variogram(median.vgm, vg.cir)
# plot(median.vgm,fit.median.cir, col = "dark grey",main = "Median semivariogram", sub= "Circular model")
# 
# 
# ####################################### Plot all fitted models ###########################################
# 
# Colores <- brewer.pal(10,"Paired") # Paleta Paired
# list_col <- c(Colores[1],Colores[2],Colores[3],Colores[7],Colores[5]) #Todos
# models_name <- c("Exponential","Spheric","Gaussian", "Linear", "Circular") #Todos
# 
# acc_par <- list(superpose.line = list(col = list_col, lwd=1.5, lty=1.5))
# 
# xyplot(gamma ~ dist, median.vgm, xlim = c(0,8), pch=20, col= "black", 
#        xlab="Distancia (mm)", ylab="semivarianza",main="Modelo global de medianas",
#        panel = function(...) {
#          panel.xyplot(...)
#          
#          # First model: exponential
#          ret = variogramLine(fit.median.exp, maxdist = 8)
#          llines(ret$dist, ret$gamma, col=Colores[1])
#          
#          # Second model: spheric
#          ret = variogramLine(fit.median.sph, maxdist = 8)
#          llines(ret$dist, ret$gamma, col=Colores[2])
#          
#          # Third model: gaussian
#          ret = variogramLine(fit.median.gauss, maxdist = 8)
#          llines(ret$dist, ret$gamma, col = Colores[3])
#          
#          # Fourth model: linear
#          ret = variogramLine(fit.median.lin, maxdist = 8)
#          llines(ret$dist, ret$gamma, col = Colores[7])
#          
#          # Fifth model: Circular
#          ret = variogramLine(fit.median.cir, maxdist = 8)
#          llines(ret$dist, ret$gamma, col = Colores[5])
#          
#        },
#        par.settings=acc_par,
#        auto.key=list(x=0.5,y=0.35, corner=c(0,1),
#                      text=models_name,
#                      points=FALSE, lines=TRUE, cex=0.8)
# )

###################################### ERRORES #############################################

# # Calidad de Ajuste -> SSErr
# models_name <- c("Exponential","Spheric","Gaussian", "Linear","Circular")
# median_SSErr <- c(attr(fit.median.exp, "SSErr"),attr(fit.median.sph, "SSErr"),attr(fit.median.gauss, "SSErr"),
#                  attr(fit.median.lin, "SSErr"),attr(fit.median.cir, "SSErr"))
# SSErr_median <- cbind(models_name,median_SSErr)
# colnames(SSErr_median) <- c("Model","SSErr")
# SSErr_median

# Save Enviroment

#save.image(file="C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Global_models.RData")


################################################################################################
#                           3. KRIGGING CROSS-VALIDATION  (Validation Data set)                #
################################################################################################

# Read global enviroment
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
setwd(root_path)
load("Global_models.RData")

# Function CV Results Kriging

# Kigriging errors with variogram fitted to model residuals

krige.results <- function(file_names,model,nmax){
  # Function: proportion of normalized above 2.5
  percent_zscore <- function(zscores){
    return((length(zscores[zscores>=2.5]))/length(zscores)*100)
  }
  #Vectores con todos los valores del conjunto 
  predicted <- vector()
  observed <- vector()
  zscores <- vector()
  residuals <- vector()
  name_camp <- vector()
  n_neighbours <- vector()
  #Vectores para guardar medias y desviaciones del resultado de cada campimetria
  #No de todo el conjunto de test
  mean_res <- vector()
  var_res <- vector()
  mean_zscore <- vector()
  var_zscore <- vector()
  mspe <- vector()
  mspe_norm <- vector()
  rmspe <- vector()
  rmspe_norm <- vector()
  cor_obs_pred <- vector()
  cor_pred_res <- vector()
  cor_pred_res_n <- vector()
  percent <- vector()
  name_list <- vector()
  neighbours_list <- vector()
  
  for (i in file_names){
    filepath <- file.path(file.path(paste0(path_sample,"/",paste0(format(i),".csv"))))
    name <- file_path_sans_ext(basename(filepath))
    camp <- read.csv2(filepath, header = TRUE, sep = ";", skipNul = TRUE)
    coordenadas <- camp[,1:2]
    puntos_pred <- SpatialPointsDataFrame(coordenadas, camp, coords.nrs = numeric(0), 
                                          proj4string = CRS(as.character(NA)), bbox = NULL)
    model.pred.krige <- krige.cv(VNorm ~  1,  puntos_pred, model = model, nmax=nmax)
    
    #Results of kriging
    pred.krige <- model.pred.krige$var1.pred
    observed.krige <- model.pred.krige$observed
    zscores.krige <- model.pred.krige$zscore
    residuals.krige <- model.pred.krige$residual
    #Means and variances (For each campimetry)
    mean_res <- c(mean_res,mean(residuals.krige))
    var_res <- c(var_res,var(residuals.krige))
    mean_zscore <- c(mean_zscore,mean(zscores.krige))
    var_zscore <- c(var_zscore,var(zscores.krige))
    #Error measures (For each campimetry)
    mspe <- c(mspe,mean(residuals.krige^2))
    mspe_norm <- c(mspe_norm,mean(zscores.krige^2))
    rmspe <- c(rmspe,sqrt(mean(residuals.krige^2)))
    rmspe_norm <- c(rmspe_norm,sqrt(mean(zscores.krige^2)))
    #correlation observed and predicted, ideally 1
    cor_obs_pred <- c(cor_obs_pred, (cor(observed.krige, pred.krige)))
    #correlation predicted and residual, ideally 0
    cor_pred_res <- c(cor_pred_res, (cor(pred.krige, residuals.krige))) 
    cor_pred_res_n <- c(cor_pred_res_n, (cor(pred.krige, zscores.krige))) 
    #Proportion zscores above 2.5
    percent <- c(percent,percent_zscore(zscores.krige))
    
    #Vector with ALL results
    name_camp <- c(name_camp,rep(name,59))
    name_list <- c(name_list,name)
    x <- rep(coordenadas[,1],length(file_names))
    y <- rep(coordenadas[,2],length(file_names))
    neighbours_list <- c(neighbours_list,nmax)
    n_neighbours <- c(n_neighbours, rep(nmax,59))
    predicted <- c(predicted, pred.krige) 
    observed <- c(observed, observed.krige)
    zscores <- c(zscores, zscores.krige)
    residuals <- c(residuals, residuals.krige)
  }

  # Devolver resultados kriging
  all.result <- data.frame(name_camp,x,y,n_neighbours,observed,predicted,residuals,zscores)
  all.errors <- data.frame(name_list,neighbours_list,mean_res,var_res,mean_zscore,var_zscore,
                           mspe,mspe_norm,rmspe,rmspe_norm,cor_obs_pred,cor_pred_res,cor_pred_res_n,percent)
  return <- list(all.result,all.errors)
  return(return)
}

# Get all results by diferent neighbourhood size in the data validation set
krig_results_by_neighbours <- function(file_names,model,seq_neighbours){
  model_results = data.frame()
  model_errors = data.frame()
    for (i in seq_neighbours){
      results_krige <- krige.results(file_names,model,i)
      values_krige <- results_krige[[1]]
      errors_krige <- results_krige[[2]]
      model_results <- rbind(model_results,values_krige)
      model_errors <- rbind(model_errors,errors_krige)
    }
  return <- list(model_results,model_errors)
  return(return)
}

# Read data validation
data_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Data_IDW"
# Read sample 3
ns <- 3
path_sample <- file.path(paste0(data_path,"/",glue("5260_Patients_Sample{ns}")))
file_names_val <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))
file_names_val_500 <- sample(file_names_val,500)
file_names_val <- file_names_val_500
# 
# media_gauss_global <- krig_results_by_neighbours(file_names_val,fit.media.gauss,c(1:10))
# write.csv2(media_gauss_global[[1]], file.path(getwd(),paste0("media_gauss_results.csv")))
# write.csv2(media_gauss_global[[2]], file.path(getwd(),paste0("media_gauss_errors.csv")))
# media_cir_global <- krig_results_by_neighbours(file_names_val,fit.media.cir,c(1:10))
# write.csv2(media_cir_global[[1]], file.path(getwd(),paste0("media_cir_results.csv")))
# write.csv2(media_cir_global[[2]], file.path(getwd(),paste0("media_cir_errors.csv")))
# media_sph_global <- krig_results_by_neighbours(file_names_val,fit.media.sph,c(1:10))
# write.csv2(media_sph_global[[1]], file.path(getwd(),paste0("media_sph_results.csv")))
# write.csv2(media_sph_global[[2]], file.path(getwd(),paste0("media_sph_errors.csv")))
# media_exp_global <- krig_results_by_neighbours(file_names_val,fit.media.exp,c(1:10))
# write.csv2(media_exp_global[[1]], file.path(getwd(),paste0("media_exp_results.csv")))
# write.csv2(media_exp_global[[2]], file.path(getwd(),paste0("media_exp_errors.csv")))
# media_lin_global <- krig_results_by_neighbours(file_names_val,fit.media.lin,c(1:10))
# write.csv2(media_lin_global[[1]], file.path(getwd(),paste0("media_lin_results.csv")))
# write.csv2(media_lin_global[[2]], file.path(getwd(),paste0("media_lin_errors.csv")))
# median_sph_global <- krig_results_by_neighbours(file_names_val,fit.median.sph,c(1:10))
# write.csv2(median_sph_global[[1]], file.path(getwd(),paste0("median_sph_results.csv")))
# write.csv2(median_sph_global[[2]], file.path(getwd(),paste0("median_sph_errors.csv")))
# median_cir_global <- krig_results_by_neighbours(file_names_val,fit.median.cir,c(1:10))
# write.csv2(median_cir_global[[1]], file.path(getwd(),paste0("median_cir_results.csv")))
# write.csv2(median_cir_global[[2]], file.path(getwd(),paste0("median_cir_errors.csv")))





################################################################################################
#                           4. KRIGGING CROSS-VALIDATION  (Test Data set)                      #
################################################################################################

# Read global enviroment (Fitted Global models)
path_models <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Global_models_500.RData"
load(path_models)

# Select fitted exponencial model and n = 6
n=6
model_exp <- fit.media.exp

# Compute Test -> Sample 4

# Read sample 4
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Data_IDW"
setwd(root_path)
ns <- 4
path_sample <- file.path(paste0(root_path,"/",glue("5260_Patients_Sample{ns}")))
file_names <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))
file_names_test <- file_names
# media_exp_global <- krige.results(file_names_test,model_exp,n)
# write.csv2(media_exp_global[[2]], file.path(getwd(),paste0("media_exp_global_errors_test.csv")))
# write.csv2(media_exp_global[[1]], file.path(getwd(),paste0("media_exp_global_results_test.csv")))


################################################################################################
#                                      5. KRIGGING EXAMPLES                                    #
################################################################################################

# Read global enviroment
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
setwd(root_path)
load("Global_models.RData")

# Eye circle shapefiles
path_eye <- "C:/Users/Agata/MASTER/TFM/Eye_Layers/Eye_cluster"
eye_cir <- readOGR(path_eye, "Eye_Cir") #Eye circular shapefile 
eye_bs <- readOGR(path_eye, "Eye_Blind") #Eye circular + Blind Spot shapefile 
blindSpot <- readOGR(path_eye, "BlindSpot") #Blind Spot shapefile 
grid <- raster(extent(eye_cir))
res(grid) <- 0.5 # Resolution
grid <- as(grid, "SpatialPixelsDataFrame") #Pixel Data Frame
eye_pixels <- grid
#plot(eye_pixels)

# Interpolate map
model <- fit.media.exp
nmax=6
sn = 4
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
data_clust_folder <- glue("Data_{n_cluster}Cluster")
data_sample_folder <- glue("5260_Patients_Sample{sn}")
path_sample <- paste0(root_path,"/",data_clust_folder,"/",data_sample_folder)
file_names_sample <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))
#ejm <- camps_test_c1[1] #"1036210_camp3"
s <- "1036210_camp3"
ejm_idx <- match(s,file_names_sample)

# Plot parameters 
scale.parameter = 1.1
original.bbox <- eye_bs@bbox
xshift = 0  # Shift to right in map units. 
yshift = 0  # Shift to left in map units.
edges = original.bbox
edges[1,] <- (edges[1,] - mean(edges[1,])) * scale.parameter + mean(edges[1,]) + xshift
edges[2,] <- (edges[2,] - mean(edges[2,])) * scale.parameter + mean(edges[2,]) + yshift
bs <- list(blindSpot, col = "blue")

for (i in file_names_sample[ejm_idx]){
  filepath <- file.path(file.path(paste0(path_sample,"/",paste0(format(i),".csv"))))
  name <- file_path_sans_ext(basename(filepath))
  camp <- read.csv2(filepath, header = TRUE, sep = ";", skipNul = TRUE)
  coordenadas <- camp[,1:2]
  puntos_pred <- SpatialPointsDataFrame(coordenadas, camp, coords.nrs = numeric(0), 
                                        proj4string = CRS(as.character(NA)), bbox = NULL)

  
  krige_map <- krige(VNorm ~  1,  puntos_pred, eye_pixels, model = model, nmax=nmax)
  krige_map_cir <- crop(krige_map, eye_cir) # Clip map
  grid.arrange(spplot(puntos_pred, "VNorm", cex = 0.7, do.log = F,
                      main=list(label="Defecto normalizado",cex=1.2),
                      #col.regions=cols,cuts=cuts,
                      sp.layout = list("sp.polygons",eye_bs,col="dark gray"),
                      xlim = edges[1,], ylim = edges[2,],
                      key.space="right"),
               spplot(krige_map_cir,"var1.pred", cex = 0.7, do.log = F,
                      main=list(label="Defecto normalizado",cex=1.2),
                      #col.regions=cols,cuts=cuts,
                      sp.layout = list(blindSpot, fill='darkgray',first = FALSE),
                      xlim = edges[1,], ylim = edges[2,],
                      key.space="right"),ncol=2, widths = c(2,1.66),
               top="Kriging global")
}


#########################################################################################################
#                             VARIOGRAM FUNCTION POINT BY POINT                                          #
#########################################################################################################
# 
# filepath_p <- file.path(file.path(getwd()),paste0(file_names[1],".csv"))
# prueba <- read.csv2(filepath_p, header = TRUE, sep = ";", skipNul = TRUE)
# 
# # Calculate distance matrix
# distancias <- as.matrix(dist(prueba[,1:2], diag = FALSE, 
#                              upper = FALSE, method = "euclidean"))
# 
# 
# # Funcion variograma para cada punto
# variograma_p <- function(p,campimetria){
#   p_value <- campimetria[p,7] #VNorm en el punto (xi)
#   distancias_p <- distancias[-p,1] #Vector distancias punto 1
#   v_values <- campimetria[-p,7] #Valores de los demás puntos
#   v_variogram <- vector() #Vector vacío datos semivariograma
#   j <- 1
#   for(i in v_values){
#     v_variogram[j] <- ((i-p_value)^2)/2
#     j <- j+1
#   }
#   # Graficar
#   points <- data.frame(distancias_p,v_variogram)
#   colnames(points) <- c("x","y")
#   coordinates(points) =  ~ x + y 
#   titulo <- paste("Semivariogram point",p)
#   plot(distancias_p, v_variogram, 
#        xlab="Distance (mm)", ylab="Semivariance", 
#        main=titulo)
#   text(coordinates(points),labels=seq_along(v_values),
#        cex = 0.5, pos = 3) #Labeled points
# }
# 
# variograma_p(1,prueba)






