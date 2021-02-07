#####################################################################################################
#                                                                                                   #
#                                               CAMPIMETERS                                         #
#                                                                                                   #
#                                            Data visualitation                                     #
#                                                                                                   #
#                                       1. Histograms                                               #
#                                         1.1 Global Sample histogram                               #
#                                         1.2 Histograms point by point                             #
#                                       2. Zonal Visualization                                      #
#                                         2.1 Measuring points map                                  #
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
library(glue)         # Format

######################################## 1. HISTOGRAMS ###############################################

# Read all campimeters 

sn = 4 # Sample number
path <- glue("C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Data_IDW/5260_Patients_Sample{sn}") # Random Sample
setwd(path)   # Fixed working directory

## 1.1 Global histogram sample data
all_values <- c()
all_means <- c()
file_names <- list.files(getwd(), pattern = "*.csv")
for(i in file_names){
  filepath <- file.path(file.path(getwd(),paste(i,sep="")))
  aux <- read.csv2(filepath,sep = ";")
  mean <- colMeans(aux)
  mean <- mean[6]
  all_values <- c(all_values,aux$VNorm)
  all_means <- c(all_means,mean)
}


# Means Histogram
# Histogram 
hist(all_means, xlab=glue("Defecto medio"), ylab="Número de campimetrías",
     col="lightblue", main=glue("Muestra {sn}"))

# ¿Percentage of values > Normal reference?
#100*length(which(all_values > 0))/length(all_values) #23.07373

# Histogram 
hist(all_values, xlab=glue("Defecto visual"), ylab="Número de puntos",
     col="lightblue", main=glue("Muestra {sn}"))
# ¿Percentage of values > Normal reference?
#100*length(which(all_values > 0))/length(all_values) #23.07373

## 1.2 Histogram point by point 

# Create dataframe:
all_values_points <- data.frame(dummy$x,dummy$y)
names(all_values_points) <- c("x","y")
file_names <- list.files(getwd(), pattern = "*.csv")
for(i in file_names){
  filepath <- file.path(file.path(getwd(),paste(i,sep="")))
  aux <- read.csv2(filepath,sep = ";")
  all_values_points <- cbind(all_values_points, aux$VNorm)
}
file_names_raw <- file_path_sans_ext(list.files(getwd()))
colnames(all_values_points) <- c("x","y",c(file_names_raw))

all_values_points_values <- all_values_points[,3:ncol(all_values_points)] #Remove x, y columns
all_values_points_values <- as.data.frame(t(all_values_points_values)) #Transposed df
#View(all_values_points_values)

#Visualize all histograms (Point by point)
attach(all_values_points_values)
par(mfrow=c(3,3))
hist(V1);hist(V2);hist(V3);hist(V4);hist(V5);hist(V6);hist(V7);hist(V8);hist(V9)
par(mfrow=c(3,3))
hist(V10);hist(V11);hist(V12);hist(V13);hist(V14);hist(V15);hist(V16);hist(V17);hist(V18)
par(mfrow=c(3,3))
hist(V19);hist(V20);hist(V21);hist(V22);hist(V23);hist(V24);hist(V25);hist(V26);hist(V27)
par(mfrow=c(3,3))
hist(V28);hist(V29);hist(V30);hist(V31);hist(V32);hist(V33);hist(V34);hist(V35);hist(V36)
par(mfrow=c(3,3))
hist(V37);hist(V38);hist(V39);hist(V40);hist(V41);hist(V42);hist(V43);hist(V44);hist(V45)
par(mfrow=c(3,3))
hist(V46);hist(V47);hist(V48);hist(V49);hist(V50);hist(V51);hist(V52);hist(V53);hist(V54)
par(mfrow=c(3,3))
hist(V55);hist(V56);hist(V57);hist(V58);hist(V59)
detach(all_values_points_values)

######################################## 2. ZONAL VISUALIZATION ######################################

## 2.1 Measuring points map
myColours <- brewer.pal(6,"Blues") # Blue palette

# Dummy campimeter to extract points
dummy <- read.csv2(file.path(file.path(getwd()),paste0(file_names[1])), 
                   header = TRUE, sep = ";", skipNul = TRUE)
coordinates <- dummy[,1:2]
points_pred <- SpatialPointsDataFrame(coordinates, dummy, coords.nrs = numeric(0), 
                                      proj4string = CRS(as.character(NA)), bbox = NULL)
x <- dummy$x
y <- dummy$y
labels <- as.character(row.names(dummy))
map <- xyplot(y~x)
xyplot(y~x,pch=16, col = myColours[3], main=list(label="Puntos de medida", cex=1))
trellis.focus("panel",1,1)
panel.text(x=map$panel.args[[1]]$x,y=map$panel.args[[1]]$y,labels = labels, pos=3, cex=0.6, col=myColours[6])
trellis.unfocus()

# Visualization by measures

#Previous dataset histogram by points (no transposed and x, y columns)
View(all_values_points)
data <- all_values_points
data$VNorm_mean <- rowMeans(data[c(-1,-2)])
data$VNorm_median <- rowMedians(as.matrix(data[c(-1,-2)]))
data$VNorm_skewness <- rowSkewness(as.matrix(data[c(-1,-2)]))
data$VNorm_kurtosis <- rowKurtosis(as.matrix(data[c(-1,-2)]))
data_measures <- data[c("x","y","VNorm_mean","VNorm_median", 
                          "VNorm_skewness", "VNorm_kurtosis")]

# 3D visualization selected measures
plot3d(data_measures$x, data_measures$y, data_measures$VNorm_mean, xlab="X", ylab="Y", zlab="Z")
plot3d(data_measures$x, data_measures$y, data_measures$VNorm_median, xlab="X", ylab="Y", zlab="Z")
plot3d(data_measures$x, data_measures$y, data_measures$VNorm_skewness, xlab="X", ylab="Y", zlab="Z")
plot3d(data_measures$x, data_measures$y, data_measures$VNorm_kurtosis, xlab="X", ylab="Y", zlab="Z")
#Se aprecia una tendencia Sureste-Noroeste (Mejor a peor visión)

# Bubble visualization
cols <- brewer.pal(5, "Spectral")
coordenadas <- data[,1:2]
data_all_sp <- SpatialPointsDataFrame(coordenadas, data, coords.nrs = numeric(0), 
                                        proj4string = CRS(as.character(NA)), bbox = NULL)
spplot(data_all_sp, "VNorm_mean", do.log=F, main=list(label="Medias de las desviaciones",cex=1.5), 
       col.regions=cols)
spplot(data_all_sp, "VNorm_median", do.log=F, main=list(label="Medianas de las desviaciones",cex=1.5),
       col.regions=cols)
spplot(data_all_sp, "VNorm_skewness", do.log=F, main=list(label="Asimetría de los valores",cex=1.5),
       col.regions=cols)
spplot(data_all_sp, "VNorm_kurtosis", do.log=F, main=list(label="Curtosis de los valores",cex=1.5),
       col.regions=cols)





