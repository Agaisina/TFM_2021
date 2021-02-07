#####################################################################################################
#                                                                                                   #
#                                               CAMPIMETERS                                         #
#                                                                                                   #
#                                           GLOBAL MODELS - IDW                                     #
#                                                                                                   #
#                         Previous: Prepare data                                                    #                                                                          #
#                         0. Read data                                                              #
#                         1. IDW Model                                                              #
#                           1.1 Function for calculation of IDW CV results and errors               #
#                           1.2 Function to get all results and errors by diferent                  #
#                               neighbourhood size and diferent power values                        #
#                           1.3 Calculations                                                        #
#                         2. Display training errors and results                                    #                          
#                           2.1 Mean of all training campimeters errors                             #
#                         3. Test errors and results                                                #
#                                                                                                   #
#####################################################################################################


library(lattice)      # Graphics
library(sp)           # Spatial data
library(fs)           # File manipulations
library(tidyverse)    # Data manipulation
library(tools)        # Files manipulation
library(stringr)      # Character manipulation
library(gstat)        # Use gstat's idw routine
library(sp)           # Used for the spsample function
library(tools)        # Files manipulation (file_path_sans_ext function)
library(spm)          # idw with k-fold cross validation
library(glue)         # Equivalent of python 'format' (str)
library(RColorBrewer) # Color palettes
library(RcmdrMisc)    # Function Hist()

#########################################################################################################

# ¡Only run once!
# Create a data frame containing all campimeters:
# Each row is a campimeter and columns are 59 Defect Values (Measuring points)

# Read Sample Folder
# sn = 1 # Sample num
# root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Data_IDW"
# data_sample_folder <- glue("5260_Patients_Sample{sn}")
# path <- paste0(root_path,"/",data_sample_folder)
# setwd(path)
# file_names <- file_path_sans_ext(list.files(file.path(getwd()), pattern = "*.csv"))
# 
# df_campimeters <- data.frame()
# file_names <- list.files(getwd(), pattern = "*.csv")
# for(i in file_names){
#         filepath <- file.path(file.path(getwd(),paste(i,sep="")))
#         aux <- read.csv2(filepath,sep = ";")
#         df_campimeters <- rbind(df_campimeters, aux$VNorm) #Df of Values
# }
# 
# file_names_raw <- file_path_sans_ext(list.files(getwd()))
# row.names(df_campimeters) <- c(file_names_raw)
# write.csv2(df_campimeters, 
#            file.path(paste0(root_path,"/","5260_Values_df","/",glue("df_5260_Values_campimeters{sn}.csv"))))



#########################################################################################################
#                                            0. READ DATA                                               #
#########################################################################################################
# Read data
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Data_IDW"
setwd(root_path)
# Read sample 2
ns <- 2
path_sample <- file.path(paste0(root_path,"/",glue("5260_Patients_Sample{ns}")))
file_names <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))

# # Split data into train, validation and test (70% 15% 15%)
# set.seed(123)
# file_names_train <- sample(file_names,(0.50*length(file_names)))
# file_names_valtest  <- file_names[!file_names %in% file_names_train]
# file_names_val <- sample(file_names_valtest,(0.5*length(file_names_valtest)))
# file_names_test <- file_names_valtest[!file_names_valtest %in% file_names_val]

# Split data into train -> 500 campimeters
set.seed(123)
file_names_train <- sample(file_names,500)


# Dummy campimeter
# Coordinates points -> Read a dummy campimeter and extract it
path <- path_sample
setwd(path)   # Fixed working directory
file_names <- list.files(getwd(), pattern = "*.csv")
dummy_camp <- read.csv2(file.path(file.path(getwd(),paste(file_names[1],sep=""))))
x <- dummy_camp$x
y <- dummy_camp$y
coordenadas <- dummy_camp[,1:2]
dummy_sp <- SpatialPointsDataFrame(coordenadas, dummy_camp, coords.nrs = numeric(0), 
                                      proj4string = CRS(as.character(NA)), bbox = NULL)


#########################################################################################################
#                                            1. IDW                                                     #
#########################################################################################################
path_sample <- path_sample

# 1.1 Function for calculate and save IDW CV errors

idw.results <- function(file_names,nmax,power){
  
  #Vectores con todos los valores del conjunto Test
  predicted <- vector()
  observed <- vector()
  #zscores <- vector() #IDW is not a geostatistical method
  residuals <- vector()
  name_camp <- vector()
  n_neighbours <- vector()
  n_power <- vector()
  #Vectores para guardar medias y desviaciones del resultado de cada campimetria
  #No de todo el conjunto de test
  mean_res <- vector()
  var_res <- vector()
  #mean_zscore <- vector()
  #var_zscore <- vector()
  mspe <- vector()
  #mspe_norm <- vector()
  rmspe <- vector()
  #rmspe_norm <- vector()
  cor_obs_pred <- vector()
  cor_pred_res <- vector()
  #cor_pred_res_n <- vector())
  #percent <- vector()
  name_list <- vector()
  neighbours_list <- vector()
  power_list <- vector()
  
  for (i in file_names){
    filepath <- file.path(file.path(paste0(path_sample,"/",paste0(format(i),".csv"))))
    name <- file_path_sans_ext(basename(filepath))
    camp <- read.csv2(filepath, header = TRUE, sep = ";", skipNul = TRUE)
    coordenadas <- camp[,1:2]
    puntos_pred <- SpatialPointsDataFrame(coordenadas, camp, coords.nrs = numeric(0), 
                                          proj4string = CRS(as.character(NA)), bbox = NULL)
    idw_camp <- gstat(formula = VNorm ~ 1, #intercept only model
                      locations=puntos_pred, 
                      nmax=nmax,
                      set=list(idp = power))
    model.pred.idw <- gstat.cv(idw_camp) #Default LOOCV
    
    #Results of IDW
    pred.idw <- model.pred.idw$var1.pred
    observed.idw <- model.pred.idw$observed
    #zscores.krige <- model.pred.krige$zscore
    residuals.idw <- model.pred.idw$residual
    # Means and variances (For each campimetry)
    mean_res <- c(mean_res,mean(residuals.idw))
    var_res <- c(var_res,var(residuals.idw))
    #mean_zscore <- c(mean_zscore,mean(zscores.krige))
    #var_zscore <- c(var_zscore,var(zscores.krige))
    # Error measures (For each campimetry)
    mspe <- c(mspe,mean(residuals.idw^2))
    #mspe_norm <- c(mspe_norm,mean(zscores.krige^2))
    rmspe <- c(rmspe,sqrt(mean(residuals.idw^2)))
    #rmspe_norm <- c(rmspe_norm,sqrt(mean(zscores.krige^2)))
    # Correlation observed and predicted, ideally 1
    cor_obs_pred <- c(cor_obs_pred, (cor(observed.idw,pred.idw)))
    # Correlation predicted and residual, ideally 0
    cor_pred_res <- c(cor_pred_res, (cor(pred.idw, residuals.idw))) 
    #cor_pred_res_n <- c(cor_pred_res_n, (cor(pred.krige, zscores.krige))) 
    # Proportion zscores above 2.5
    #percent <- c(percent,percent_zscore(zscores.krige))
    
    # Vector with ALL results
    name_camp <- c(name_camp,rep(name,59))
    name_list <- c(name_list,name)
    x <- rep(coordenadas[,1],length(file_names))
    y <- rep(coordenadas[,2],length(file_names))
    neighbours_list <- c(neighbours_list,nmax)
    n_neighbours <- c(n_neighbours, rep(nmax,59))
    power_list <- c(power_list,power)
    n_power <- c(n_power,rep(power,59))
    predicted <- c(predicted, pred.idw)
    observed <- c(observed, observed.idw)
    #zscores <- c(zscores, zscores.krige)
    residuals <- c(residuals, residuals.idw)
  }
  
  # Devolver resultados IDW
  #all.result <- data.frame(name_camp,n_neighbours,n_power,observed,predicted,residuals)
  all.result <- data.frame(name_camp,x,y,n_neighbours,n_power,
                           observed,predicted,residuals)
  all.errors <- data.frame(name_list,neighbours_list,power_list,mean_res,var_res,
                           mspe,rmspe,cor_obs_pred,cor_pred_res)
  return <- list(all.result,all.errors)
  return(return)
}

# 1.2 Function to get all results and errors by diferent neighbourhood size and 
# diferent power values

idw_results_by_neighbours_power <- function(file_names,list_neighbors,list_power){
  model_results = data.frame()
  model_errors = data.frame()
  for (i in list_neighbors){
    for (j in list_power){
      results_idw <- idw.results(file_names,i,j)
      values_idw <- results_idw[[1]]
      errors_idw <- results_idw[[2]]
      model_results <- rbind(model_results,values_idw)
      model_errors <- rbind(model_errors,errors_idw)
    }
  }
  return <- list(model_results,model_errors)
  return(return)
}

# 1.3 Calculations

# list_n <- seq(from = 1, to = 15, by = 1)
# list_b <- seq(from = 1, to = 5, by = 0.25 )
# 
# idw_train <- idw_results_by_neighbours_power(file_names_train,list_n,list_b)
# 
# write.csv2(idw_train[[1]], file.path(root_path,paste0("idw_train_results.csv")))
# write.csv2(idw_train[[2]], file.path(root_path,paste0("idw_train_errors.csv")))


#########################################################################################################
#                                            2. Visualize training errors and results                   #
#########################################################################################################

root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/RESULTS/IDW"
setwd(root_path)
errors_train <- read.csv2(file.path(root_path,"idw_train_errors.csv"), sep = ";")
results_train <-read.csv2(file.path(root_path,"idw_train_results.csv"), sep = ";")

# 2.1 Mean of all training campimeters errors 

errors_train_mean <- errors_train %>% group_by(neighbours_list,power_list) %>%
                                      summarize_all(mean)
errors_train_mean_n <- errors_train %>% group_by(neighbours_list) %>%
                                      summarize_all(mean)
View(errors_train_mean_n)
# Plot mean errors by number of neighbours

blues <- brewer.pal(5,"Blues")

par(mfrow=c(1,3))
plot(x=errors_train_mean_n$neighbours_list, y=errors_train_mean_n$mean_res, 
     col = blues[4], type = "b", pch = 19,
     xlab = "Número vecinos", ylab = "Error medio")
plot(x=errors_train_mean_n$neighbours_list, y=errors_train_mean_n$mspe, 
     col = blues[4], type = "b", pch = 19,
     xlab = "Número vecinos", ylab = "MSE")
plot(x=errors_train_mean_n$neighbours_list, y=errors_train_mean_n$cor_obs_pred, 
     col = blues[4], type = "b", pch = 19,
     xlab = "Número vecinos", ylab = "Correlación observado-predicho")
mtext("Evolución de las medidas de error", side = 3, line = -2.5, outer = TRUE)
dev.off()

# n = 5 y p = 2


#########################################################################################################
#                                            3. Test errors and results                                 #
#########################################################################################################


# 3.1 Compute Test 
# 
# root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Data_IDW"
# setwd(root_path)
# # Read sample 4
# ns <- 4
# path_sample <- file.path(paste0(root_path,"/",glue("5260_Patients_Sample{ns}")))
# file_names <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))
# file_names_test <- file_names # Sample 4 
# 
# idw_test <- idw_results_by_neighbours_power(file_names_test,5,2)
# write.csv2(idw_test[[1]], file.path(root_path,paste0("idw_test_results.csv")))
# write.csv2(idw_test[[2]], file.path(root_path,paste0("idw_test_errors.csv")))

# 3.2 Display test errors

root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/RESULTS/IDW"
setwd(root_path)
errors_test <- read.csv2(file.path(root_path,"idw_test_errors.csv"), sep = ";")
results_test <-read.csv2(file.path(root_path,"idw_test_results.csv"), sep = ";")

# Mean of all test campimeters errors 

errors_test_mean <- errors_test %>% group_by(neighbours_list,power_list) %>%
                                    summarize_all(mean)
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
title("Modelo IDW Global", line = 0.5, outer = TRUE, cex.main=1.3)
dev.off()

################################################################################################
#                                      4. IDW EXAMPLES                                         #
################################################################################################


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
p=2
nmax=5
sn = 4
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
data_clust_folder <- glue("Data_{n_cluster}Cluster")
data_sample_folder <- glue("5260_Patients_Sample{sn}")
path_sample <- paste0(root_path,"/",data_clust_folder,"/",data_sample_folder)
file_names_sample <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))
ejm <- camps_test_c1[1] #"1036210_camp3"
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
  
  
  idw_map <- idw(formula = VNorm ~ 1, #intercept only model
                   locations=puntos_pred, eye_pixels,
                   nmax=nmax,
                   idp = p)
  idw_map_cir <- crop(idw_map, eye_cir) # Clip map
  grid.arrange(spplot(puntos_pred, "VNorm", cex = 0.7, do.log = F,
                      main=list(label="Defecto normalizado",cex=1.2),
                      #col.regions=cols,cuts=cuts,
                      sp.layout = list("sp.polygons",eye_bs,col="dark gray"),
                      xlim = edges[1,], ylim = edges[2,],
                      key.space="right"),
               spplot(idw_map_cir,"var1.pred", cex = 0.7, do.log = F,
                      main=list(label="Defecto normalizado",cex=1.2),
                      #col.regions=cols,cuts=cuts,
                      sp.layout = list(blindSpot, fill='darkgray',first = FALSE),
                      xlim = edges[1,], ylim = edges[2,],
                      key.space="right"),ncol=2, widths = c(2,1.66),
               top=glue("IDW global"))
}


