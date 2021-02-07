#####################################################################################################
#                                                                                                   #
#                                               CAMPIMETERS                                         #
#                                                                                                   #
#                                        (8) (10) CLUSTER VISUAL FIELD                              #
#                                                                                                   #
#                           1. Measuring points and (8) (10) cluster visual fields visualization    #
#                           2. Creation of new variables - 'Global Indices'  (Only run once)        #
#                           3. Clustering of patients                                               #
#                              3.1 Determining optimal number of clusters                           #
#                              3.2 Compute clustering                                               #
#                              3.3 Visualization of clusters                                        #
#                                  3.3.2 Visualization of all clusters (Unique plot)                #
#                                  3.3.2 Visualization by clusters (Multiple plots)                 #
#                           4. Assign training dataset to clusters                                  #
#                              4.1 Display some prediction examples                                 #
#                           5. Fitting 'Average' variograms models   (Training data set)            # 
#                              5.1 Compute means variograms                                         #
#                              5.2 Fit variograms models                                            #
#                           6. Kriging CV (Validation data set)                                     #
#                              6.1 Assign test data set into cluster                                #
#                              6.2 Compute kriging CV results                                       #
#                           7. Compare variogram models (Kriging CV - Validation dataset)           #
#                           8. Kriging CV (Test data set)                                           #
#                              8.1 Assign test data set into cluster                                #
#                              8.2 Compute kriging CV results                                       #
#                              8.3 Display error measurements                                       #
#                           9. Kriging examples                                                     #
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
library(stringi)      # Function stri_sub
library(gridExtra)    # Tables as images

#########################################################################################################
#     0. Load data and define number of visual field clusters                                           #
#########################################################################################################

n_cluster <- 10 # Number of visual field areas:
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
#     2. Calculation of 'Global Indices' as new variables -> Only run once                              #
#########################################################################################################

# Function for calculating global indices of a campimeter
global_indices <- function(campimeter,n_clust){
        camp_c <- campimeter %>%
                group_by(Cluster) %>%
                summarize(N = n(),
                          MD_C = sum(VNorm)/N, #Mean defect by cluster
                          sLV_C = sqrt(sum((VNorm - MD_C)^2)/(N - 1))) #Square root of loss of variance by cluster
        camp_t <- campimeter %>%
                summarize(N = n(),
                          MD_T = sum(VNorm)/N, #Mean defect (global)
                          sLV_T = sqrt(sum((VNorm - MD_T)^2)/(N - 1))) #Square root of loss of variance (global)
        camp_c <- na.omit(camp_c)
        camp_MD_T <- data.frame(matrix(unlist(camp_t$MD_T), 1, byrow=T))
        colnames(camp_MD_T) <- paste("MD_T")
        camp_sLV_T <- data.frame(matrix(unlist(camp_t$sLV_T), 1, byrow=T))
        colnames(camp_sLV_T) <- paste("sLV_T")
        camp_MD_C <- data.frame(matrix(unlist(camp_c$MD_C), 1, byrow=T))
        colnames(camp_MD_C) <- paste("MD_C",seq(1:n_clust),sep="")
        camp_sLV_C <- data.frame(matrix(unlist(camp_c$sLV_C), 1, byrow=T))
        colnames(camp_sLV_C) <- paste("sLV_C",seq(1:n_clust),sep="")
        camp <- cbind(camp_MD_T,camp_sLV_T,camp_MD_C,camp_sLV_C)
        return(camp)
}

# # ¡Only run once!
# # Create a data frame containing all campimeters:
# # 1. Each row is a campimeter and columns are Global Indices (22 variables)
# # 2. Each row is a campimeter and columns are 59 Defect Values (Measuring points)
# n_cluster = 10
# sn = 4 # Sample num
# root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
# data_cluster_folder <- glue("Data_{n_cluster}Cluster")
# data_sample_folder <- glue("5260_Patients_Sample{sn}")
# path <- paste0(root_path,"/",data_cluster_folder,"/",data_sample_folder)
# setwd(path)
# file_names <- file_path_sans_ext(list.files(file.path(getwd()), pattern = "*.csv"))
# 
# df_campimeters_v <- data.frame()
# df_campimeters_GI <- data.frame()
# file_names <- list.files(getwd(), pattern = "*.csv")
# for(i in file_names){
#         filepath <- file.path(file.path(getwd(),paste(i,sep="")))
#         aux <- read.csv2(filepath,sep = ";")
#         camp <- global_indices(aux,n_cluster)   #Calculation Global Indices
#         df_campimeters_v <- rbind(df_campimeters_v, aux$VNorm)   #Df of Values
#         df_campimeters_GI <- rbind(df_campimeters_GI, camp)      #Df of Global Indices
# }
# 
# file_names_raw <- file_path_sans_ext(list.files(getwd()))
# row.names(df_campimeters_v) <- c(file_names_raw)
# row.names(df_campimeters_GI) <- c(file_names_raw)
# 
# write.csv2(df_campimeters_v,
#            file.path(paste0(root_path,"/",data_cluster_folder,"/",
#                             "5260_Values_df","/",glue("df_5260_Values_campimeters{sn}.csv"))))
# write.csv2(df_campimeters_GI,
#            file.path(paste0(root_path,"/",data_cluster_folder,"/",
#                             "5260_Global_Indices_df","/",glue("df_5260_GI_campimeters{sn}.csv"))))

#########################################################################################################
#     3. Clustering of patients                                                                         #
#########################################################################################################

# 3.1 Determining optimal number of clusters

# # Load data and prepare:
# n_cluster <- n_cluster # Number of visual field areas:
# eye_layout <- eye_layout #"Eye_clusters10" #"Eye_clusters8"
# sn = 1 # Sample num (1 = 5260 campimeters, 2 = 10520 campimeters...)
# data_clust_folder <- glue("Data_{n_cluster}Cluster")
# root_path = "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/"
# df <- read.csv2(paste0(root_path,data_clust_folder,"/5260_Global_Indices_df",glue("/df_5260_GI_campimeters{sn}.csv")))
# df <- subset(df, select = -c(MD_T, sLV_T))
# rownames(df) <- df$X
# df$X <- NULL
# # Data preparation: (No) Standardize the data
# #data <- scale(df)
# data <- df
# 
# # Elbow method
# fviz_nbclust(data, kmeans, method="wss", k.max=15) +
#         geom_vline(xintercept = 4, linetype = 2)+
#         labs(subtitle = "Elbow method")
# 
# # Silhouette method
# fviz_nbclust(data, kmeans, method = "silhouette",k.max=15)+
#        labs(subtitle = "Silhouette method")
#
# # Gap statistic
# # nboot = 50 to keep the function speedy.
# # recommended value: nboot= 500 for your analysis.
# # Use verbose = FALSE to hide computing progression.
# set.seed(123)
# fviz_nbclust(data, FUN=kmeans, nstart = 25,  method = "gap_stat", k.max = 15, nboot = 100)+
#        labs(subtitle = "Gap statistic method") #14
#
# # clvalid - Internal validation and stability validation
# # https://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
# intern <- clValid(data, 3:15, clMethods=c("kmeans"), 
#                   method = c("average","ward", validation="internal"))
# stability <- clValid(data, 3:15, clMethods=c("kmeans"),
#                      method = c("average","ward", validation="stabilty"))
# 
# # NbClust(): provides 30 indices for choosing the best number of clusters
# nbclust_out <- NbClust(
#         data = data,
#         distance = "euclidean",
#         min.nc = 10, # minimum number of clusters
#         max.nc = 15, # maximum number of clusters
#         method = "kmeans" # one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans"
# )
# # create a dataframe of the optimal number of clusters
# nbclust_plot <- data.frame(clusters = nbclust_out$Best.nc[1, ])
# write.csv2(nbclust_plot, file.path(paste0(root_path,data_clust_folder),
#                                    paste0("Nbclust_8Cluster.csv")))
# # create plot
# nbclust_plot <- subset(nbclust_plot, clusters >= 8 & clusters <= 15)
# ggplot(nbclust_plot) + geom_bar(aes(x = clusters), fill=myColours[5]) +
#         labs(x = "Número de clusters", y = "Frecuencia entre índices", title = "Número óptimo de clusters") +
#         theme(axis.text.x = element_text(size = 10))

# 3.2 Compute clustering

# Load and prepare data:
n_cluster <- n_cluster 
eye_layout <- eye_layout 
sn = 1 # Sample num (1 = 5260 campimeters, 2 = 10520 campimeters...)
data_clust_folder <- glue("Data_{n_cluster}Cluster")
root_path = "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
df <- read.csv2(paste0(root_path,"/",data_clust_folder,"/5260_Global_Indices_df",
                       glue("/df_5260_GI_campimeters{sn}.csv")))
# Variables: MD and sLV:
df <- subset(df, select = -c(MD_T, sLV_T))
# Variables: only MD
#MD_variables <- paste0("MD_C",1:n_cluster)
#df <- subset(df, select = -c(MD_T, sLV_T))
#df <- subset(df, select = MD_variables)
rownames(df) <- df$X
df$X <- NULL
# Data preparation: (NO) Standardize the data
#data <- scale(df)
data <- df

# Compute kmeans clustering 
set.seed(123)
n_clust_p <- 12   #Number of patient clusters 
km <- kmeans(data, n_clust_p, nstart = 50, iter.max=50)
#pam.10 <- pam(data, 10, metric = "euclidean", stand = FALSE)
clust <- km
camps_c1 <- rownames(as.data.frame(clust$cluster[clust$cluster==1]))
camps_c2 <- rownames(as.data.frame(clust$cluster[clust$cluster==2]))
camps_c3 <- rownames(as.data.frame(clust$cluster[clust$cluster==3]))
camps_c4 <- rownames(as.data.frame(clust$cluster[clust$cluster==4]))
camps_c5 <- rownames(as.data.frame(clust$cluster[clust$cluster==5]))
camps_c6 <- rownames(as.data.frame(clust$cluster[clust$cluster==6]))
camps_c7 <- rownames(as.data.frame(clust$cluster[clust$cluster==7]))
camps_c8 <- rownames(as.data.frame(clust$cluster[clust$cluster==8]))
camps_c9 <- rownames(as.data.frame(clust$cluster[clust$cluster==9]))
camps_c10 <- rownames(as.data.frame(clust$cluster[clust$cluster==10]))
camps_c11 <- rownames(as.data.frame(clust$cluster[clust$cluster==11]))
camps_c12 <- rownames(as.data.frame(clust$cluster[clust$cluster==12]))
#camps_c13 <- rownames(as.data.frame(clust$cluster[clust$cluster==13]))
#camps_c14 <- rownames(as.data.frame(clust$cluster[clust$cluster==14]))
#camps_c15 <- rownames(as.data.frame(clust$cluster[clust$cluster==15]))

# 3.3 Visualization of clusters

# Subsetting data of different clusters of original data
data_values <- read.csv2(paste0(root_path,"/",data_clust_folder,"/5260_Values_df",
                                glue("/df_5260_Values_campimeters{sn}.csv")))
rownames(data_values) <- data_values$X
data_values$X <- NULL
data_c1 <- data_values[camps_c1,]
data_c2 <- data_values[camps_c2,]
data_c3 <- data_values[camps_c3,]
data_c4 <- data_values[camps_c4,]
data_c5 <- data_values[camps_c5,]
data_c6 <- data_values[camps_c6,]
data_c7 <- data_values[camps_c7,]
data_c8 <- data_values[camps_c8,]
data_c9 <- data_values[camps_c9,]
data_c10 <- data_values[camps_c10,]
data_c11 <- data_values[camps_c11,]
data_c12 <- data_values[camps_c12,]
#data_c13 <- data_values[camps_c13,]
#data_c14 <- data_values[camps_c14,]
#data_c15 <- data_values[camps_c15,]

# 3.3.1 Visualization of all clusters (Unique plot)

# Create list of data clusters
cluster_list <- c()      # List of data clusters
for (n in (1:n_clust_p)){
        data <- paste0('data_c',n)
        cluster_list <- append(cluster_list,data)
}

# Function: Transpose clustering data, add coordinates and calculate visual defect statistics 
# (All clusters)
stats_points_all <- function(data_list){
        i = 1
        data_all <- data.frame(coordinates)
        coords <- coordinates  # From dummy campimeter 
        for (data in data_list){
                aux <- get(data)
                data_t <- as.data.frame(t(aux)); data_t <- cbind(y = y, data_t); data_t <- cbind(x = x, data_t)
                data_t$VNorm_mean <- rowMeans(data_t)
                data_t$VNorm_median <- rowMedians(as.matrix(data_t))
                data_t$VNorm_skw <- rowSkewness(as.matrix(data_t))
                data_t$VNorm_kurt <- rowKurtosis(as.matrix(data_t))
                data_t <- select(data_t, VNorm_mean, VNorm_median, VNorm_skw, VNorm_kurt)
                col_names <- c(glue("VNorm_mean{i}"),glue("VNorm_medan{i}"),
                               glue("VNorm_skw{i}"),glue("VNorm_kurt{i}"))
                colnames(data_t) <- col_names
                data_all <- cbind(data_all,data_t)
                i = i + 1
        }
        data_all_sp <- SpatialPointsDataFrame(coords, data_all, coords.nrs = numeric(0), 
                                              proj4string = CRS(as.character(NA)), bbox = NULL)
        return(data_all_sp)
}

# Calculate statistics for all patient clusters
data_clust_stats <- stats_points_all(cluster_list)

# Function: Plot (average) defect for all clusters
eye_plot_all <- function(data_all_stats,statistic,title,cols,cuts){
        # Select desire statistic
        vars_names <- paste0(glue("VNorm_{statistic}"), 1:n_clust_p) 
        data_all_stats@data <- data_all_stats@data %>% 
                select(all_of(vars_names))
        clust_names = paste0("C",(1:n_clust_p))
        colnames(data_all_stats@data) = clust_names
        # Plot parameters 
        scale.parameter = 1.1
        original.bbox <- eye_clusters@bbox
        xshift = 0  # Shift to right in map units. 
        yshift = 0  # Shift to left in map units.
        edges = original.bbox
        edges[1,] <- (edges[1,] - mean(edges[1,])) * scale.parameter + mean(edges[1,]) + xshift
        edges[2,] <- (edges[2,] - mean(edges[2,])) * scale.parameter + mean(edges[2,]) + yshift
        spplot(data_all_stats, clust_names, cex = 0.6, do.log = F,
               main=list(label=title,cex=1.2),
               col.regions=cols,cuts=cuts,
               #color.key=TRUE,
               sp.layout = list("sp.polygons",eye_clusters,col="dark gray"),
               xlim = edges[1,], ylim = edges[2,],
               key.space="right",
               as.table = TRUE)
}

cuts <- c(-50,10,50,100,150,200,300)
#cuts <- c(-50,10,60,110,160,250,300,350,400,500)
cols <- brewer.pal(length(cuts), "Spectral")
cols <- rev(cols)

eye_plot_all(data_clust_stats,'mean','Defecto medio por cluster',cols,cuts)


# Distance between centroids

dist_cent <- dist(clust$centers, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
print(dist_cent)

# 3.3.2 Visualization by clusters (Multiple plots)

# # Define plot propierties
# cuts <- c(-50,10,50,100,150,200,300)
# cols <- brewer.pal(length(cuts), "Spectral")
# cols <- rev(cols)
# 
# # Function for plotting points and 10 cluster visual field
# eye_plot <- function(data_sp,field,title,cols,cuts){
#         scale.parameter = 1.1
#         original.bbox <- eye_clusters@bbox
#         xshift = 0  # Shift to right in map units. 
#         yshift = 0  # Shift to left in map units.
#         edges = original.bbox
#         edges[1,] <- (edges[1,] - mean(edges[1,])) * scale.parameter + mean(edges[1,]) + xshift
#         edges[2,] <- (edges[2,] - mean(edges[2,])) * scale.parameter + mean(edges[2,]) + yshift
#         spplot(data_sp, field, do.log = F,
#                main=list(label=title,cex=1.2),
#                col.regions=cols,cuts=cuts,
#                sp.layout = list("sp.polygons",eye_clusters,col="gray"),
#                xlim = edges[1,], ylim = edges[2,],
#                key.space="right")}
# 
# # Function: Transpose clustering data, add coordinates and calculate visual defect statistics 
# stats_points <- function(data){
#         data_t <- as.data.frame(t(data)); data_t <- cbind(y = y, data_t); data_t <- cbind(x = x, data_t)
#         data_t$VNorm_mean <- rowMeans(data_t)
#         data_t$VNorm_median <- rowMedians(as.matrix(data_t))
#         data_t$VNorm_skw <- rowSkewness(as.matrix(data_t))
#         data_t$VNorm_kurt <- rowKurtosis(as.matrix(data_t))
#         coords <- coordinates
#         # Keep only variables
#         data_t <- select(data_t, VNorm_mean, VNorm_median, VNorm_skw, VNorm_kurt)
#         data_sp <- SpatialPointsDataFrame(coords, data_t, coords.nrs = numeric(0), 
#                                              proj4string = CRS(as.character(NA)), bbox = NULL)
#         return(data_sp)
# }
# 
# # Cluster 1
# data_c1_sp <- stats_points(data_c1)
# eye_plot(data_c1_sp,"VNorm_mean","C1 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c1_sp,"VNorm_median","C1 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 2
# data_c2_sp <- stats_points(data_c2)
# eye_plot(data_c2_sp,"VNorm_mean","C2 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c2_sp,"VNorm_median","C2 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 3
# data_c3_sp <- stats_points(data_c3)
# eye_plot(data_c3_sp,"VNorm_mean","C3 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c3_sp,"VNorm_median","C3 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 4
# data_c4_sp <- stats_points(data_c4)
# eye_plot(data_c4_sp,"VNorm_mean","C4 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c4_sp,"VNorm_median","C4 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 5
# data_c5_sp <- stats_points(data_c5)
# eye_plot(data_c5_sp,"VNorm_mean","C5 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c5_sp,"VNorm_median","C5 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 6
# data_c6_sp <- stats_points(data_c6)
# eye_plot(data_c6_sp,"VNorm_mean","C6 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c6_sp,"VNorm_median","C6 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 7
# data_c7_sp <- stats_points(data_c7)
# eye_plot(data_c7_sp,"VNorm_mean","C7 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c7_sp,"VNorm_median","C7 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 8
# data_c8_sp <- stats_points(data_c8)
# eye_plot(data_c8_sp,"VNorm_mean","C8 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c8_sp,"VNorm_median","C8 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 9
# data_c9_sp <- stats_points(data_c9)
# eye_plot(data_c9_sp,"VNorm_mean","C9 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c9_sp,"VNorm_median","C9 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 10
# data_c10_sp <- stats_points(data_c10)
# eye_plot(data_c10_sp,"VNorm_mean","C10 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c10_sp,"VNorm_median","C10 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 11
# data_c11_sp <- stats_points(data_c11)
# eye_plot(data_c11_sp,"VNorm_mean","C11 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c11_sp,"VNorm_median","C11 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 12
# #data_c12_sp <- stats_points(data_c12)
# #eye_plot(data_c12_sp,"VNorm_mean","C12 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c12_sp,"VNorm_median","C12 - Defecto visual (mediana)",cols,cuts)
# 
# # Cluster 13
# #data_c13_sp <- stats_points(data_c13)
# #eye_plot(data_c10_sp,"VNorm_mean","C13 - Defecto visual (media)",cols,cuts)
# #eye_plot(data_c10_sp,"VNorm_median","C13 - Defecto visual (mediana)",cols,cuts)


#########################################################################################################
#     4. Assign training dataset to clusters                                                            #                                                                       #
#########################################################################################################
set.seed(123)

# Load data: 
sn = 2
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
data_clust_folder <- glue("Data_{n_cluster}Cluster")
data_sample_folder <- glue("5260_Patients_Sample{sn}")
df <- read.csv2(paste0(root_path,"/",data_clust_folder,"/5260_Global_Indices_df",
                       glue("/df_5260_GI_campimeters{sn}.csv")))
rownames(df) <- df$X
df$X <- NULL
df <- subset(df, select = -c(MD_T, sLV_T)) # Variables: MD and sLV

# Data preparation: (NO) Standardize the data
data_train <- df

# Assign dataset of validation to clusters using cluster model
predict.kmeans <- function(object, newdata){
        centers <- object$centers
        n_centers <- nrow(centers)
        dist_mat <- as.matrix(dist(rbind(centers, newdata)))
        dist_mat <- dist_mat[-seq(n_centers), seq(n_centers)]
        max.col(-dist_mat)
}
km_predict <- predict.kmeans(clust,data_train)
data_train['Cluster'] <- as.vector(km_predict)
camps_train_c1 <- row.names(data_train[which(data_train$Cluster==1),])
camps_train_c2 <- row.names(data_train[which(data_train$Cluster==2),])
camps_train_c3 <- row.names(data_train[which(data_train$Cluster==3),])
camps_train_c4 <- row.names(data_train[which(data_train$Cluster==4),])
camps_train_c5 <- row.names(data_train[which(data_train$Cluster==5),])
camps_train_c6 <- row.names(data_train[which(data_train$Cluster==6),])
camps_train_c7 <- row.names(data_train[which(data_train$Cluster==7),])
camps_train_c8 <- row.names(data_train[which(data_train$Cluster==8),])
camps_train_c9 <- row.names(data_train[which(data_train$Cluster==9),])
camps_train_c10 <- row.names(data_train[which(data_train$Cluster==10),])
camps_train_c11 <- row.names(data_train[which(data_train$Cluster==11),])
camps_train_c12 <- row.names(data_train[which(data_train$Cluster==12),])

data_values_train <- read.csv2(paste0(root_path,"/",data_clust_folder,"/5260_Values_df",
                                    glue("/df_5260_Values_campimeters{sn}.csv")))
rownames(data_values_train) <- data_values_train$X
data_values_train$X <- NULL

data_train_c1 <- data_values_train[camps_train_c1,]
data_train_c2 <- data_values_train[camps_train_c2,]
data_train_c3 <- data_values_train[camps_train_c3,]
data_train_c4 <- data_values_train[camps_train_c4,]
data_train_c5 <- data_values_train[camps_train_c5,]
data_train_c6 <- data_values_train[camps_train_c6,]
data_train_c7 <- data_values_train[camps_train_c7,]
data_train_c8 <- data_values_train[camps_train_c8,]
data_train_c9 <- data_values_train[camps_train_c9,]
data_train_c10 <- data_values_train[camps_train_c10,]
data_train_c11 <- data_values_train[camps_train_c11,]
data_train_c12 <- data_values_train[camps_train_c12,]

# Calculate statistics for all patient clusters
# Create list of data clusters
cluster_list_train <- c()      # List of data clusters
for (n in (1:n_clust_p)){
        data <- paste0('data_train_c',n)
        cluster_list_train <- append(cluster_list_train,data)
}
cluster_list_train <- c(cluster_list_train)
data_clust_stats <- stats_points_all(cluster_list_train)
eye_plot_all(data_clust_stats,'mean','Defecto medio por cluster',cols,cuts)

# 4.1 Display some prediction examples

path_train <- paste0(root_path,"/",data_clust_folder,"/",data_sample_folder)
file_names_train <- file_path_sans_ext(list.files(path_train, pattern = "*.csv"))
set.seed(444) 
#n_ejm = 20
#random_ejms <- data_train[sample(nrow(data_train), n_ejm), ] 
#random_ejms_names <- row.names(random_ejms)
#prediction <- random_ejms$Cluster 
#ejm_idx <- match(random_ejms_names,file_names_train)
n_ejm_c <- 5
random_ejms_grouped <- data_train %>% rownames_to_column('Camp') %>% group_by(Cluster) %>% sample_n(size=n_ejm_c)
random_ejms_names <- random_ejms_grouped$Camp
prediction <- random_ejms_grouped$Cluster
ejm_idx <- match(random_ejms_names,file_names_train)
n_ejm <- n_ejm_c * n_clust_p

# Merge all sample values
data_all <- data.frame(coordinates)
coords <- coordinates  # From dummy campimeter
j <- 1
for (i in file_names_train[ejm_idx]){
        name <- i
        filepath <- file.path(file.path(paste0(path_train,"/",paste0(format(i),".csv"))))
        name <- file_path_sans_ext(basename(filepath))
        camp <- read.csv2(filepath, header = TRUE, sep = ";", skipNul = TRUE)
        camp_VNorm <- as.data.frame(as.integer(camp$VNorm))
        colnames(camp_VNorm) <- paste0("VNorm",j)
        #colnames(camp_VNorm) <- name
        data_all <- cbind(data_all,camp_VNorm)
        j <- j + 1
}
data_ejm_sp <- SpatialPointsDataFrame(data_all[,1:2], data_all, coords.nrs = numeric(0), 
                                      proj4string = CRS(as.character(NA)), bbox = NULL)
#View(data_ejm_sp@data)
# Plot parameters 
scale.parameter = 1.1
original.bbox <- eye_clusters@bbox
xshift = 0  # Shift to right in map units. 
yshift = 0  # Shift to left in map units.
edges = original.bbox
edges[1,] <- (edges[1,] - mean(edges[1,])) * scale.parameter + mean(edges[1,]) + xshift
edges[2,] <- (edges[2,] - mean(edges[2,])) * scale.parameter + mean(edges[2,]) + yshift
seq_names_ejm <- paste(paste("Ejm",rep(seq(1:n_ejm_c),n_clust_p/2)),"-")
panel_names <- paste0(paste0(paste(seq_names_ejm,paste0("C",prediction[1:(n_ejm/2)]))))
spplot(data_ejm_sp[,1:((n_ejm/2)+2)], paste0("VNorm",1:(n_ejm/2)), names.attr=panel_names,
       cex = 0.7, do.log = F,
       main=list(label="Predicción de campimetrías",cex=1.2),
       col.regions=cols,cuts=cuts,
       par.settings = list(fontsize = list(text = 8)),
       #color.key=FALSE,
       sp.layout = list("sp.polygons",eye_clusters,col="dark gray"),
       xlim = edges[1,], ylim = edges[2,],
       key.space="right",
       as.table = TRUE)
panel_names <- paste0(paste0(paste(seq_names_ejm,paste0("C",prediction[((n_ejm/2)+1):n_ejm]))))
spplot(data_ejm_sp[,c(1,2,c(((n_ejm/2)+3):(n_ejm+2)))], paste0("VNorm",((n_ejm/2)+1):(n_ejm)), names.attr=panel_names,
       cex = 0.7, do.log = F,
       main=list(label="Predicción de campimetrías",cex=1.2),
       col.regions=cols,cuts=cuts,
       par.settings = list(fontsize = list(text = 8)),
       #color.key=TRUE,
       sp.layout = list("sp.polygons",eye_clusters,col="dark gray"),
       xlim = edges[1,], ylim = edges[2,],
       key.space="right",
       as.table = TRUE)


#########################################################################################################
#     5. Fitting 'Average' variograms models   (Training data set)                                      #                                                                       #
#########################################################################################################

# Training data set (Sample 2)
sn = 2
root_path_train <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
data_clust_folder <- glue("Data_{n_cluster}Cluster")
path_train <- paste0(root_path_train,"/",data_clust_folder,glue("/5260_Patients_Sample{sn}"))

# 5.1 Compute means variograms

# Dummy variogram (Extract lags and np)
file_names_train <- file_path_sans_ext(list.files(path_train, pattern = "*.csv"))
path_dummy <- file.path(file.path(path_train,paste0(file_names_train[1],".csv")))
dummy <- read.csv2(path_dummy, header = TRUE, sep = ";", skipNul = FALSE)
coordinates <- dummy[,1:2]
dummy_sp <- SpatialPointsDataFrame(coordinates, dummy, coords.nrs = numeric(0), 
                                   proj4string = CRS(as.character(NA)), bbox = NULL)
dummy_v <- variogram(VNorm ~ 1 , dummy_sp, boundaries=c(1,2,3,4,5,6,7,8,9)) #Boundaries
lags <- dummy_v$dist
np <- dummy_v$np

# Calculate Mean of the values in each distance lag. 
# One point for each lag. Do it for each cluster
cluster_list_train <- c(cluster_list_train)
means_variogram_clusts <- data.frame(lags)
clust_names <- paste0('C',1:length(cluster_list_train))
col_names <- append(c('Distance'), clust_names)

for (c in cluster_list_train){
        variograms_values <- vector()
        file_names_c <- row.names(get(c))
        for(i in file_names_c){
                filepath <- file.path(file.path(path_train,paste0(format(i),".csv")))
                camp <- read.csv2(filepath, header = TRUE, sep = ";", skipNul = TRUE)
                coordinates(camp) = ~ x + y
                v_i <- variogram(VNorm ~ 1 , camp, boundaries=c(1,2,3,4,5,6,7,8,9))
                variograms_values <- rbind(variograms_values,as.vector(v_i[,3]))
                #variance <- colVars(variograms_values)
                #std <- colStdevs(variograms_values)
                means <- data.frame(colMeans(variograms_values))
                #medianas <- colMedians(variograms_values)
                #sum <- colSums(variograms_values)
        }
        means_variogram_clusts <- cbind(means_variogram_clusts,means)
}
colnames(means_variogram_clusts) <- col_names

# Create variogram class for each cluster (assign into variables)

for(i in names(means_variogram_clusts)[2:length(means_variogram_clusts)]){
        media.vgm = data.frame(dist=lags,gamma=means_variogram_clusts[[i]])
        media.vgm$np = np #Number of point pairs (Calculado antes)
        media.vgm$dir.hor = rep(0,length(np))
        media.vgm$dir.ver = rep(0,length(np))
        class(media.vgm) = c("gstatVariogram","data.frame")
        assign(paste0('media.vgm.',i),media.vgm)
}

cluster_list_train_var <- c()      # List of mean variograms by cluster
for (n in (1:n_clust_p)){
        data <- paste0('media.vgm.C',n)
        cluster_list_train_var <- append(cluster_list_train_var,data)
}


# 5.2 Fit variograms models 

# Exponencial
models_errors <- c()
n <- 1
for (i in cluster_list_train_var){
        var <- get(i)
        vg.exp <- vgm(psill=max(var$gamma)*0.9, model = "Exp", range=max(var$dist)/2,
                      nugget = mean(var$gamma)/4)
        fit.media.exp <- fit.variogram(var, vg.exp,fit.method = 7)
        assign(paste0('fit.media.exp.C',n),fit.media.exp)
        if ((fit.media.exp$range[1] < 0) | (fit.media.exp$range[2] < 0)){
                models_errors <- c(models_errors,i)
                }
        n <- n+1
}
models_errors


# # VF = 8 -> Arreglar fit.media.exp.C1
# var <- get(models_errors[1])
# vg.exp <- vgm(psill=max(var$gamma)*0.7, model = "Exp", range=max(var$dist))
# fit.media.exp.C1 <- fit.variogram(var, vg.exp,fit.method = 7)

# VF = 10 -> Arreglar fit.media.exp.C12
var <- get(models_errors[1])
vg.exp <- vgm(psill=max(var$gamma)*1, model = "Exp", range=max(var$dist)/4)
fit.media.exp.C12 <- fit.variogram(var, vg.exp,fit.method = 7)
plot(media.vgm.C12,fit.media.exp.C12, col = "red",main = "Mean semivariogram", sub= "Exponential model")


# Spheric
models_errors <- c()
n <- 1
for (i in cluster_list_train_var){
        var <- get(i)
        vg.sph <- vgm(psill=max(var$gamma)*0.9, model = "Sph", range=max(var$dist)/2,
                      nugget = mean(var$gamma)/4)
        fit.media.sph <- fit.variogram(var, vg.sph,fit.method = 7)
        assign(paste0('fit.media.sph.C',n),fit.media.sph)
        if ((fit.media.sph$range[1] < 0) | (fit.media.sph$range[2] < 0)){
                models_errors <- c(models_errors,i)
        }
        n <- n+1
}
models_errors
#plot(media.vgm.C3,fit.media.sph.C3, col = "red",main = "Mean semivariogram", sub= "Spheric model")

# Gaussian
models_errors <- c()
n <- 1
for (i in cluster_list_train_var){
        var <- get(i)
        vg.gau <- vgm(psill=max(var$gamma)*0.9, model = "Gau", range=max(var$dist)/2,
                      nugget = mean(var$gamma)/4)
        fit.media.gau <- fit.variogram(var, vg.gau,fit.method = 7)
        assign(paste0('fit.media.gau.C',n),fit.media.gau)
        if ((fit.media.gau$range[1] < 0) | (fit.media.gau$range[2] < 0)){
                models_errors <- c(models_errors,i)
        }
        n <- n+1
}
models_errors

# # VF 8 -> Arreglar C1, C2, C9 y C10
# var <- get(models_errors[1])
# vg.gau <- vgm(psill=max(var$gamma)*0.9, model = "Gau", range=max(var$dist)/4)
# fit.media.gau.C1 <- fit.variogram(var, vg.gau,fit.method = 6)
# plot(media.vgm.C1,fit.media.gau.C1, col = "red",main = "Mean semivariogram", sub= "Gaussian model")
# var <- get(models_errors[2])
# vg.gau <- vgm(psill=max(var$gamma)*1.2, model = "Gau", range=max(var$dist)/4, nugget=2000)
# fit.media.gau.C2 <- fit.variogram(var, vg.gau,fit.method = 6)
# plot(media.vgm.C2,fit.media.gau.C2, col = "red",main = "Mean semivariogram", sub= "Gaussian model")
# var <- get(models_errors[3])
# vg.gau <- vgm(psill=max(var$gamma)*1, model = "Gau", range=max(var$dist)/4, nugget = 1500)
# fit.media.gau.C9 <- fit.variogram(var, vg.gau,fit.method = 6)
# plot(media.vgm.C9,fit.media.gau.C9, col = "red",main = "Mean semivariogram", sub= "Gaussian model")
# var <- get(models_errors[4])
# vg.gau <- vgm(psill=max(var$gamma)*1, model = "Gau", range=max(var$dist)/4)
# fit.media.gau.C10 <- fit.variogram(var, vg.gau,fit.method = 6)
# plot(media.vgm.C10,fit.media.gau.C10, col = "red",main = "Mean semivariogram", sub= "Gaussian model")

# VF 10 -> Arreglar C4, C11, C12
var <- get(models_errors[1])
vg.gau <- vgm(psill=max(var$gamma)*0.8, model = "Gau", range=max(var$dist)/2,nugget=500)
fit.media.gau.C4 <- fit.variogram(var, vg.gau,fit.method = 7)
plot(media.vgm.C4,fit.media.gau.C4, col = "red",main = "Mean semivariogram", sub= "Gaussian model")
var <- get(models_errors[2])
vg.gau <- vgm(psill=max(var$gamma)*1.2, model = "Gau", range=max(var$dist)/4, nugget=2000)
fit.media.gau.C11 <- fit.variogram(var, vg.gau,fit.method = 6)
plot(media.vgm.C11,fit.media.gau.C11, col = "red",main = "Mean semivariogram", sub= "Gaussian model")
var <- get(models_errors[3])
vg.gau <- vgm(psill=max(var$gamma)*1, model = "Gau", range=max(var$dist)/4, nugget = 1000)
fit.media.gau.C12 <- fit.variogram(var, vg.gau,fit.method = 7)
plot(media.vgm.C12,fit.media.gau.C12, col = "red",main = "Mean semivariogram", sub= "Gaussian model")


# Lineal
models_errors <- c()
n <- 1
for (i in cluster_list_train_var){
        var <- get(i)
        vg.lin <- vgm(psill=max(var$gamma)*0.9, model = "Lin", range=max(var$dist)/2,
                      nugget = mean(var$gamma)/4)
        fit.media.lin <- fit.variogram(var, vg.lin,fit.method = 7)
        assign(paste0('fit.media.lin.C',n),fit.media.lin)
        if ((fit.media.lin$range[1] < 0) | (fit.media.lin$range[2] < 0)){
                models_errors <- c(models_errors,i)
        }
        n <- n+1
}
models_errors

# Circular
models_errors <- c()
n <- 1
for (i in cluster_list_train_var){
        var <- get(i)
        vg.cir <- vgm(psill=max(var$gamma)*0.9, model = "Cir", range=max(var$dist)/2, 
                      nugget = mean(var$gamma)/4)
        fit.media.cir <- fit.variogram(var, vg.cir,fit.method = 7)
        assign(paste0('fit.media.cir.C',n),fit.media.cir)
        if ((fit.media.cir$range[1] < 0) | (fit.media.cir$range[2] < 0)){
                models_errors <- c(models_errors,i)
        }
        n <- n+1
}
models_errors
#plot(media.vgm.C2,fit.media.cir.C2, main = "Mean semivariogram", sub= "Circular model")

#save.image(file=paste0(root_path,"/Global_models_VF10_unique.RData"))


#########################################################################################################
#     6. Kriging CV (Validation data set)                                                               #     
#########################################################################################################
root_path = "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
load(paste0(root_path,"/Global_models_VF8_unique.RData"))

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

# 6.1 Assign validation data set into clusters (Sample 3)

set.seed(123)
# Load data: 
sn = 3
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
data_clust_folder <- glue("Data_{n_cluster}Cluster")
data_sample_folder <- glue("5260_Patients_Sample{sn}")
path_sample <- paste0(root_path,"/",data_clust_folder,"/",data_sample_folder)
df <- read.csv2(paste0(root_path,"/",data_clust_folder,"/5260_Global_Indices_df",
                       glue("/df_5260_GI_campimeters{sn}.csv")))
rownames(df) <- df$X
df$X <- NULL
df <- subset(df, select = -c(MD_T, sLV_T)) # Variables: MD and sLV

# Data preparation: (NO) Standardize the data
data_val <- df

km_predict <- predict.kmeans(clust,data_val)
data_val['Cluster'] <- as.vector(km_predict)
camps_val_c1 <- row.names(data_val[which(data_val$Cluster==1),])
camps_val_c2 <- row.names(data_val[which(data_val$Cluster==2),])
camps_val_c3 <- row.names(data_val[which(data_val$Cluster==3),])
camps_val_c4 <- row.names(data_val[which(data_val$Cluster==4),])
camps_val_c5 <- row.names(data_val[which(data_val$Cluster==5),])
camps_val_c6 <- row.names(data_val[which(data_val$Cluster==6),])
camps_val_c7 <- row.names(data_val[which(data_val$Cluster==7),])
camps_val_c8 <- row.names(data_val[which(data_val$Cluster==8),])
camps_val_c9 <- row.names(data_val[which(data_val$Cluster==9),])
camps_val_c10 <- row.names(data_val[which(data_val$Cluster==10),])
camps_val_c11 <- row.names(data_val[which(data_val$Cluster==11),])
camps_val_c12 <- row.names(data_val[which(data_val$Cluster==12),])

data_values_val <- read.csv2(paste0(root_path,"/",data_clust_folder,"/5260_Values_df",
                                       glue("/df_5260_Values_campimeters{sn}.csv")))
rownames(data_values_val) <- data_values_val$X
data_values_val$X <- NULL

data_val_c1 <- data_values_train[camps_val_c1,]
data_val_c2 <- data_values_train[camps_val_c2,]
data_val_c3 <- data_values_train[camps_val_c3,]
data_val_c4 <- data_values_train[camps_val_c4,]
data_val_c5 <- data_values_train[camps_val_c5,]
data_val_c6 <- data_values_train[camps_val_c6,]
data_val_c7 <- data_values_train[camps_val_c7,]
data_val_c8 <- data_values_train[camps_val_c8,]
data_val_c9 <- data_values_train[camps_val_c9,]
data_val_c10 <- data_values_train[camps_val_c10,]
data_val_c11 <- data_values_train[camps_val_c11,]
data_val_c12 <- data_values_train[camps_val_c12,]

# Calculate statistics for all patient clusters
# Create list of data clusters
cluster_list_val <- c()      # List of data clusters
for (n in (1:n_clust_p)){
        data <- paste0('data_val_c',n)
        cluster_list_val <- append(cluster_list_val,data)
}

# Create list of fitted models 
cluster_list_models_exp <- c()  
cluster_list_models_sph <- c() 
cluster_list_models_gau <- c()  
cluster_list_models_lin <- c()  
cluster_list_models_cir <- c()   


for (n in (1:n_clust_p)){
        model_exp <- paste0('fit.media.exp.C',n)
        cluster_list_models_exp <- append(cluster_list_models_exp,model_exp)
        model_sph <- paste0('fit.media.sph.C',n)
        cluster_list_models_sph <- append(cluster_list_models_sph,model_sph)
        model_gau <- paste0('fit.media.gau.C',n)
        cluster_list_models_gau <- append(cluster_list_models_gau,model_gau)
        model_lin <- paste0('fit.media.lin.C',n)
        cluster_list_models_lin <- append(cluster_list_models_lin,model_lin)
        model_cir <- paste0('fit.media.cir.C',n)
        cluster_list_models_cir <- append(cluster_list_models_cir,model_cir)
}
cluster_list_models_exp <- c(cluster_list_models_exp)
cluster_list_models_sph <- c(cluster_list_models_sph)
cluster_list_models_gau <- c(cluster_list_models_gau)
cluster_list_models_lin <- c(cluster_list_models_lin)
cluster_list_models_cir <- c(cluster_list_models_cir)

studied_models <- c('cluster_list_models_exp','cluster_list_models_sph',
                    'cluster_list_models_gau','cluster_list_models_lin',
                    'cluster_list_models_cir')


# Save CV Kriging by model and their respective cluster validation data set
path_results <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/RESULTS"
for (j in studied_models){
        name_model <- stri_sub(j,-3,-1)
        cluster_list_models <- get(j)
        c <- 1
        for (m in cluster_list_models){
                model_c <- get(m)
                file_names_c <- get(paste0('camps_val_c',c))
                krigeCV <- krig_results_by_neighbours(file_names_c,model_c,seq(2,8))
                write.csv2(krigeCV[[1]], file.path(path_results,glue("media_{name_model}_results_C{c}_VF{n_cluster}.csv")))
                write.csv2(krigeCV[[2]], file.path(path_results,glue("media_{name_model}_errors_C{c}_VF{n_cluster}.csv")))
                c <- c+1
        }
}


#########################################################################################################
#     7. Compare variogram models (Kriging CV - Validation dataset)                                     #                                                                       #
#########################################################################################################
n_cluster <- 8 # Number of Visual Field zones
cluster <- 12 # Select patient's cluster
# Read data:
path_errors <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/RESULTS"
folder_VF <- glue("{n_cluster}VF_kriging")
folder_cluster <- glue("C{cluster}")
path_cluster <- file.path(path_errors,folder_VF,folder_cluster)
path_errors <- path_cluster

#media_lin_results <- read.csv2(paste0(path_errors,glue("/media_lin_results_C{cluster}_VF{n_cluster}.csv")))
media_lin_errors <- read.csv2(paste0(path_errors,glue("/media_lin_errors_C{cluster}_VF{n_cluster}.csv")))
#media_gauss_results <- read.csv2(paste0(path_errors,glue("/media_gau_results_C{cluster}_VF{n_cluster}.csv")))
media_gauss_errors <- read.csv2(paste0(path_errors,glue("/media_gau_errors_C{cluster}_VF{n_cluster}.csv")))
#media_cir_results <- read.csv2(paste0(path_errors,glue("/media_cir_results_C{cluster}_VF{n_cluster}.csv")))
media_cir_errors <- read.csv2(paste0(path_errors,glue("/media_cir_errors_C{cluster}_VF{n_cluster}.csv")))
#media_sph_results <- read.csv2(paste0(path_errors,glue("/media_sph_results_C{cluster}_VF{n_cluster}.csv")))
media_sph_errors <- read.csv2(paste0(path_errors,glue("/media_sph_errors_C{cluster}_VF{n_cluster}.csv")))
#media_exp_results <- read.csv2(paste0(path_errors,glue("/media_exp_results_C{cluster}_VF{n_cluster}.csv")))
media_exp_errors <- read.csv2(paste0(path_errors,glue("/media_exp_errors_C{cluster}_VF{n_cluster}.csv")))

# 7.1 Mean of error measurements by campimeter and by total

# 7.1.1 Mean of error measurements by campimetry
media_lin_errors_grouped    <- media_lin_errors %>%
        group_by(neighbours_list) %>%
        summarize_all(mean)
media_gauss_errors_grouped  <- media_gauss_errors %>%
        group_by(neighbours_list) %>%
        summarize_all(mean)
media_cir_errors_grouped  <- media_cir_errors %>%
        group_by(neighbours_list) %>%
        summarize_all(mean)
media_sph_errors_grouped  <- media_sph_errors %>%
        group_by(neighbours_list) %>%
        summarize_all(mean)
media_exp_errors_grouped  <- media_exp_errors %>%
        group_by(neighbours_list) %>%
        summarize_all(mean)

# 7.1.2 Mean of error measurements by total
media_lin_results_grouped    <- media_lin_results %>%
        group_by(n_neighbours) %>%
        summarize_all(c(mean,var))
media_gauss_results_grouped  <- media_gauss_results %>%
        group_by(n_neighbours) %>%
        summarize_all(c(mean,var))
media_cir_results_grouped  <- media_cir_results %>%
        group_by(n_neighbours) %>%
        summarize_all(c(mean,var))
media_sph_results_grouped  <- media_sph_results %>%
        group_by(n_neighbours) %>%
        summarize_all(c(mean,var))
media_exp_results_grouped  <- media_exp_results %>%
        group_by(n_neighbours) %>%
        summarize_all(c(mean,var))

# 7.2. Display of error measurements

# 7.2.1 Evolution of error measures with respect to the number of neighbors 

# Colores <- brewer.pal(10,"Paired") # Paleta Paired
# list_col <- c(Colores[1],Colores[2],Colores[3],Colores[7],Colores[5]) #Todos
# dev.off()
# par(mfcol = c(2, 3))
# # Media de la media de los residuos por campimetría
# plot(media_exp_errors_grouped$neighbours_list,media_exp_errors_grouped$mean_res,type="l",
#      col=Colores[1], xlab= "Número de vecinos", ylab="Error medio",
#      ylim=c(-0.5,1.5), main="Media del error medio", cex.main=0.8)
# lines(media_sph_errors_grouped$neighbours_list,media_sph_errors_grouped$mean_res,type="l",
#       col=Colores[2], xlab= "Número de vecinos", ylab="Error medio")
# lines(media_gauss_errors_grouped$neighbours_list,media_gauss_errors_grouped$mean_res,type="l",
#       col=Colores[3], xlab= "Número de vecinos", ylab="Error medio")
# lines(media_lin_errors_grouped$neighbours_list,media_lin_errors_grouped$mean_res,type="l",
#       col=Colores[7], xlab= "Número de vecinos", ylab="Error medio")
# lines(media_cir_errors_grouped$neighbours_list,media_cir_errors_grouped$mean_res,type="l",
#       col=Colores[5], xlab= "Número de vecinos", ylab="Error medio")
# legend(2, 1.49, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
#                           "Modelo gausiano de medias", "Modelo lineal de medias",
#                           "Modelo circular de medias"),
#        col=list_col, lty=1, cex=0.3, box.lty=0)
# 
# 
# # # Media de la media de los residuos normalizados por campimetría
# # plot(media_exp_errors_grouped$neighbours_list,media_exp_errors_grouped$mean_zscore,type="l",
# #      col=Colores[1], xlab= "Número de vecinos", ylab="Error medio normalizado",
# #      ylim=c(-0.4,0.4), main="Media del error medio normalizado", cex.main=0.8)
# # lines(media_sph_errors_grouped$neighbours_list,media_sph_errors_grouped$mean_zscore,type="l",
# #       col=Colores[2], xlab= "Número de vecinos", ylab="Error medio")
# # lines(media_gauss_errors_grouped$neighbours_list,media_gauss_errors_grouped$mean_zscore,type="l",
# #       col=Colores[3], xlab= "Número de vecinos", ylab="Error medio")
# # lines(media_lin_errors_grouped$neighbours_list,media_lin_errors_grouped$mean_zscore,type="l",
# #       col=Colores[7], xlab= "Número de vecinos", ylab="Error medio")
# # lines(media_cir_errors_grouped$neighbours_list,media_cir_errors_grouped$mean_zscore,type="l",
# #       col=Colores[5], xlab= "Número de vecinos", ylab="Error medio")
# # legend(2, -0.25, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
# #                             "Modelo gausiano de medias", "Modelo lineal de medias",
# #                             "Modelo circular de medias"),
# #        col=list_col, lty=1, cex=0.6)
# 
# 
# # Media de los errores cuadráticos medios por campimetría
# plot(media_exp_errors_grouped$neighbours_list,media_exp_errors_grouped$mspe,type="l",
#      col=Colores[1], xlab= "Número de vecinos", ylab="Error cuadrático medio",
#      ylim=c(500,700), main="Media del error cuadrático medio", cex.main=0.8)
# lines(media_sph_errors_grouped$neighbours_list,media_sph_errors_grouped$mspe,type="l",
#       col=Colores[2], xlab= "Número de vecinos", ylab="Error cuadrático medio")
# lines(media_gauss_errors_grouped$neighbours_list,media_gauss_errors_grouped$mspe,type="l",
#       col=Colores[3], xlab= "Número de vecinos", ylab="Error cuadrático medio")
# lines(media_lin_errors_grouped$neighbours_list,media_lin_errors_grouped$mspe,type="l",
#       col=Colores[7], xlab= "Número de vecinos", ylab="Error cuadrático medio")
# lines(media_cir_errors_grouped$neighbours_list,media_cir_errors_grouped$mspe,type="l",
#       col=Colores[5], xlab= "Número de vecinos", ylab="Error cuadrático medio")
# legend(2, 699, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
#                           "Modelo gausiano de medias", "Modelo lineal de medias",
#                           "Modelo circular de medias"),
#        col=list_col, lty=1, cex=0.3, box.lty=0)
# 
# # Media de los errores cuadráticos normalizados medios por campimetría
# plot(media_exp_errors_grouped$neighbours_list,media_exp_errors_grouped$mspe_norm,type="l",
#      col=Colores[1], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado",
#      ylim=c(0.6,1), main="Media del error cuadrático medio normalizado", cex.main=0.8)
# lines(media_sph_errors_grouped$neighbours_list,media_sph_errors_grouped$mspe_norm,type="l",
#       col=Colores[2], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado")
# lines(media_gauss_errors_grouped$neighbours_list,media_gauss_errors_grouped$mspe_norm,type="l",
#       col=Colores[3], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado")
# lines(media_lin_errors_grouped$neighbours_list,media_lin_errors_grouped$mspe_norm,type="l",
#       col=Colores[7], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado")
# lines(media_cir_errors_grouped$neighbours_list,media_cir_errors_grouped$mspe_norm,type="l",
#       col=Colores[5], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado")
# legend(2, 0.99, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
#                          "Modelo gausiano de medias", "Modelo lineal de medias",
#                          "Modelo circular de medias"),
#        col=list_col, lty=1, cex=0.3, box.lty=0)
# 
# # Media de la correlación entre valores observados y predichos
# plot(media_exp_errors_grouped$neighbours_list,media_exp_errors_grouped$cor_obs_pred,type="l",
#      col=Colores[1], xlab= "Número de vecinos", ylab="Correlación observado-predicho",
#      ylim=c(0.7,0.75), main="Media de la correlación observado-predicho", cex.main=0.8)
# lines(media_sph_errors_grouped$neighbours_list,media_sph_errors_grouped$cor_obs_pred,type="l",
#       col=Colores[2], xlab= "Número de vecinos", ylab="Correlación observado-predicho")
# lines(media_gauss_errors_grouped$neighbours_list,media_gauss_errors_grouped$cor_obs_pred,type="l",
#       col=Colores[3], xlab= "Número de vecinos", ylab="Correlación observado-predicho")
# lines(media_lin_errors_grouped$neighbours_list,media_lin_errors_grouped$cor_obs_pred,type="l",
#       col=Colores[7], xlab= "Número de vecinos", ylab="Correlación observado-predicho")
# lines(media_cir_errors_grouped$neighbours_list,media_cir_errors_grouped$cor_obs_pred,type="l",
#       col=Colores[5], xlab= "Número de vecinos", ylab="Correlación observado-predicho")
# legend(2, 0.74999, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
#                          "Modelo gausiano de medias", "Modelo lineal de medias",
#                          "Modelo circular de medias"),
#        col=list_col, lty=1, cex=0.3, box.lty=0)
# 
# # Media de la correlación entre valores predichos y errores normalizados
# plot(media_exp_errors_grouped$neighbours_list,media_exp_errors_grouped$cor_pred_res_n,type="l",
#      col=Colores[1], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho",
#      ylim=c(-0.2,0.25), main="Media de la correlación error normalizado-predicho", cex.main=0.8)
# lines(media_sph_errors_grouped$neighbours_list,media_sph_errors_grouped$cor_pred_res_n,type="l",
#       col=Colores[2], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho")
# lines(media_gauss_errors_grouped$neighbours_list,media_gauss_errors_grouped$cor_pred_res_n,type="l",
#       col=Colores[3], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho")
# lines(media_lin_errors_grouped$neighbours_list,media_lin_errors_grouped$cor_pred_res_n,type="l",
#       col=Colores[7], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho")
# lines(media_cir_errors_grouped$neighbours_list,media_cir_errors_grouped$cor_pred_res_n,type="l",
#       col=Colores[5], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho")
# legend(2, 0.249, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
#                           "Modelo gausiano de medias", "Modelo lineal de medias",
#                           "Modelo circular de medias"),
#        col=list_col, lty=1, cex=0.3, box.lty=0)
# 
# # Proportion of normalized errors above 2.5
# plot(media_exp_errors_grouped$neighbours_list,media_exp_errors_grouped$percent,type="l",
#      col=Colores[1], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5",
#      ylim=c(0.75,1.5), main="Media del porcentaje errores normalizados > 2.5", cex.main=0.8)
# lines(media_sph_errors_grouped$neighbours_list,media_sph_errors_grouped$percent,type="l",
#       col=Colores[2], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5")
# lines(media_gauss_errors_grouped$neighbours_list,media_gauss_errors_grouped$percent,type="l",
#       col=Colores[3], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5")
# lines(media_lin_errors_grouped$neighbours_list,media_lin_errors_grouped$percent,type="l",
#       col=Colores[7], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5")
# lines(media_cir_errors_grouped$neighbours_list,media_cir_errors_grouped$percent,type="l",
#       col=Colores[5], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5")
# legend(5.5, 0.9, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
#                         "Modelo gausiano de medias", "Modelo lineal de medias",
#                         "Modelo circular de medias"),
#        col=list_col, lty=1, cex=0.3, box.lty=0)
# 
# mtext(glue("Medidas de error por número de vecinos: C{cluster}"), side = 3, line = -1.5, outer = TRUE, cex=0.8)

# 7.2.2 Table of mean error measurements by model. Fixed n 

n = 5

# -Media linear
#media_lin_res_n <- media_lin_results[which(media_lin_results$n_neighbours==n),]
media_lin_err_n <- media_lin_errors[which(media_lin_errors$neighbours_list==n),]
# -Media gaussiano
#media_gauss_res_n <- media_gauss_results[which(media_gauss_results$n_neighbours==n),]
media_gauss_err_n <- media_gauss_errors[which(media_gauss_errors$neighbours_list==n),]
# -Media circular
#media_cir_res_n <- media_cir_results[which(media_cir_results$n_neighbours==n),]
media_cir_err_n <- media_cir_errors[which(media_cir_errors$neighbours_list==n),]
# -Media esférico
#media_sph_res_n <- media_sph_results[which(media_sph_results$n_neighbours==n),]
media_sph_err_n <- media_sph_errors[which(media_sph_errors$neighbours_list==n),]
# -Media exponencial
#media_exp_res_n <- media_exp_results[which(media_exp_results$n_neighbours==n),]
media_exp_err_n <- media_exp_errors[which(media_exp_errors$neighbours_list==n),]

errores_medios <- c(round(mean(media_exp_err_n$mean_res),4),round(mean(media_sph_err_n$mean_res),4),
                    round(mean(media_gauss_err_n$mean_res),4),round(mean(media_lin_err_n$mean_res),4),
                    round(mean(media_cir_err_n$mean_res),4))
errores_cuad_medios <- c(round(mean(media_exp_err_n$mspe),4),round(mean(media_sph_err_n$mspe),4),
                         round(mean(media_gauss_err_n$mspe),4),round(mean(media_lin_err_n$mspe),4),
                         round(mean(media_cir_err_n$mspe),4))
errores_cuad_medios_n <- c(round(mean(media_exp_err_n$mspe_norm),4),round(mean(media_sph_err_n$mspe_norm),4),
                           round(mean(media_gauss_err_n$mspe_norm),4),round(mean(media_lin_err_n$mspe_norm),4),
                           round(mean(media_cir_err_n$mspe_norm),4))
errores_cor_obs_pred <- c(round(mean(media_exp_err_n$cor_obs_pred),4),round(mean(media_sph_err_n$cor_obs_pred),4),
                          round(mean(media_gauss_err_n$cor_obs_pred),4),round(mean(media_lin_err_n$cor_obs_pred),4),
                          round(mean(media_cir_err_n$cor_obs_pred),4))
errores_cor_pred_errn <- c(round(mean(media_exp_err_n$cor_pred_res_n),4),round(mean(media_sph_err_n$cor_pred_res_n),4),
                           round(mean(media_gauss_err_n$cor_pred_res_n),4),round(mean(media_lin_err_n$cor_pred_res_n),4),
                           round(mean(media_cir_err_n$cor_pred_res_n),4))
errores_porcentaje <- c(round(mean(media_exp_err_n$percent),4),round(mean(media_sph_err_n$percent),4),
                        round(mean(media_gauss_err_n$percent),4),round(mean(media_lin_err_n$percent),4),
                        round(mean(media_cir_err_n$percent),4))

errores_modelos <- matrix(c(errores_medios,errores_cuad_medios,errores_cuad_medios_n,
                            errores_cor_obs_pred,errores_cor_pred_errn,errores_porcentaje),
                          ncol=5,byrow=TRUE)
colnames(errores_modelos) <- c("Exponencial","Esférico","Gausiano","Lineal","Circular")
rownames(errores_modelos) <- c("Error medio","Error cuadrático medio","Error cuadrático medio normalizado",
                               "Correlación observado-predicho","Correlación predicho-error normalizado",
                               "Porcentaje errores normalizados >2.5")
errores_modelos <- as.table(errores_modelos)
errores_modelos

dev.off()
grid.table(errores_modelos)


#########################################################################################################
#     8. Kriging CV (Test data set)                                                                     #     
#########################################################################################################
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
load(paste0(root_path,"/Global_models_VF8_unique.RData"))

# 8.1 Assign test data set into cluster
set.seed(123)
# Load data: 4 Sample
sn = 4
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
data_clust_folder <- glue("Data_{n_cluster}Cluster")
data_sample_folder <- glue("5260_Patients_Sample{sn}")
path_sample <- file.path(root_path,data_clust_folder,data_sample_folder)
df <- read.csv2(paste0(root_path,"/",data_clust_folder,"/5260_Global_Indices_df",
                        glue("/df_5260_GI_campimeters{sn}.csv")))
rownames(df) <- df$X
df$X <- NULL
df <- subset(df, select = -c(MD_T, sLV_T)) # Variables: MD and sLV

# Data preparation: (NO) Standardize the data
data_test <- df

# Assign dataset of validation to clusters using cluster model
predict.kmeans <- function(object, newdata){
        centers <- object$centers
        n_centers <- nrow(centers)
        dist_mat <- as.matrix(dist(rbind(centers, newdata)))
        dist_mat <- dist_mat[-seq(n_centers), seq(n_centers)]
        max.col(-dist_mat)
}
km_predict <- predict.kmeans(clust,data_test)
data_test['Cluster'] <- as.vector(km_predict)
camps_test_c1 <- row.names(data_test[which(data_test$Cluster==1),])
camps_test_c2 <- row.names(data_test[which(data_test$Cluster==2),])
camps_test_c3 <- row.names(data_test[which(data_test$Cluster==3),])
camps_test_c4 <- row.names(data_test[which(data_test$Cluster==4),])
camps_test_c5 <- row.names(data_test[which(data_test$Cluster==5),])
camps_test_c6 <- row.names(data_test[which(data_test$Cluster==6),])
camps_test_c7 <- row.names(data_test[which(data_test$Cluster==7),])
camps_test_c8 <- row.names(data_test[which(data_test$Cluster==8),])
camps_test_c9 <- row.names(data_test[which(data_test$Cluster==9),])
camps_test_c10 <- row.names(data_test[which(data_test$Cluster==10),])
camps_test_c11 <- row.names(data_test[which(data_test$Cluster==11),])
camps_test_c12 <- row.names(data_test[which(data_test$Cluster==12),])

data_values_test <- read.csv2(paste0(root_path,"/",data_clust_folder,"/5260_Values_df",
                                       glue("/df_5260_Values_campimeters{sn}.csv")))
rownames(data_values_test) <- data_values_test$X
data_values_test$X <- NULL

data_test_c1 <- data_values_test[camps_test_c1,]
data_test_c2 <- data_values_test[camps_test_c2,]
data_test_c3 <- data_values_test[camps_test_c3,]
data_test_c4 <- data_values_test[camps_test_c4,]
data_test_c5 <- data_values_test[camps_test_c5,]
data_test_c6 <- data_values_test[camps_test_c6,]
data_test_c7 <- data_values_test[camps_test_c7,]
data_test_c8 <- data_values_test[camps_test_c8,]
data_test_c9 <- data_values_test[camps_test_c9,]
data_test_c10 <- data_values_test[camps_test_c10,]
data_test_c11 <- data_values_test[camps_test_c11,]
data_test_c12 <- data_values_test[camps_test_c12,]

# Calculate statistics for all patient clusters
# Create list of data clusters
cluster_list_test <- c()      # List of data clusters
for (n in (1:n_clust_p)){
        data <- paste0('data_test_c',n)
        cluster_list_test <- append(cluster_list_test,data)
}
cluster_list_test <- c(cluster_list_test)
data_clust_stats <- stats_points_all(cluster_list_test)
eye_plot_all(data_clust_stats,'mean','Defecto medio por cluster',cols,cuts)

# 8.2 Compute kriging CV results 

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


# Create list of fitted models

# selected_models_8 <- c("fit.media.cir.C1","fit.media.exp.C2","fit.media.gau.C3",
#                        "fit.media.gau.C4","fit.media.sph.C5","fit.media.exp.C6",
#                        "fit.media.exp.C7","fit.media.exp.C8","fit.media.exp.C9",
#                        "fit.media.exp.C10","fit.media.sph.C11","fit.media.exp.C12")
# n_neighbors_8 <- c(6,8,3,5,4,7,4,8,4,4,8,8)

selected_models_10 <- c("fit.media.sph.C1","fit.media.sph.C2","fit.media.gau.C3",
                        "fit.media.sph.C4","fit.media.gau.C5","fit.media.exp.C6",
                        "fit.media.exp.C7","fit.media.exp.C8","fit.media.exp.C9",
                        "fit.media.exp.C10","fit.media.lin.C11","fit.media.exp.C12")
n_neighbors_10 <- c(6,8,5,4,7,8,8,8,4,8,5,5)


# Save CV Kriging by model and their respective cluster test data set

cluster_list_models <- selected_models_10
n_neighbors <- n_neighbors_10
c <- 1
for (m in cluster_list_models){
        model_c <- get(m)
        file_names_c <- get(paste0('camps_test_c',c))
        krigeCV <- krige.results(file_names_c,model_c,n_neighbors[c])
        write.csv2(krigeCV[[1]], file.path(root_path,glue("media_results_C{c}_VF{n_cluster}_test.csv")))
        write.csv2(krigeCV[[2]], file.path(root_path,glue("media_errors_C{c}_VF{n_cluster}_test.csv")))
        c <- c+1
}


# 8.3 Display error measurements
n_cluster = 10
n_clust_p = 12

path_results <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/RESULTS"
results_folder <- glue("{n_cluster}VF_kriging")
path_results <- file.path(path_results,results_folder,"TEST")

for (i in 1:n_clust_p){
        errors <- read.csv2(file.path(path_results,glue("media_errors_C{i}_VF{n_cluster}_test.csv")), sep = ";")
        results <- read.csv2(file.path(path_results,glue("media_results_C{i}_VF{n_cluster}_test.csv")), sep = ";")
        assign(glue('errors_test_C{i}'),errors)
        assign(glue('results_test_C{i}'),results)
}

# Mean of all test campimeters errors by cluster

# errors_test_mean_C1 <- errors_test_C1 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C1)
# errors_test_mean_C2 <- errors_test_C2 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C2)
# errors_test_mean_C3 <- errors_test_C3 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C3)
# errors_test_mean_C4 <- errors_test_C4 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C4)
# errors_test_mean_C5 <- errors_test_C5 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C5)
# errors_test_mean_C6 <- errors_test_C6 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C6)
# errors_test_mean_C7 <- errors_test_C7 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C7)
# errors_test_mean_C8 <- errors_test_C8 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C8)
# errors_test_mean_C9 <- errors_test_C9 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C9)
# errors_test_mean_C10 <- errors_test_C10 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C10)
# errors_test_mean_C11 <- errors_test_C11 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C11)
# errors_test_mean_C12 <- errors_test_C12 %>% group_by(neighbours_list) %>%
#         summarize_all(mean)
# View(errors_test_mean_C12)

# Mean of all test campimeters errors

errors_test_total <- rbind(errors_test_C1,errors_test_C2,errors_test_C3,
                           errors_test_C4,errors_test_C5,errors_test_C6,
                           errors_test_C7,errors_test_C8,errors_test_C9,
                           errors_test_C10,errors_test_C11,errors_test_C12)
errors_test_mean_total <- errors_test_total %>%
                          summarize_all(mean)
#View(errors_test_mean_total)


# ERRORS histogram (total)
errors_test <- errors_test_total

blues <- brewer.pal(5,"Blues")
breaks_me <- seq(-11,15,0.5)
breaks_mspe <- seq(0,5000,100)
#breaks_mspe_n <- seq(0,8,0.5)
breaks_cor_obs_pred <- seq(-0.1,1,0.05)
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
# Hist(errors_test$var_res, main="Varianza entre errores",ylim=c(0,3000),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_var, col=blues[4], scale="percent")
title("Kriging global: Segmentación 2", line = 0, outer = TRUE, cex.main=1.3)
dev.off()


# RESULTS:

results_test_total <- rbind(results_test_C1,results_test_C2,results_test_C3,
                            results_test_C4,results_test_C5,results_test_C6,
                            results_test_C7,results_test_C8,results_test_C9,
                            results_test_C10,results_test_C11,results_test_C12)
results_test_mean_total <- results_test_total %>%
        summarize_all(mean)

#View(results_test_mean_total)
#View(results_test_total)

dev.off()
blues <- brewer.pal(5,"Blues")
breaks_res <- seq(-300,300,5)

par(mfrow=c(1,1),oma=c(0,0,2,0))
par(mar=c(3,2,2.4,1))
Hist(results_test_total$residuals, main="Residuos totales",ylim=c(0,100000),cex.main=0.7,
     ylab="Porcentaje",xlab=NULL, breaks=breaks_res, col=blues[4], scale="percent")
title("Kriging global: Segmentación 1", line = 0, outer = TRUE, cex.main=1.3)




# Spplot residuals 

eye_plot_all_res <- function(data_all_stats,title,cols,cuts){
        # Select desire statistic
        vars_names <- paste0("residuals_C", 1:n_clust_p) 
        data_all_stats@data <- data_all_stats@data %>% 
                select(all_of(vars_names))
        clust_names = paste0("C",(1:n_clust_p))
        colnames(data_all_stats@data) = clust_names
        # Plot parameters 
        scale.parameter = 1.1
        original.bbox <- eye_clusters@bbox
        xshift = 0  # Shift to right in map units. 
        yshift = 0  # Shift to left in map units.
        edges = original.bbox
        edges[1,] <- (edges[1,] - mean(edges[1,])) * scale.parameter + mean(edges[1,]) + xshift
        edges[2,] <- (edges[2,] - mean(edges[2,])) * scale.parameter + mean(edges[2,]) + yshift
        spplot(data_all_stats, clust_names, cex = 0.7, do.log = F,
               main=list(label=title,cex=1.2),
               col.regions=cols,cuts=cuts,
               #color.key=TRUE,
               sp.layout = list("sp.polygons",eye_clusters,col="dark gray"),
               xlim = edges[1,], ylim = edges[2,],
               key.space="right",
               as.table = TRUE)
}

# Mean of all test campimeters errors by cluster

results_test_mean_C1 <- results_test_C1 %>% group_by(x,y) %>% summarise(residuals_C1 = mean(residuals))
results_test_mean_C2 <- results_test_C2 %>% group_by(x,y) %>% summarise(residuals_C2 = mean(residuals))
results_test_mean_C3 <- results_test_C3 %>% group_by(x,y) %>% summarise(residuals_C3 = mean(residuals))
results_test_mean_C4 <- results_test_C4 %>% group_by(x,y) %>% summarise(residuals_C4 = mean(residuals))
results_test_mean_C5 <- results_test_C5 %>% group_by(x,y) %>% summarise(residuals_C5 = mean(residuals))
results_test_mean_C6 <- results_test_C6 %>% group_by(x,y) %>% summarise(residuals_C6 = mean(residuals))
results_test_mean_C7 <- results_test_C7 %>% group_by(x,y) %>% summarise(residuals_C7 = mean(residuals))
results_test_mean_C8 <- results_test_C8 %>% group_by(x,y) %>% summarise(residuals_C8 = mean(residuals))
results_test_mean_C9 <- results_test_C9 %>% group_by(x,y) %>% summarise(residuals_C9 = mean(residuals))
results_test_mean_C10 <- results_test_C10 %>% group_by(x,y) %>% summarise(residuals_C10 = mean(residuals))
results_test_mean_C11 <- results_test_C11 %>% group_by(x,y) %>% summarise(residuals_C11 = mean(residuals))
results_test_mean_C12 <- results_test_C12 %>% group_by(x,y) %>% summarise(residuals_C12 = mean(residuals))

results_test_residuals_all <- Reduce(function(x, y) merge(x, y, by=c("x","y"), all=TRUE), list(results_test_mean_C1,results_test_mean_C2,results_test_mean_C3,
                                                  results_test_mean_C4,results_test_mean_C5,results_test_mean_C6,
                                                  results_test_mean_C7,results_test_mean_C8,results_test_mean_C9,
                                                  results_test_mean_C10,results_test_mean_C11,results_test_mean_C12))



View(results_test_residuals_all)

cuts <- c(-40,-25,-10,-5,-2.5,-1,-0.5,0.5,1,2.5,5,10,25,40)
cols <- brewer.pal(length(cuts)/2, "OrRd")
cols1 <- rev(cols)
cols2 <- cols
cols <- c(cols1,cols2)
coordinates <- results_test_residuals_all[,1:2]
results_test_residuals_all_sp <- SpatialPointsDataFrame(coordinates, results_test_residuals_all, coords.nrs = numeric(0), 
                                               proj4string = CRS(as.character(NA)), bbox = NULL)

eye_plot_all_res(results_test_residuals_all_sp,"Residuos medios",cols,cuts)


################################################################################################
#                                      9. KRIGGING EXAMPLES                                    #
################################################################################################

# Read global enviroment
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"
setwd(root_path)
load("Global_models_VF8_unique.RData")
n_seg=1

# Eye circle shapefiles
path_eye <- "C:/Users/Agata/MASTER/TFM/Eye_Layers/Eye_cluster"
eye_cir <- readOGR(path_eye, "Eye_Cir") #Eye circular shapefile 
eye_bs <- readOGR(path_eye, "Eye_Blind") #Eye circular + Blind Spot shapefile 
blindSpot <- readOGR(path_eye, "BlindSpot") #Blind Spot shapefile 
grid <- raster(extent(eye_cir))
res(grid) <- 0.5 # Resolution
#proj4string(grid)<-proj4string(eye_cir)
grid <- as(grid, "SpatialPixelsDataFrame") #Pixel Data Frame
eye_pixels <- grid
#plot(eye_pixels)

# Interpolate map
model <- fit.media.cir.C1  #fit.media.sph.C2   #fit.media.cir.C1 
nmax=6
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
                     top=glue("Kriging global por clusters: Segmentación {n_seg}"))
}

