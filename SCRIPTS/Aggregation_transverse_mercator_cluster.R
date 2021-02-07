#####################################################################################################
#                                                                                                   #
#                                               CAMPIMETERS                                         #
#                                                                                                   #
#                             Aggregation and transformation of original data                       #
#                                                                                                   #
#                         1. Folder contains campimeters for each patient                           #
#                         2. Folder contains a random sample (1 campimeter for each patient)        #
#                         3. Dataframe that contains in each columns a different                    #
#                            1-campimeter-patient sample                                            #
#                         4. Folder contains all campimeters                                        #
#                                                                                                   #
#####################################################################################################

library(lattice)    # Graphics
library(sp)         # Spatial data
library(fs)         # File manipulations
library(tidyverse)  # Data manipulation
library(tools)      # Files manipulations
library(stringr)    # Character manipulation
library(geosphere)  # Usage of mercator projection

#####################################################################################################

# 1. Folder contains campimeters for each patient 
# For each patient, create a folder with all their campimetries in tranversal mercator projection
# Reclassification of measure points into internal clusters. 

path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/patient_by_campimetries" # Original data
setwd(path)   # Fixed working directory
file_names <- list.files(getwd(), pattern = "*.csv")

# Function  10 cluster visual fields assignation (Only for left eye campimeters) 
cluster_assign10 <- function(campimeter){
  c0 <- c(43)
  c1 <- c(7,8,21,52)
  c2 <- c(9,10,24,55)
  c3 <- c(22,53,56,57,41,5,27)
  c4 <- c(23,54,3,58,59,4,42,6,28)
  c5 <- c(1,2,37,38,30,49)
  c6 <- c(50,31,39,40)
  c7 <- c(34,33,18,45,12,11,44,29)
  c8 <- c(35,36,19,46,13,14,47,32)
  c9 <- c(15,25,48,17)
  c10 <- c(16,26,51,20)  
  campimeter$label <- as.character(row.names(campimeter))
  campimeter$cluster = NA
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c0, NA, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c1, 1, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c2, 2, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c3, 3, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c4, 4, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c5, 5, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c6, 6, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c7, 7, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c8, 8, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c9, 9, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c10, 10, cluster))
  return(campimeter)
}

# Function  5 cluster visual fields assignation (Only for left eye campimeters) 
cluster_assign5 <- function(campimeter){
  c1 <- c(34,18,45,12,49,30,27,5)
  c2 <- c(33,11,44,17,29,48,25,15)
  c3 <- c(36,14,47,20,32,51,26,16)
  c4 <- c(35,19,46,13,50,31,28,6)
  c5 <- c(43,7,8,9,10,21,22,23,24,52,53,54,55,56,57,58,59,1,2,3,4,41,42,37,38,39,40)
  campimeter$label <- as.character(row.names(campimeter))
  campimeter$cluster = NA
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c1, 1, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c2, 2, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c3, 3, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c4, 4, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c5, 5, cluster))
  return(campimeter)
}

# Function 8 cluster visual fields assignation (Only for left eye campimeters) 
cluster_assign8 <- function(campimeter){
  c0 <- c(43)
  c1 <- c(34,18,45,12,49,30,27,5)
  c2 <- c(33,11,44,17,29,48,25,15)
  c3 <- c(36,14,47,20,32,51,26,16)
  c4 <- c(35,19,46,13,50,31,28,6)
  c5 <- c(8,22,57,53,2,41,38)
  c6 <- c(7,21,56,52,1,37)
  c7 <- c(10,24,59,55,4,40)
  c8 <- c(9,23,58,54,3,42,39)
  campimeter$label <- as.character(row.names(campimeter))
  campimeter$cluster = NA
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c0, NA, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c1, 1, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c2, 2, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c3, 3, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c4, 4, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c5, 5, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c6, 6, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c7, 7, cluster))
  campimeter <- mutate(campimeter, cluster = ifelse(label %in% c8, 8, cluster))
  return(campimeter)
}

for(i in file_names){
  filepath <- file.path(file.path(getwd(),paste(i,sep="")))
  aux <- assign(i, read.delim(filepath,sep = ";"))
  data_info <- aux[c(2,3,4,26,27,28)]
  data_info$ncamp <- seq.int(nrow(aux)) # Añadir número campimetría
  colnames(data_info) <- c("id","birthdate","sex","duration","date","time","ncamp")
  id <- data_info[1,1]
  dir.create(format(id)) # New folder Id patient name
  write.csv2(data_info, file.path(format(id),paste0((format(id)),"_info.csv")), row.names = FALSE)
  data_camp <- aux[,44:ncol(aux)]
  for (j in 1:nrow(data_camp)) {
    v_camp <- as.numeric(data_camp[j,]) # jth campimetry
    l_camp <- split(v_camp, ceiling(seq_along(v_camp)/5)) # Data list by coordinates
    camp <- data.frame(matrix(unlist(l_camp), nrow=length(l_camp), byrow=T)) #Df ith campimetry
    colnames(camp) <- c("X","Y","R1","R2","Ref")
    # 'Standarize' values (Real value - Reference Value) -> Visual Defect
    camp$VNorm <- (camp$Ref-camp$R1)  # New column Value Reference Value - Real value
    # -Assign points into 10 cluster visual fields
    #camp <- cluster_assign10(camp)
    # -Assign points into 5 cluster visual fields
    #camp <- cluster_assign5(camp)
    # -Assign points into 8 cluster visual fields
    #camp <- cluster_assign8(camp)
    # Projected mercator coordinates (r = 1 cm) -> Coordinates in cm
    camp_mercator <- mercator(as.data.frame(c((camp["X"]/10),(camp["Y"]/10))), inverse=FALSE, r=1)
    camp <- cbind(camp,camp_mercator)
    #camp <- data.frame(camp$x,camp$y,camp$R1, camp$R2, camp$Ref, camp$VNorm, camp$cluster) #Cluster
    camp <- data.frame(camp$x,camp$y,camp$R1, camp$R2, camp$Ref, camp$VNorm)                #No Cluster
    #colnames(camp) <- c("x","y","R1","R2","Ref","VNorm","Cluster")                         #Cluster
    colnames(camp) <- c("x","y","R1","R2","Ref","VNorm")                                    #No Cluster        
    #  cm to mm  #
    camp$x <- camp$x * 10
    camp$y <- camp$y * 10
    write.csv2(camp, file.path(format(id),paste0(format(id),"_camp", format(j), ".csv")), row.names = FALSE)
  }
}

#####################################################################################################

# 2. Folder contains a random sample (1 or more campimeter for each patient)
# Extract (3) random campimeter of each patient and put all the campimetries together 
# in the same folder -> 5261_Patients_Sample_1, 5261_Patients_Sample_2 ...

set.seed(123)
#Origin folder 
model_name <- "Data_IDW"
root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA"   
path_model <- file.path(paste0(root_path,"/",model_name))
path_patients <- file.path(paste0(root_path,"/",model_name,"/","5260_Patients")) 
#Destiny folders
to_path1 <- file.path(paste0(path_model,"/","5260_Patients_Sample1"))  
to_path2 <- file.path(paste0(path_model,"/","5260_Patients_Sample2"))  
to_path3 <- file.path(paste0(path_model,"/","5260_Patients_Sample3")) 
to_path4 <- file.path(paste0(path_model,"/","5260_Patients_Sample4"))

folder_names <- list.files(path_patients)
for (i in folder_names){
  files_names <- list.files(file.path(path_patients,paste(i,sep="")), pattern = "*.csv")
  files_names <- files_names[!grepl('info', files_names)] # Remove "info" file 
  random_campimeter <- sample(files_names,4) # Select (4) random campimeter of a patient
  random_campimeter1 <- random_campimeter[1]
  random_campimeter2 <- random_campimeter[2]
  random_campimeter3 <- random_campimeter[3]
  random_campimeter4 <- random_campimeter[4]
   
  # Copy random campimeter to the sample folder
  file.copy(file.path(path_patients,paste(i,sep=""),paste(random_campimeter1,sep="")), to_path1)
  file.copy(file.path(path_patients,paste(i,sep=""),paste(random_campimeter2,sep="")), to_path2)
  file.copy(file.path(path_patients,paste(i,sep=""),paste(random_campimeter3,sep="")), to_path3)
  file.copy(file.path(path_patients,paste(i,sep=""),paste(random_campimeter4,sep="")), to_path4)
}


#####################################################################################################
# 3. Folder contains all campimeters  
# Put all the campimetries together in the same folder

path <- "C:/Users/Agata/MASTER/TFM/Scripts/ts"  # Origin folder (Last step)
setwd(path)   # Fixed working directory
to_path <- "C:/Users/Agata/MASTER/TFM/Scripts/All_campimeters_c"  # Destiny Folder
folder_names <- list.files(getwd())
for (i in folder_names){
  files_names <- list.files(file.path(getwd(),paste(i,sep="")), pattern = "*.csv")
  for (j in files_names){
    file.copy(file.path(getwd(),paste(i,sep=""),paste(j,sep="")), to_path, overwrite = TRUE )
  }
}
setwd(to_path)
file_names  <- list.files(getwd(), pattern = "*.csv")
# Delete *csv contains info patient
for (i in file_names){
  if (str_detect(format(i),"info")){     # If file name contains "info" -> Delete
    file.remove(i)
  }
}

#####################################################################################################

# 5. Read all campimetries for a patient folder and assing into a variable

#file_names <- file_path_sans_ext(list.files(file.path(getwd(),format(id)), pattern = "*.csv"))
#for(i in file_names){
#  filepath <- file.path(file.path(getwd(),format(id)),paste(i,sep=""))
#  assign(i, read.csv2(filepath,sep = ";"))
#}