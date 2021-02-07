#####################################################################################################
#                                                                                                   #
#                                               CAMPIMETERS                                         #
#                                                                                                   #
#                                           "GLOBAL" VARIOGRAMS                                     #
#                                                                                                   #
#                                       COMPARISON OF GLOBAL MEANS MODEL                            #
#                                                                                                   #
#                                                                                                   #
#                           1. Validation set errors                                                #
#                             1.1 Mean of error measurements by campimeter and by total             #          
#                             1.2.Display of error measurements                                     #
#                                1.2.1 Evolution of error measures with respect to                  #
#                                      the number of neighbors                                      #
#                                1.2.2 Histograms of error measurements. Fixed n = 8                #  
#                                1.2.3 Table of mean error measurements by model                    #
#                           2. Test set errors                                                      #
#                             2.1 Load fitted models' (train) enviorment                            #
#                             2.2 Compute Test                                                      #
#                             2.3.Display test errors                                               #
#                                                                                                   #
#####################################################################################################

library(lattice)       # Graphics
library(RColorBrewer)  # Color Palettes for plots
library(tools)         # Files manipulation (file_path_sans_ext function)
library(sp)            # Spatial data
library(gstat)         # Geostatistics
library(dplyr)         # Sumarise and group by
library(RcmdrMisc)     # Function Hist()
library(gridExtra)     # Tables as images
library(glue)          # Format str

################################################################################################
#                                  1. VALIDATION SET ERRORS                                    #
################################################################################################

# Read data:
path_errors <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/RESULTS/Global_kriging_500"
setwd(path_errors)

media_lin_results <- read.csv2(paste0(path_errors,"/media_lin_results_val.csv"))
media_lin_errors <- read.csv2(paste0(path_errors,"/media_lin_errors_val.csv"))
media_gauss_results <- read.csv2(paste0(path_errors,"/media_gauss_results_val.csv"))
media_gauss_errors <- read.csv2(paste0(path_errors,"/media_gauss_errors_val.csv"))
media_cir_results <- read.csv2(paste0(path_errors,"/media_cir_results_val.csv"))
media_cir_errors <- read.csv2(paste0(path_errors,"/media_cir_errors_val.csv"))
media_sph_results <- read.csv2(paste0(path_errors,"/media_sph_results_val.csv"))
media_sph_errors <- read.csv2(paste0(path_errors,"/media_sph_errors_val.csv"))
media_exp_results <- read.csv2(paste0(path_errors,"/media_exp_results_val.csv"))
media_exp_errors <- read.csv2(paste0(path_errors,"/media_exp_errors_val.csv"))

# 1.1 Mean of error measurements by campimeter and by total

# 1.1.1 Mean of error measurements by campimetry
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

# 1.1.2 Mean of error measurements by total
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

# 1.2. Display of error measurements

# 1.2.1 Evolution of error measures with respect to the number of neighbors 

dev.off()
Colores <- brewer.pal(10,"Paired") # Paleta Paired
list_col <- c(Colores[1],Colores[2],Colores[3],Colores[7],Colores[5]) #Todos

par(mfcol = c(2, 3))
# Media de la media de los residuos por campimetría
plot(row.names(media_exp_errors_grouped),media_exp_errors_grouped$mean_res,type="l",
      col=Colores[1], xlab= "Número de vecinos", ylab="Error medio",
      ylim=c(-1,1.5), main="Media del error medio", cex.main=0.8)
lines(row.names(media_sph_errors_grouped),media_sph_errors_grouped$mean_res,type="l",
      col=Colores[2], xlab= "Número de vecinos", ylab="Error medio")
lines(row.names(media_gauss_errors_grouped),media_gauss_errors_grouped$mean_res,type="l",
      col=Colores[3], xlab= "Número de vecinos", ylab="Error medio")
lines(row.names(media_lin_errors_grouped),media_lin_errors_grouped$mean_res,type="l",
     col=Colores[7], xlab= "Número de vecinos", ylab="Error medio")
lines(row.names(media_cir_errors_grouped),media_cir_errors_grouped$mean_res,type="l",
      col=Colores[5], xlab= "Número de vecinos", ylab="Error medio")
legend(1.1, -0.55, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
                          "Modelo gausiano de medias", "Modelo lineal de medias",
                          "Modelo circular de medias"),
                          col=list_col, lty=1, cex=0.3, box.lty=0)


# # Media de la media de los residuos normalizados por campimetría
# plot(row.names(media_exp_errors_grouped),media_exp_errors_grouped$mean_zscore,type="l",
#      col=Colores[1], xlab= "Número de vecinos", ylab="Error medio normalizado",
#      ylim=c(-0.038,0.015), main="Media del error medio normalizado", cex.main=0.8)
# lines(row.names(media_sph_errors_grouped),media_sph_errors_grouped$mean_zscore,type="l",
#       col=Colores[2], xlab= "Número de vecinos", ylab="Error medio")
# lines(row.names(media_gauss_errors_grouped),media_gauss_errors_grouped$mean_zscore,type="l",
#       col=Colores[3], xlab= "Número de vecinos", ylab="Error medio")
# lines(row.names(media_lin_errors_grouped),media_lin_errors_grouped$mean_zscore,type="l",
#       col=Colores[7], xlab= "Número de vecinos", ylab="Error medio")
# lines(row.names(media_cir_errors_grouped),media_cir_errors_grouped$mean_zscore,type="l",
#       col=Colores[5], xlab= "Número de vecinos", ylab="Error medio")
# legend(22, -0.018, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
#                           "Modelo gausiano de medias", "Modelo lineal de medias",
#                           "Modelo circular de medias"),
#        col=list_col, lty=1, cex=0.6)


# Media de los errores cuadráticos medios por campimetría
plot(row.names(media_exp_errors_grouped),media_exp_errors_grouped$mspe,type="l",
     col=Colores[1], xlab= "Número de vecinos", ylab="Error cuadrático medio",
     ylim=c(550,1300), main="Media del error cuadrático medio", cex.main=0.8)
lines(row.names(media_sph_errors_grouped),media_sph_errors_grouped$mspe,type="l",
      col=Colores[2], xlab= "Número de vecinos", ylab="Error cuadrático medio")
lines(row.names(media_gauss_errors_grouped),media_gauss_errors_grouped$mspe,type="l",
      col=Colores[3], xlab= "Número de vecinos", ylab="Error cuadrático medio")
lines(row.names(media_lin_errors_grouped),media_lin_errors_grouped$mspe,type="l",
      col=Colores[7], xlab= "Número de vecinos", ylab="Error cuadrático medio")
lines(row.names(media_cir_errors_grouped),media_cir_errors_grouped$mspe,type="l",
      col=Colores[5], xlab= "Número de vecinos", ylab="Error cuadrático medio")
legend(7.5, 1290, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
                            "Modelo gausiano de medias", "Modelo lineal de medias",
                            "Modelo circular de medias"),
       col=list_col, lty=1, cex=0.3, box.lty=0)
       
# Media de los errores cuadráticos normalizados medios por campimetría
plot(row.names(media_exp_errors_grouped),media_exp_errors_grouped$mspe_norm,type="l",
     col=Colores[1], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado",
     ylim=c(0.4,1), main="Media del error cuadrático medio normalizado", cex.main=0.8)
lines(row.names(media_sph_errors_grouped),media_sph_errors_grouped$mspe_norm,type="l",
      col=Colores[2], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado")
lines(row.names(media_gauss_errors_grouped),media_gauss_errors_grouped$mspe_norm,type="l",
      col=Colores[3], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado")
lines(row.names(media_lin_errors_grouped),media_lin_errors_grouped$mspe_norm,type="l",
      col=Colores[7], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado")
lines(row.names(media_cir_errors_grouped),media_cir_errors_grouped$mspe_norm,type="l",
      col=Colores[5], xlab= "Número de vecinos", ylab="Error cuadrático medio normalizado")
legend(7.5, 0.99, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
                          "Modelo gausiano de medias", "Modelo lineal de medias",
                          "Modelo circular de medias"),
       col=list_col, lty=1, cex=0.3, box.lty=0)

# Media de la correlación entre valores observados y predichos
plot(row.names(media_exp_errors_grouped),media_exp_errors_grouped$cor_obs_pred,type="l",
     col=Colores[1], xlab= "Número de vecinos", ylab="Correlación observado-predicho",
     ylim=c(0.59,0.8), main="Media de la correlación observado-predicho", cex.main=0.8)
lines(row.names(media_sph_errors_grouped),media_sph_errors_grouped$cor_obs_pred,type="l",
      col=Colores[2], xlab= "Número de vecinos", ylab="Correlación observado-predicho")
lines(row.names(media_gauss_errors_grouped),media_gauss_errors_grouped$cor_obs_pred,type="l",
      col=Colores[3], xlab= "Número de vecinos", ylab="Correlación observado-predicho")
lines(row.names(media_lin_errors_grouped),media_lin_errors_grouped$cor_obs_pred,type="l",
      col=Colores[7], xlab= "Número de vecinos", ylab="Correlación observado-predicho")
lines(row.names(media_cir_errors_grouped),media_cir_errors_grouped$cor_obs_pred,type="l",
      col=Colores[5], xlab= "Número de vecinos", ylab="Correlación observado-predicho")
legend(7.5, 0.63, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
                         "Modelo gausiano de medias", "Modelo lineal de medias",
                         "Modelo circular de medias"),
       col=list_col, lty=1, cex=0.3, box.lty=0)

# Media de la correlación entre valores predichos y errores normalizados
plot(row.names(media_exp_errors_grouped),media_exp_errors_grouped$cor_pred_res_n,type="l",
     col=Colores[1], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho",
     ylim=c(-0.45,0.2), main="Media de la correlación error normalizado-predicho", cex.main=0.8)
lines(row.names(media_sph_errors_grouped),media_sph_errors_grouped$cor_pred_res_n,type="l",
      col=Colores[2], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho")
lines(row.names(media_gauss_errors_grouped),media_gauss_errors_grouped$cor_pred_res_n,type="l",
      col=Colores[3], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho")
lines(row.names(media_lin_errors_grouped),media_lin_errors_grouped$cor_pred_res_n,type="l",
      col=Colores[7], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho")
lines(row.names(media_cir_errors_grouped),media_cir_errors_grouped$cor_pred_res_n,type="l",
      col=Colores[5], xlab= "Número de vecinos", ylab="Correlación error normalizado-predicho")
legend(7.5, -0.34, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
                         "Modelo gausiano de medias", "Modelo lineal de medias",
                         "Modelo circular de medias"),
       col=list_col, lty=1, cex=0.3, box.lty=0)

# Proportion of normalized errors above 2.5
plot(row.names(media_exp_errors_grouped),media_exp_errors_grouped$percent,type="l",
     col=Colores[1], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5",
     ylim=c(0.4,1.2), main="Media del porcentaje errores normalizados > 2.5", cex.main=0.8)
lines(row.names(media_sph_errors_grouped),media_sph_errors_grouped$percent,type="l",
      col=Colores[2], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5")
lines(row.names(media_gauss_errors_grouped),media_gauss_errors_grouped$percent,type="l",
      col=Colores[3], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5")
lines(row.names(media_lin_errors_grouped),media_lin_errors_grouped$percent,type="l",
      col=Colores[7], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5")
lines(row.names(media_cir_errors_grouped),media_cir_errors_grouped$percent,type="l",
      col=Colores[5], xlab= "Número de vecinos", ylab="Porcentaje errores normalizados > 2.5")
legend(7.5, 0.54, legend=c("Modelo exponencial de medias", "Modelo esférico de medias",
                          "Modelo gausiano de medias", "Modelo lineal de medias",
                          "Modelo circular de medias"),
       col=list_col, lty=1, cex=0.3, box.lty=0)

mtext("Medidas de error por número de vecinos: Kriging Global", side = 3, line = -1.5, outer = TRUE, cex=0.8)

# 1.2.2 Histograms of error measurements. Fixed n = 6

# 1.2.2.1. Preparar datos
# -Media linear
media_lin_res_n <- media_lin_results[which(media_lin_results$n_neighbours==6),]
media_lin_err_n <- media_lin_errors[which(media_lin_errors$neighbours_list==6),]
# -Media gaussiano
media_gauss_res_n <- media_gauss_results[which(media_gauss_results$n_neighbours==6),]
media_gauss_err_n <- media_gauss_errors[which(media_gauss_errors$neighbours_list==6),]
# -Media circular
media_cir_res_n <- media_cir_results[which(media_cir_results$n_neighbours==6),]
media_cir_err_n <- media_cir_errors[which(media_cir_errors$neighbours_list==6),]
# -Media esférico
media_sph_res_n <- media_sph_results[which(media_sph_results$n_neighbours==6),]
media_sph_err_n <- media_sph_errors[which(media_sph_errors$neighbours_list==6),]
# -Media exponencial
media_exp_res_n <- media_exp_results[which(media_exp_results$n_neighbours==6),]
media_exp_err_n <- media_exp_errors[which(media_exp_errors$neighbours_list==6),]

# # HISTOGRAMS
# blues <- brewer.pal(5,"Blues")
# color <- blues[4]
# breaks_me <- seq(-5,10,1.25)
# breaks_mspe <- seq(0,2600,200)
# breaks_mspe_n <- seq(0,8,0.5)
# breaks_cor_obs_pred <- seq(0,1,0.1)
# breaks_cor_pred_zscore <- seq(-0.5,1,0.1)
# breaks_percent <- seq(0,15,0.5)
# View(media_lin_err_8)
# 
# # - Media linear
# par(mfrow=c(2,3),oma=c(0,0,2,0))
# par(mar=c(3,2,2.4,1))
# Hist(media_lin_err_8$mean_res, main="Errores medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_me, col=color, scale="percent")
# Hist(media_lin_err_8$mspe, main="Errores cuadráticos medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe, col=color, scale="percent")
# Hist(media_lin_err_8$mspe_norm, main="Errores cuadráticos medios normalizados",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe_n, col=color, scale="percent")
# Hist(media_lin_err_8$cor_obs_pred, main="Correlaciones observado-predicho",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_obs_pred, col=color, scale="percent")
# Hist(media_lin_err_8$cor_pred_resN, main="Correlaciones predicho-error normalizado",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_pred_zscore, col=color, scale="percent")
# Hist(media_lin_err_8$zscores_percent, main="Porcentaje errores normalizados > 2.5",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_percent, col=color, scale="percent")
# title("Modelo lineal de medias n=8", line = 0.5, outer = TRUE, cex.main=1.3)
# # - Media gaussiano
# par(mfrow=c(2,3),oma=c(0,0,2,0))
# par(mar=c(3,2,2.4,1))
# Hist(media_gauss_err_8$mean_res, main="Errores medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_me, col=color, scale="percent")
# Hist(media_gauss_err_8$mspe, main="Errores cuadráticos medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe, col=color, scale="percent")
# Hist(media_gauss_err_8$mspe_norm, main="Errores cuadráticos medios normalizados",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe_n, col=color, scale="percent")
# Hist(media_gauss_err_8$cor_obs_pred, main="Correlaciones observado-predicho",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_obs_pred, col=color, scale="percent")
# Hist(media_gauss_err_8$cor_pred_resN, main="Correlaciones predicho-error normalizado",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_pred_zscore, col=color, scale="percent")
# Hist(media_gauss_err_8$zscores_percent, main="Porcentaje errores normalizados > 2.5",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_percent, col=color, scale="percent")
# title("Modelo gaussiano de medias n=8", line = 0.5, outer = TRUE, cex.main=1.3)
# # - Media circular
# par(mfrow=c(2,3),oma=c(0,0,2,0))
# par(mar=c(3,2,2.4,1))
# Hist(media_cir_err_8$mean_res, main="Errores medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_me, col=color, scale="percent")
# Hist(media_cir_err_8$mspe, main="Errores cuadráticos medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe, col=color, scale="percent")
# Hist(media_cir_err_8$mspe_norm, main="Errores cuadráticos medios normalizados",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe_n, col=color, scale="percent")
# Hist(media_cir_err_8$cor_obs_pred, main="Correlaciones observado-predicho",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_obs_pred, col=color, scale="percent")
# Hist(media_cir_err_8$cor_pred_resN, main="Correlaciones predicho-error normalizado",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_pred_zscore, col=color, scale="percent")
# Hist(media_cir_err_8$zscores_percent, main="Porcentaje errores normalizados > 2.5",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_percent, col=color, scale="percent")
# title("Modelo circular de medias n=8", line = 0.5, outer = TRUE, cex.main=1.3)
# # - Media esférico
# par(mfrow=c(2,3),oma=c(0,0,2,0))
# par(mar=c(3,2,2.4,1))
# Hist(media_sph_err_8$mean_res, main="Errores medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_me, col=color, scale="percent")
# Hist(media_sph_err_8$mspe, main="Errores cuadráticos medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe, col=color, scale="percent")
# Hist(media_sph_err_8$mspe_norm, main="Errores cuadráticos medios normalizados",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe_n, col=color, scale="percent")
# Hist(media_sph_err_8$cor_obs_pred, main="Correlaciones observado-predicho",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_obs_pred, col=color, scale="percent")
# Hist(media_sph_err_8$cor_pred_resN, main="Correlaciones predicho-error normalizado",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_pred_zscore, col=color, scale="percent")
# Hist(media_sph_err_8$zscores_percent, main="Porcentaje errores normalizados > 2.5",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_percent, col=color, scale="percent")
# title("Modelo esférico de medias n=8", line = 0.5, outer = TRUE, cex.main=1.3)
# # - Media exponencial
# par(mfrow=c(2,3),oma=c(0,0,2,0))
# par(mar=c(3,2,2.4,1))
# Hist(media_exp_err_8$mean_res, main="Errores medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_me, col=color, scale="percent")
# Hist(media_exp_err_8$mspe, main="Errores cuadráticos medios",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe, col=color, scale="percent")
# Hist(media_exp_err_8$mspe_norm, main="Errores cuadráticos medios normalizados",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_mspe_n, col=color, scale="percent")
# Hist(media_exp_err_8$cor_obs_pred, main="Correlaciones observado-predicho",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_obs_pred, col=color, scale="percent")
# Hist(media_exp_err_8$cor_pred_resN, main="Correlaciones predicho-error normalizado",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_cor_pred_zscore, col=color, scale="percent")
# Hist(media_exp_err_8$zscores_percent, main="Porcentaje errores normalizados > 2.5",ylim=c(0,520),cex.main=0.7,
#      ylab="Porcentaje",xlab=NULL, breaks=breaks_percent, col=color, scale="percent")
# title("Modelo exponencial de medias n=8", line = 0.5, outer = TRUE, cex.main=1.3)

# 1.2.3 Table of mean error measurements by model. Fixed n = 6

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

grid.table(errores_modelos)



################################################################################################
#                                  3. TEST GLOBAL MODEL                                        #
################################################################################################

# Function: Kigriging errors with variogram fitted to model residuals

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

# 2.1 Load fitted models' (train) enviorment 

# Read global enviroment (Fitted Global models)
path_models <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Global_models.RData"
load(path_models)

# Select fitted circular model 
model_cir <- fit.media.cir

# 2.2 Compute Test -> Sample 4

# Read sample 4
# root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/DATA/Data_IDW"
# setwd(root_path)
# ns <- 4
# path_sample <- file.path(paste0(root_path,"/",glue("5260_Patients_Sample{ns}")))
# file_names <- file_path_sans_ext(list.files(path_sample, pattern = "*.csv"))
# file_names_test <- file_names
# media_cir_global <- krige.results(file_names_test,model_cir,8)
# write.csv2(media_cir_global[[1]], file.path(getwd(),paste0("media_cir_global_results_test.csv")))
# write.csv2(media_cir_global[[2]], file.path(getwd(),paste0("media_cir_global_errors_test.csv")))

# 2.3. Display test errors

root_path <- "C:/Users/Agata/MASTER/TFM/Data_Scripts/RESULTS/Global_kriging"
setwd(root_path)
errors_test <- read.csv2(file.path(root_path,"media_exp_global_errors_test.csv"), sep = ";")
results_test <-read.csv2(file.path(root_path,"media_exp_global_results_test.csv"), sep = ";")

# Mean of all test campimeters errors 

errors_test_mean <- errors_test %>% group_by(neighbours_list) %>%
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
title("Modelo Kriging Global", line = 0.5, outer = TRUE, cex.main=1.3)
dev.off()


# RESULTS

results_test_total <- results_test
dev.off()
blues <- brewer.pal(5,"Blues")
breaks_res <- seq(-300,300,5)

par(mfrow=c(1,1),oma=c(0,0,2,0))
par(mar=c(3,2,2.4,1))
Hist(results_test_total$residuals, main="Residuos totales",ylim=c(0,100000),cex.main=0.7,
     ylab="Porcentaje",xlab=NULL, breaks=breaks_res, col=blues[4], scale="percent")
title("Kriging global", line = 0, outer = TRUE, cex.main=1.3)




# Spplot mean residuals

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
