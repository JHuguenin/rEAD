library(baseline)
library(dygraphs)
library(magrittr)
library(MALDIquant)
library(pracma)
library(RColorBrewer)
library(readr)
library(rmarkdown)
library(stringr)
library(xts)

library(rEAD)
#########################
######  TEST  ###########
#########################

##### MERGE ####
wd = "C:/Users/huguenin/Documents/R/GCEAD" # working directory (repertoire de travail)
setwd(wd)  # changement du repertoire dans R

b24 <- import.GC.EAD("data/Abeilles/2022-07-28_lavandeO3pool_bee24_Speed_025_spitless_20Hz_EAD.csv",delay = 438)
b26 <- import.GC.EAD("data/Abeilles/2022-07-28_lavandeO3pool_bee26_Speed_025_spitless_20Hz_EAD.csv",delay = 438)
b29 <- import.GC.EAD("data/Abeilles/2022-07-29_lavandeO3pool_bee29_Speed_025_spitless_20Hz_EAD.csv",delay = 438)
b30 <- import.GC.EAD("data/Abeilles/2022-07-29_lavandeO3pool_bee30_Speed_025_spitless_20Hz_EAD.csv",delay = 438)
b31 <- import.GC.EAD("data/Abeilles/2022-07-29_lavandeO3pool_bee31_Speed_025_spitless_20Hz_EAD.csv",delay = 438)

list_test <- gcead.merge(b24, b26, b29, b30, b31,
                         gcead_names = c("b24","b26","b29","b30","b31"),
                         RT_limits = c(4.5,16))
# list_save <- list_test
# list_test <- list_save
list_test <- m.EAD.norm(shift = -0.5, amplitude = c(1,-1,1,1,1), overlay = FALSE, pk_mat = list_test)
list_test <- m.EAD.norm(shift = -0.5, amplitude = 2, overlay = TRUE, pk_mat = list_test)
list_test <- m.EAD.baseline(cr_shift = -0.05, pk_mat = list_test)
list_test <- m.EAD.average(list_test,"mean",)
list_test <- m.EAD.average(list_test,"median")
mead.graph(list_test, save_htlm = TRUE)
mead.graph(list_test, view_all =  FALSE, save_htlm = TRUE)

list_test <- FID.delete.window(X2 = 7.1, pk_mat = list_test)
list_test <- FID.delete.window(X1 = 15.5, pk_mat = list_test)

list_test <- FID.create.peak(11.12,11.21,pk_mat = list_test)
list_test <- FID.create.peak(12.35,12.43,pk_mat = list_test)
# great =)

list_test <- FID.modify.peak(pk = 11.63, bm_sup = 11.686, pk_mat = list_test)
list_test <- FID.modify.peak(pk = 8.13, bm_sup = 8.16, pk_mat = list_test)
list_test <- FID.modify.peak(pk = 9.28, bm_inf = 9.24, bm_sup = 9.296, pk_mat = list_test)
list_test <- m.EAD.baseline(cr_shift = -0.05, pk_mat = list_test)
list_test <- m.EAD.average(list_test,"median", shift = -0.5, overlay = TRUE)
list_test <- m.EAD.average(list_test,"mean", shift = -0.5, overlay = TRUE)
mead.graph(list_test,view_all = TRUE)
mead.graph(list_test,view_all = FALSE)

#              !!
#             !!!!
#            !!  !!
#           !!    !!
#          !!  !!  !!
#         !!   !!   !!
#        !!    !!    !!
#       !!            !!
#      !!      !!      !!
#     !!                !!
#    !!!!!!!!!!!!!!!!!!!!!!

pk_mat <- list_test
