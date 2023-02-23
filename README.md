
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rEAD

<!-- badges: start -->
<!-- badges: end -->

rEAD is a R package for import and analyse GC-EAD data.

## Installation

You can install the development version of rEAD from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JHuguenin/rEAD")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(rEAD)

## gestion de repertoir de travail # 
wd = "C:/Users/huguenin/Documents/R/rEAD_moustiques" # working directory
setwd(wd)  # changement du repertoire dans R

## 1er import #
a04 <- import.GC.EAD(file_csv = "data/2023-01-25_hela_aedes04_Speed_025_split4_20Hz_EAD.csv")
# check la figure dynamique + la figure de calibration dans ~/figures
gcead.graph(a04)

## calcul du decalage entre les deux FID #
recalc.delay(pk_EAD = 8.09559, pk_FID = 8.1579,pk_mat = a04) # 374
recalc.delay(pk_EAD = 9.64715, pk_FID = 9.703967,pk_mat = a04) # 341
recalc.delay(pk_EAD = 10.2951, pk_FID = 10.3486,pk_mat = a04) # 321
recalc.delay(pk_EAD = 10.641198, pk_FID = 10.697349,pk_mat = a04) # 337

## importation de tous les echantillons #
a04 <- import.GC.EAD(file_csv = "data/2023-01-25_hela_aedes04_Speed_025_split4_20Hz_EAD.csv",delay = 340)
a05 <- import.GC.EAD(file_csv = "data/2023-01-25_hela_aedes05_Speed_025_split4_20Hz_EAD.csv",delay = 340)
a07 <- import.GC.EAD(file_csv = "data/2023-01-26_hela_aedes07_Speed_025_split4_20Hz_EAD.csv",delay = 340)
a08 <- import.GC.EAD(file_csv = "data/2023-01-26_hela_aedes08_Speed_025_split4_20Hz_EAD.csv",delay = 340)
a09 <- import.GC.EAD(file_csv = "data/2023-01-26_hela_aedes09_Speed_025_split4_20Hz_EAD.csv",delay = 340)
a10 <- import.GC.EAD(file_csv = "data/2023-01-26_hela_aedes10_Speed_025_split4_20Hz_EAD.csv",delay = 340)
a11 <- import.GC.EAD(file_csv = "data/2023-01-26_hela_aedes11_Speed_025_split4_20Hz_EAD.csv",delay = 340)
a12 <- import.GC.EAD(file_csv = "data/2023-01-26_hela_aedes12_Speed_025_split4_20Hz_EAD.csv",delay = 340)
a14 <- import.GC.EAD(file_csv = "data/2023-01-27_hela_aedes14_Speed_025_split4_20Hz_EAD.csv",delay = 340)
a15 <- import.GC.EAD(file_csv = "data/2023-01-27_hela_aedes15_Speed_025_split4_20Hz_EAD.csv",delay = 340)
# regarder les figures comparatives de calibration dans ~/figures

## rassembler les echantillons #
list_test <- gcead.merge(a04, a05, a07, a08, a09, a10, a11, a12, a14, a15, # tous les echantillons importes
                         gcead_names = c("a04","a05","a07","a08","a09","a10","a11","a12","a14","a15"), # leurs noms reduits
                         RT_limits = c(6,18), # le temps tronques des parties inutiles
                         gap_FID = 1) # un facteur multiplicatif pour le signal FID

list_test <- m.EAD.norm(shift = -0.5, amplitude = 5, overlay = FALSE, pk_mat = list_test) # une mise en forme des signaux

## modification du pic 11.2268 #
list_test <- FID.modify.peak(pk = 11.2268, bm_sup = 11.23785, pk_mat = list_test)
list_test <- FID.create.peak(11.2379, 11.269,pk_mat = list_test)

## suppression du pic 9.74254 #
list_test <- FID.delete.peak(9.74,pk_mat = list_test)

## baseline et calcule des statistiques #
list_test <- m.EAD.baseline(cr_shift = -0.05, double_baseline = TRUE, pk_mat = list_test) # baseline

list_test <- m.EAD.average(list_test,"mean", shift = -0.5, overlay = TRUE) # calcul de la moyenne
#list_test <- m.EAD.average(list_test,"median", shift = -0.5, overlay = TRUE) # calcule de la mediane

## graphes #
mead.graph(list_test, view_raw =  TRUE) # differents graphes
mead.graph(list_test, view_raw =  FALSE)
mead.graph(list_test, view_raw =  FALSE, view_selection = "mean")
mead.graph(list_test, view_raw =  FALSE, view_selection = c("mean","sd"))
mead.graph(list_test, view_raw =  FALSE, view_selection = c("mean","sd"), view_snr = TRUE)
mead.graph(list_test, view_raw =  FALSE, view_snr = TRUE)

## measure #
depol_param <- pk.measure(list_test)  # calcul des intensites et de la duree de depolarisation des pics EAD
# View(depol_param[[1]])
```
