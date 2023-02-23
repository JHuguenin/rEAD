# FUNCTIONS #

#' rEAD : read EAD data
#'
#' import, pretraitement and tools for analyze data of GC-EAD
#'
#' @docType package
#' @name rEAD
#'
#' @import dygraphs
#' @import graphics
#' @import grDevices
#' @import methods
#' @import stringr
#' @import utils
#' @importFrom baseline baseline
#' @importFrom htmlwidgets saveWidget
#' @importFrom magrittr add
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr subtract
#' @importFrom MALDIquant createMassSpectrum
#' @importFrom MALDIquant detectPeaks
#' @importFrom plotly add_trace
#' @importFrom plotly plot_ly
#' @importFrom pracma savgol
#' @importFrom RColorBrewer brewer.pal
#' @importFrom readr read_delim
#' @importFrom stats lm
#' @importFrom stats median
NULL

# 00 Library ####

# library(dygraphs)
# library(grDevices)
# library(utils)
# library(baseline)
# library(magrittr)
# library(MALDIquant)
# library(plotly)
# library(pracma)
# library(RColorBrewer)
# library(readr)
# library(stats)
# library(stringr)

# 01 Import ####

#' Import data
#'
#' La fonction importe le file_csv correspondant a la sorti csv de l'instrument
#' de syntech. Le signal EAD est en colonne 3 et le signal FID en colonne 4. Si
#' delay est precise, le signal EAD est decale d'autant de pixel.
#'
#' Le signal FID du syntech est bien souvent sature pour les pics les plus intenses.
#' Cela oblige a mettre en place plusieurs etapes pour corriger le probleme.
#' Premierement, il faut rehausse le signal pour que la ligne de base (la mediane)
#' soit superieur a zero. Puis il faut ramener les parties tronquees a l'intensite maximum.
#'
#' Ceci fait, la partie FID est chargee. Il y a quatre colonnes le temps en ms,
#' en min, le RI et le FID. Si skip_time est precise, le debut du FID est tronque
#' dans cette partie.
#'
#' Suite a ça, les cinq pics les plus intenses qui ne sont pas satures sont detectes
#' sur le signal FID de Syntech (GC-EAD). Les cinq meme pics sont detectes sur le
#' signal FID d'agilent. Cela permet de passer a la phase de d'etalonnage afin de
#' synchroniser l'echelle de temps avec le signal EAD de syntech.
#'
#' Enfin, le signal du FID d'agilent est recalcule pour avoir le meme nombre de
#' point que le signal EAD de syntech.
#' Dans une derniere partie, les pics sont recalcule et les donnees misent en forme
#' pour obtenir un seul objet de classe gcead.
#'
#'
#' @param file_csv name of the file to be imported
#' @param num the spectrum number. Only if several spectra have been concatenated to a single file
#' @param wd working directory
#' @param delay delays between FID and EAD spectrum
#' @param skip_time delete x firsts minutes (X = number of minutes)
#' @param FID_calibration calibration file
#' @param param parameters of EAD HalfWindowsSize, EAD SNR, FID HalfWindowsSize and FID SNR
#' @param skip_pk_cal_GCEAD numeric. Skip peak, higher in first, during the calibration phase.
#' @param skip_pk_cal_FID numeric. Skip peak, higher in first, during the calibration phase.
#' @param save_fid_cal logical. Save the plot html of calibration transfert from FID to syntech format
#' @param FID_cal numeric. Value depend of window width filter for calibration transfert
#'
#' @return a gcead object with spectra and metadata
#' @export
#'
#' @examples
#' # phase 1
#' # gcead <- import.GC.EAD(file_csv = "data/sbe18.csv", num = 1, delay = 0, FID_calibration = "data/sbe21_FID2.csv")
#' # gcead <- print.GCEAD()
import.GC.EAD <- function(file_csv = "name.csv", num = 1, wd = NULL, delay = 0,
                          skip_time = 6, FID_calibration = NULL, param = c(40,7,20,7),
                          skip_pk_cal_GCEAD = NA, skip_pk_cal_FID = NA,
                          save_fid_cal = TRUE, FID_cal = -10^6){
  # check ####
  if(is.null(wd) == TRUE) wd <- getwd()
  if (!is.character(wd)) stop("'wd' must be character")
  if(("figures" %in% dir(wd))==FALSE) dir.create(paste0(wd,"/figures"))
  if (!is.character(file_csv)) stop("'file_csv' must be character")
  if (!is.numeric(num)) stop("'num' must be numeric")
  if (!is.numeric(delay)) stop("'delay' must be numeric")
  if (!is.numeric(skip_time)) stop("'skip_time' must be numeric")
  if (!is.numeric(param)) stop("'param' must be numeric")
  if (length(param) != 4) stop("'param' must have for values (EAD HalfWindowsSize, EAD SNR, FID HalfWindowsSize, FID SNR)")
  if (is.null(FID_calibration)) FID_calibration <- str_replace(file_csv,"EAD.csv","FID.csv")
  if (!is.character(FID_calibration)) stop("'FID_calibration' must be character")
  if (!file.exists(FID_calibration)) stop("FID calibration file doesn't exist")
  if (str_ends(file_csv,"EAD.csv")) title_file <- str_remove(file_csv,"_EAD.csv")
  if (!is.numeric(skip_pk_cal_GCEAD)) if(!is.na(skip_pk_cal_GCEAD)) stop("'skip_pk_cal_GCEAD' must be a integer (length(1)) or NA")
  if (!is.numeric(skip_pk_cal_FID)) if(!is.na(skip_pk_cal_FID)) stop("'skip_pk_cal_FID' must be a integer (length(1)) or NA")

  # import GC-EAD ####
  gcead <- read.csv(paste0(wd,"/",file_csv))[,(3*num):(3*num+1)]
  colnames(gcead) <- c("EAD","FID")
  nrG <- nrow(gcead)

  # modify delay
  if(delay > 0) gcead$EAD[(delay+1):nrG] <- gcead$EAD[1:(nrG-delay)]
  if(delay < 0) gcead$EAD[1:(nrG+delay)] <- gcead$EAD[(1-delay):nrG]

  # FID saturation corection ####
  # corrected by the median
  if(median(gcead$FID) <0) gcead$FID <- gcead$FID + abs(median(gcead$FID))*1.2

  # corrected the saturation windows
  fmr <- which(gcead$FID <0)
  ind_sat <- unique(c(fmr,fmr-20, fmr+20)) %>% sort()
  gcead$FID[which(gcead$FID <0)] <- max(gcead$FID)*0.99

  # # scale EAD and FID
  # gcead <- scale(gcead) %>% as.data.frame()

  # import FID calibration ####
  fid <- read.csv(paste0(wd,"/",FID_calibration))

  colnames(fid) <- c("time_ms","time","RI","FID")
  nrF <- nrow(fid)
  lowlimF <- match(fit.exact.value(skip_time, fid$time),fid$time)
  xFid <- lowlimF:nrF
  vF <- fid$FID[xFid]

  # find top five peaks GC-EAD ####
  xEad <- round(range(fid$time_ms[xFid])/10,0)
  xEad <- xEad[1]:min(xEad[2],nrG)
  vGn <- gcead$FID[xEad]
  vG <- vGn - min(vGn)
  pkGC <- createMassSpectrum(1:length(vG), vG)
  pkGC <- detectPeaks(pkGC, halfWindowSize = param[1], method = "MAD", SNR = param[2])

  if(is.na(skip_pk_cal_GCEAD)){
    maxPK <- which(pkGC@mass %in% (ind_sat-xEad[1]))
    n_skip <- 0
  }
  if(!is.na(skip_pk_cal_GCEAD)){
    maxPK <- 0
    n_skip <- skip_pk_cal_GCEAD
  }

  if(length(maxPK)>0){
    pkGC <- data.frame(mass = pkGC@mass[-maxPK],intensity =  pkGC@intensity[-maxPK])
  } else {
    pkGC <- data.frame(mass = pkGC@mass,intensity =  pkGC@intensity)
  }
  top_pkG <- sort(pkGC$intensity,decreasing = TRUE)[(1:5)+n_skip] %>% match(pkGC$intensity)

  # find top five peaks FID ####
  if(!is.na(skip_pk_cal_FID)) n_skip <- skip_pk_cal_FID
  if(is.na(skip_pk_cal_FID)){
    if(is.na(skip_pk_cal_GCEAD)) n_skip <- length(maxPK) %>% divide_by(2) %>% ceiling()
    if(!is.na(skip_pk_cal_GCEAD)) n_skip <- skip_pk_cal_GCEAD
  }

  pkFID <- createMassSpectrum(1:length(vF), vF)
  pkFID <- detectPeaks(pkFID, halfWindowSize = param[3], method = "MAD", SNR = param[4])
  top_pkF <- sort(pkFID@intensity,decreasing = TRUE)[(1:5)+n_skip] %>% match(pkFID@intensity)

  # calibration phase ####
  cal <- data.frame(time = fid$time[(nrF/3-1)+pkFID@mass[top_pkF]],
                    index = round(xEad[2]+pkGC$mass[top_pkG],0))

  res_reg <- lm(time ~ index,data = cal)

  gcead$time <- (1:nrG)*res_reg$coefficients[2]+res_reg$coefficients[1]

  # export figure control ####
  fmr <- str_split(title_file,pattern = "/")[[1]]
  title_nm <- fmr[length(fmr)]

  tiff(paste0(wd,"/figures/", title_nm,"n",num,".tiff"), width = 600, height = 900, units = "px", res=NA)
   par(mfrow=c(3,1), oma = c(0,0,2,0), mar = c(3,3,2,0), mgp = c(2,0.5,0),
       cex = 1.5)
      matplot(fid$time[xFid],vF, type ="l",main = "calibrated FID",
              xlab = "time (min)", ylab = "intensity (u.a.)")
      points(fid$time[(lowlimF-1)+pkFID@mass[top_pkF]],pkFID@intensity[top_pkF],
             pch=16,col="blue")
      mtext(text = paste(title_file,"n#",num), side = 3,line = 2.5, cex = 2)

      matplot(xEad, vG, type ="l",main = "uncalibrated GC-EAD",
              xlab = "index", ylab = "intensity (u.a.)")
      points(xEad[1]+pkGC$mass[top_pkG]-1,pkGC$intensity[top_pkG],pch=16,col="red")

      rX <- range(fid$time[(nrF/3):nrF])
      rY <- c(median(gcead$FID),max(pkGC$intensity))
      matplot(gcead$time,gcead$FID, type ="l",main = "calibrated GC-EAD",
              xlab = "time (min)", ylab = "intensity (u.a.)", xlim = rX, ylim = rY)
      points(gcead$time[xEad[2]+pkGC$mass[top_pkG]],pkGC$intensity[top_pkG]+min(vGn),
             pch=16,col="green")
  dev.off()

  # FID calibration transfert ####

  ## a filter centered on L0, and with a FWHM of B is given by:
  ##     y/ymax = exp( -4*ln(2)*( (L-L0)/B )^2 )
  ##     FID_cal = -4*ln(2) / B^2

  # initialize vector
  X0 <- gcead$time * 0
  rt <- fid$time[xFid]

  # detect limit
  zi <- range(rt) %>% sapply(fit.exact.value,gcead$time) %>% sapply(match,gcead$time)

  # calculate filter :
  for(i in zi[1]:zi[2]){
    ind <- fit.exact.value(gcead$time[i], rt) %>% sapply(match,rt) # index of rt signal
    X <- subtract(rt,gcead$time[i]) # centred
    X <- exp(FID_cal*X^2) * vF # gaussian windows applied on signal
    X0[i] <- sum(X) # value
  }
  X0 <- X0 * max(vF)/max(X0) # normalized

  if(save_fid_cal == TRUE){
    df_cal <- data.frame(rt, vF)
    df_fid <- data.frame(rt2 = gcead$time, fid = X0)

    fig <- plot_ly(df_cal, x = ~rt, y = ~vF, type = 'scatter', mode = 'lines', name = "FID cal")
    fig <- fig %>% add_trace(data = df_fid, x = ~rt2, y = ~fid, mode = 'lines', name = "FID flt")

    htmlwidgets::saveWidget(fig, paste0(wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
    file.rename(from = paste0(wd,"/figures/options_graph_dy.html"),
                to = paste0(wd,"/figures/",title_nm,".html"))
  }

  # scale new FID
  X0 <- X0*10/max(X0)
  gcead$FID <- as.numeric(X0) - min(X0)

  # recalc peaks ####
  pkGC <- createMassSpectrum(1:nrG, gcead$FID)
  pkGC <- detectPeaks(pkGC, halfWindowSize = param[3], method = "MAD", SNR = param[4])
  pkGC <- data.frame(mass = pkGC@mass,intensity =  pkGC@intensity)

  # happy end ####
  cal_use <- rep(TRUE,10)
  cal_use[1:n_skip] <- FALSE
  dec_ead <- order(pkGC$intensity,decreasing = TRUE)  %>% head(10)
  dec_fid <- order(pkFID@intensity,decreasing = TRUE) %>% head(10)

  pk_res <- data.frame(EADmass = pkGC$mass[dec_ead],
                       EAD_int = pkGC$intensity[dec_ead],
                       FIDmass = pkFID@mass[dec_fid],
                       FID_int = pkFID@intensity[dec_fid],
                       pk_time = fid$time[(nrF/3-1)+pkFID@mass[dec_fid]],
                       cal_use = cal_use)

  return(new(Class = "gcead", GC_EAD = gcead,
                              delay = delay,
                              file = title_file,
                              num = num,
                              wd = wd,
                              pk_res = pk_res,
                              type = "raw"))
}

#' allows to recalculate delay between FID and EAD signal
#'
#' @param pk_EAD the time of a depolarizing EAD peak
#' @param pk_FID the time of a trigger FID peak
#' @param pk_mat a gcead object
#'
#' @return the delay between FID and EAD signal.
#' @export
#'
#' @examples
#' # recalc.delay(pk_FID = 13.688, pk_EAD = 13.623, dl_init = 0) # 386 # dl_init = delai initial
#' # 386
#'
#' # phase 2
#' # gcead <- import.GC.EAD(file_csv = "data/sbe18.csv", num = 1, delay = 386, FID_calibration = "data/sbe21_FID2.csv")
#' # gcead <- print.GCEAD()
recalc.delay <- function(pk_EAD, pk_FID, pk_mat = gc_ead){

  # check ####
  if (class(pk_mat)[1] != "gcead") stop("'pk_mat' must be a gcead S4 object")
  if (!is.numeric(pk_EAD)) stop("'pk_EAD' must be a numeric")
  if (!is.numeric(pk_FID)) stop("'pk_FID' must be a numeric")

  pk_EAD <- fit.exact.value(pk_EAD, pk_mat@GC_EAD$time)
  pk_FID <- fit.exact.value(pk_FID, pk_mat@GC_EAD$time)

  fmr <- subtract(match(pk_EAD, pk_mat@GC_EAD$time), match(pk_FID, pk_mat@GC_EAD$time))
  return(pk_mat@delay - fmr)
}

# 02 Print ####

#' print gcead object
#'
#' @param pk_mat a gcead object with spectra and metadata
#'
#' @return a dynamic graph
#' @export
#'
#' @examples
#' # gcead <- gcead.graph()
gcead.graph <- function(pk_mat = gc_ead){

  # check #
  if (class(pk_mat)[1] != "gcead") stop("'pk_mat' must be a gcead S4 object")

  # graph
  dyG  <- pk_mat@GC_EAD[,c("time","EAD","FID")] %>%
            dygraph(main = paste(pk_mat@file,"n#",pk_mat@num)) %>%
            dyAxis("x", label = "Time and mass") %>%
            dyAxis("y", label = "Intensity (a.u.)") %>%
            dyOptions(axisLineColor = "navy",
                      gridLineColor = "lightblue",
                      useDataTimezone = FALSE) %>%
            dyRangeSelector() %>%
            dyUnzoom()
  print(dyG)
}

#' print mead object
#'
#' @param pk_mat a mead object with spectra and metadata
#' @param view_raw logical. Add raw EAD if TRUE
#' @param save_htlm logical. save the graphe in 'Figures' folder
#' @param view_snr logical. view the SNR (signal noise ratio) associted with mean or median.
#' @param view_selection characters. Selected samples or treatments. Argument can be unique (="all") or multiple (= c("mean","median")).
#' @param prefixe characters. A prefixe for the title
#'
#' @return a dynamic graph
#' @export
#'
#' @examples
#' # list_ead <- gcead.graph(list_ead)
mead.graph <- function(pk_mat = list_ead, view_raw = TRUE, view_selection = c("all"),
                       view_snr = FALSE, save_htlm = FALSE, prefixe = ""){

  # check ####
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")
  if (!is.logical(view_raw)) stop("'view_raw' must be a logical")
  if (!is.logical(view_snr)) stop("'view_snr' must be a logical")
  if ((view_raw == TRUE)&(view_snr == TRUE)) stop("'view_raw' and 'view_snr' can't be TRUE in same time")
  if (!is.logical(save_htlm)) stop("'save_htlm' must be a logical")
  if (!is.character(view_selection[1])) stop("'view_selection' must be a character")
  if (!is.character(prefixe)) stop("'prefixe' must be a character")
  if(("figures" %in% dir(pk_mat@wd))==FALSE) dir.create(paste0(pk_mat@wd,"/figures"))
  view_selection <- str_replace(view_selection, "mean","ave")
  view_selection <- str_replace(view_selection, "average","ave")
  view_selection <- str_replace(view_selection, "median","med")

  # delete unuses SNR
  if(view_snr == FALSE){ # supprime tous snr de cr_EAD
    fmr <- grep("snr",names(pk_mat@cr_EAD))
    if(length(fmr) >0) pk_mat@cr_EAD <- pk_mat@cr_EAD[,-fmr]
  }

  if((view_snr == TRUE) & view_selection[1] != "all"){
    fmr <- sapply(view_selection, grep, names(pk_mat@cr_EAD)) %>% unlist()
    fmr <- grep("snr",names(pk_mat@cr_EAD)) %>% setdiff(fmr)
    if(length(fmr) >0) pk_mat@cr_EAD <- pk_mat@cr_EAD[,-fmr]
  } # supprime ceux non selection

  # color init ####
  ns <- max(ncol(pk_mat@EAD),ncol(pk_mat@cr_EAD))
  col_init <- c(brewer.pal(7,name = "Dark2"),brewer.pal(8,name = "Set1")[-6],brewer.pal(7,name = "Accent")[-4])
  col_init <- divide_by(ns,20) %>% ceiling() %>% rep(col_init,.)

  # initiation
  if(ncol(pk_mat@bl_EAD) == 0){
    if (view_raw == FALSE) stop("'view_raw' can't is FALSE whithout EAD corrected")
    dy_df <- data.frame(time = pk_mat@time, FID = pk_mat@FID, pk = pk_mat@pk_windows) %>%
              cbind(pk_mat@EAD)
    col <- c(col_init[1:ns],brewer.pal(3,name = "Paired")[1:2])
  }

  if(ncol(pk_mat@bl_EAD) > 0){
    dy_df <- data.frame(time = pk_mat@time,
                        FID = pk_mat@FID,
                        pk = pk_mat@pk_windows)
    if(view_raw == TRUE){
      dy_df <- cbind(dy_df, pk_mat@EAD, pk_mat@bl_EAD, pk_mat@cr_EAD)
      col <- c(rep(col_init[1:ns],3),brewer.pal(3,name = "Paired")[1:2])
    }
    if(view_raw == FALSE){
      dy_df <- cbind(dy_df, pk_mat@cr_EAD)
      col <- c(col_init[1:ns], brewer.pal(3,name = "Paired")[1:2])
    }
  }

  # le titre
  title <- str_flatten(pk_mat@name,collapse = ", ")

  # selection
  if(view_selection[1] != "all"){ # view_selection = c("mean","median")
    for(i in 1:length(view_selection)) if(length(grep(view_selection[i],names(dy_df)))==0) stop(" 'selection' is uncorrect")
    fmr <- names(dy_df) %>% str_remove_all("EAD_")
    if(view_snr == TRUE) view_selection <- paste0("snr_",view_selection) %>% intersect(fmr) %>% c(view_selection,.)
    keep_that <- sapply(view_selection, match, table = fmr)
    dy_df <- dy_df[,c(1:3,keep_that)]
    lc <- length(col)
    col <- col[c(keep_that, lc-1, lc)]
  }

  # Creation du graphe
  {
  dyG  <- dygraph(dy_df, main = title) %>%
    dySeries("FID", label = "FID_signal", strokeWidth = 2, brewer.pal(3,name = "Paired")[1]) %>%
    dySeries("pk", label = "FID_peak", strokeWidth = 2, brewer.pal(3,name = "Paired")[2]) %>%
    dyAxis("x", label = "Time and mass") %>%
    dyAxis("y", label = "Intensity (a.u.)") %>%
    dyOptions(axisLineColor = "navy", gridLineColor = "lightblue",
              useDataTimezone = FALSE, colors = col, strokeWidth =1) %>%
    dyLegend(show = "onmouseover", hideOnMouseOut = TRUE) %>%
    dyRangeSelector() %>% dyUnzoom() } # initiation

  {
    # Print SD correctly
    if(("raw_EAD_sd" %in% dyG$x$attrs$labels) == TRUE){
      dyG <- dySeries(dyG,"raw_EAD_sd", label = "raw_EAD_sd", strokeWidth = 2, color = "#188935")
    }
    if(("EAD_sd" %in% dyG$x$attrs$labels) == TRUE){
      dyG <- dySeries(dyG,"EAD_sd", label = "EAD_sd", strokeWidth = 2, color = "#188935")
    }
    if(("bl_EAD_sd" %in% dyG$x$attrs$labels) == TRUE){
      dyG <- dySeries(dyG,"bl_EAD_sd", label = "bl_EAD_sd", strokeWidth = 2, color = "#188935")
    }

    # Print SNR med correctly
    if(("EAD_snr_med" %in% dyG$x$attrs$labels) == TRUE) dyG <- dySeries(dyG,"EAD_snr_med",
                                                                        label = "EAD_snr_med",
                                                                        strokeWidth = 1.5,
                                                                        color = "#fbb4ae")

    # Print SNR med correctly
    if(("EAD_snr_ave" %in% dyG$x$attrs$labels) == TRUE) dyG <- dySeries(dyG,"EAD_snr_ave",
                                                                        label = "EAD_snr_ave",
                                                                        strokeWidth = 1.5,
                                                                      color = "#fed9a6")
    # Print average/mean correctly
    if(("raw_EAD_ave" %in% dyG$x$attrs$labels) == TRUE){
      dyG <- dySeries(dyG,"raw_EAD_ave", label = "raw_EAD_ave", strokeWidth = 3, color = "#1A1A1A")
    }
    if(("EAD_ave" %in% dyG$x$attrs$labels) == TRUE){
      dyG <- dySeries(dyG,"EAD_ave", label = "EAD_ave", strokeWidth = 3, color = "#1A1A1A")
    }
    if(("bl_EAD_ave" %in% dyG$x$attrs$labels) == TRUE){
      dyG <- dySeries(dyG,"bl_EAD_ave", label = "bl_EAD_ave", strokeWidth = 3, color = "#1A1A1A")
    }

    # Print median correctly
    if(("raw_EAD_med" %in% dyG$x$attrs$labels) == TRUE){
      dyG <- dySeries(dyG,"raw_EAD_med", label = "raw_EAD_med", strokeWidth = 3, color = "#666666")
    }
    if(("EAD_med" %in% dyG$x$attrs$labels) == TRUE){
      dyG <- dySeries(dyG,"EAD_med", label = "EAD_med", strokeWidth = 3, color = "#666666")
    }
    if(("bl_EAD_med" %in% dyG$x$attrs$labels) == TRUE){
      dyG <- dySeries(dyG,"bl_EAD_med", label = "bl_EAD_med", strokeWidth = 3, color = "#666666")
    }

    } # gestion des series

  {
  tooltips <- list()
  for(i in 1:ncol(pk_mat@peak)){
    fmr <- list(list(x = pk_mat@time[pk_mat@peak[1,i]],
                     text = pk_mat@time[pk_mat@peak[1,i]],
                     series = "FID_signal",
                     tooltip = colnames(pk_mat@peak)[i], width=50, tickHeight = 0.1))
    tooltips <- c(tooltips, fmr)
  }

  annotator <- function(x,y){
    d = do.call(dyAnnotation,modifyList(list(dygraph=x),y))
    return(d)
  }

  dyG <- Reduce(annotator, tooltips, init=dyG )

  dyG$x$css = ".dygraphDefaultAnnotation {
                background-color: transparent;
                border: none;
                tickColor: transparent !important;
                color: black !important;
                width: initial !important;
                font-size: 60% !important;
              }"
  } # gestion de l'affiche de pic

  print(dyG)

  if(save_htlm == TRUE){

    title <- str_replace_all(title," ","_") %>% str_replace_all(",","_") %>% str_replace_all("__","_")
    if(prefixe != "") title <- paste0(prefixe,"_",title)
    fmr <- ""
    if(view_raw == TRUE) fmr <- "_raw"
    if(view_selection[1] != "all"){
      if(view_snr == TRUE) view_selection <- view_selection[-grep("snr",view_selection)]
      fmr <- str_flatten(c("",view_selection),"_") %>% c(fmr) %>% str_flatten()
    }

    htmlwidgets::saveWidget(dyG, paste0(pk_mat@wd,"/figures/options_graph_dy.html"), selfcontained = TRUE)
    file.rename(from = paste0(pk_mat@wd,"/figures/options_graph_dy.html"),
                to = paste0(pk_mat@wd,"/figures/",title,fmr,"_dyGraph.html"))
  } # save htlm

}

# 03 Pretreatments ####

#' modify an existing FID peak
#'
#' @param pk time of the peak
#' @param bm_inf the new lower limit
#' @param bm_sup the new upper limit
#' @param pk_mat a mead object
#'
#' @return a mead object
#' @export
#'
#' @examples
#' # list_ead <- FID.modify.peak(pk = 11.804, bm_inf = 11.715)
#' # list_ead <- FID.modify.peak(pk = 11.612, bm_inf = 11.526)
#' # list_ead <- FID.modify.peak(pk = 11.453 , bm_inf = 11.382)
#' # list_ead <- FID.modify.peak(pk = 11.292, bm_inf = 11.226)
FID.modify.peak <- function(pk, bm_inf = NA, bm_sup = NA, pk_mat = list_ead){

  # check ####
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")
  if (!is.numeric(pk)) stop("'pk' must be a numeric")
  if (!is.numeric(bm_inf)) if(!is.na(bm_inf)) stop("'bm_inf' must be a numeric or NA")
  if (!is.numeric(bm_sup)) if(!is.na(bm_sup)) stop("'bm_sup' must be a numeric or NA")

  # find exact value for pk, bm_inf and bm_sup
  pk <- colnames(pk_mat@peak) %>% as.numeric() %>% fit.exact.value(pk,.)
  if(is.na(bm_inf) == FALSE) bm_inf <- fit.exact.value(bm_inf, pk_mat@time)
  if(is.na(bm_sup) == FALSE) bm_sup <- fit.exact.value(bm_sup, pk_mat@time)

  # find the index, the real power of the matrix in R, and replace new limit
  fmr <- match(pk, colnames(pk_mat@peak))
  if(is.na(bm_inf) == FALSE) pk_mat@peak[2,fmr] <- match(bm_inf,pk_mat@time)
  if(is.na(bm_sup) == FALSE) pk_mat@peak[3,fmr] <- match(bm_sup,pk_mat@time)

  # peak windows
  pkz <- rep(0,length(pk_mat@time))
  for (i in 1:ncol(pk_mat@peak)) pkz[pk_mat@peak[2,i]:pk_mat@peak[3,i]] <- pk_mat@FID[pk_mat@peak[1,i]]
  pk_mat@pk_windows <- pkz

  return(pk_mat)
}

#' create new FID peak
#'
#' @param bm_inf the lower limit
#' @param bm_sup the upper limit
#' @param pk_mat a mead object
#'
#' @return a mead object
#' @export
#'
#' @examples
#' # list_ead <- FID.create.peak(pk = 10.756, bm_inf = 10.73, bm_sup = 10.8)
FID.create.peak <- function(bm_inf, bm_sup, pk_mat = gc_ead){

  # check ####
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")
  if (!is.numeric(bm_inf)) stop("'bm_inf' must be a numeric")
  if (!is.numeric(bm_sup)) stop("'bm_sup' must be a numeric")

  # find exact value
  bm_inf <- fit.exact.value(bm_inf, pk_mat@time)
  bm_sup <- fit.exact.value(bm_sup, pk_mat@time)

  # find the index
  ind_inf <- match(bm_inf,pk_mat@time)
  ind_sup <- match(bm_sup,pk_mat@time)
  ind_pk <- mean(c(ind_inf,ind_sup)) %>% round(0)

  # find the peak
  int_pk <- max(pk_mat@FID[ind_inf:ind_sup])

  pk_mat@peak <- cbind(pk_mat@peak, c(ind_pk,ind_inf,ind_sup))
  colnames(pk_mat@peak)[ncol(pk_mat@peak)] <- pk_mat@time[ind_pk]
  pk_mat@peak <- pk_mat@peak[,order(pk_mat@peak[1,])]

  # peak windows
  pkz <- rep(0,length(pk_mat@time))
  for (i in 1:ncol(pk_mat@peak)) pkz[pk_mat@peak[2,i]:pk_mat@peak[3,i]] <- pk_mat@FID[pk_mat@peak[1,i]]
  pk_mat@pk_windows <- pkz

  return(pk_mat)
}

#' delete FID peak
#'
#' @param pk RT of peak
#' @param pk_mat a mead object
#'
#' @return a mead object
#' @export
#'
#' @examples
#' # list_ead <- FID.delete.peak(pk = 9.74254)
FID.delete.peak <- function(pk, pk_mat = gc_ead){

  # check ####
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")
  if (!is.numeric(pk)) stop("'pk' must be a numeric")

  # find exact value
  pk <- names(pk_mat@peak) %>% as.numeric() %>% fit.exact.value(pk,.)

  # find the index
  ind_pk <- match(pk,names(pk_mat@peak))
  pk_mat@peak <- pk_mat@peak[,-ind_pk]

  # peak windows
  pkz <- rep(0,length(pk_mat@time))
  for (i in 1:ncol(pk_mat@peak)) pkz[pk_mat@peak[2,i]:pk_mat@peak[3,i]] <- pk_mat@FID[pk_mat@peak[1,i]]
  pk_mat@pk_windows <- pkz

  return(pk_mat)
}

#' delete FID window. Allows you to delete the start or end of a run
#'
#' @param X1 the lower limit or "start" (=T0)
#' @param X2 the upper limit or "end" (=Tmax)
#' @param pk_mat a mead object
#'
#' @return a mead object
#' @export
#'
#' @examples
#' # list_ead <- FID.delete.window(X2 = 5)
#' # list_ead <- FID.delete.window(X1 = 18)
FID.delete.window <- function(X1 = "start", X2 = "end", pk_mat = list_ead){

  # check ####
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")
  if(X1 == "start" & X2 == "end") stop("Poor boy ! I can't delete all. Re-try with another parameters")

  if (!is.numeric(X1)) if (X1 != "start") stop("'X1' must be a numeric")
  if(X1 == "start") X1 <- pk_mat@time[1]
  if (!is.numeric(X2)) if (X2 != "end") stop("'X2' must be a numeric")
  if(X2 == "end") X2 <- pk_mat@time[length(pk_mat@time)]

  # find index
  X1 <- fit.exact.value(X1, pk_mat@time) %>% match(pk_mat@time)
  X2 <- fit.exact.value(X2, pk_mat@time) %>% match(pk_mat@time)

  # cut the window
  pk_mat@EAD <- pk_mat@EAD[-(X1:X2),]
  if(ncol(pk_mat@bl_EAD)!=0) pk_mat@bl_EAD <- pk_mat@bl_EAD[-(X1:X2),]
  if(ncol(pk_mat@cr_EAD)!=0) pk_mat@cr_EAD <- pk_mat@cr_EAD[-(X1:X2),]

  pk_mat@FID <- pk_mat@FID[-(X1:X2)]
  pk_mat@time <- pk_mat@time[-(X1:X2)]
  pk_mat@pk_windows <- pk_mat@pk_windows[-(X1:X2)]

  # deletee peak in the window
  ind <- pk_mat@peak[1,pk_mat@peak[1,] %in% X1:X2]
  if(length(ind) >0) pk_mat@peak <- pk_mat@peak[,-match(ind, pk_mat@peak[1,])]

  # recalc index of peak and peak limits
  fmr <- match(colnames(pk_mat@peak), pk_mat@time) %>% subtract(pk_mat@peak[1,])
  pk_mat@peak[1,] <- add(fmr,pk_mat@peak[1,])
  pk_mat@peak[2,] <- add(fmr,pk_mat@peak[2,])
  pk_mat@peak[3,] <- add(fmr,pk_mat@peak[3,])

  return(pk_mat)
}

# 05 Multiple EAD ####

#' merge a lot of gcead objects in a m(ultiple-gc-)ead
#'
#' @param ... a suit of gcead objects
#' @param gcead_names a vector with new name (optionnal)
#' @param RT_limits a couple of decimals for delimit the boundaries of RT (egg : RT_limits = c(4.5,16))
#' @param width_smooth width smotthing for "fl" procma::savgol function
#' @param wm_bl_fid Width of local window for minimization/maximization for baseline::baseline.rollingBall function
#' @param ws_bl_fid Width of local window for smoothing for baseline::baseline.rollingBall function
#' @param th_pk_fid threshold of FID peak intensity
#' @param hws_pk_fid numeric, half window size for MALDIquant::detectPeaks
#' @param snr_pk_fid single numeric value for MALDIquant::detectPeaks. SNR is an abbreviation for signal-to-noise-ratio. A local maximum has to be higher than SNR*noise to be recognize as peak.
#' @param gap_FID numeric. Multiplied intensity of FID signal.
#'
#' @return
#' @export
#'
#' @examples
#' # list_test <- gcead.merge(b23, b24, b26, b29, b30, b31, gcead_names = c("b23","b24","b26","b29","b30","b31"))
gcead.merge <- function(..., gcead_names = TRUE, RT_limits = NULL,
                        width_smooth = 125, wm_bl_fid = 400,
                        ws_bl_fid = 10, th_pk_fid = 0.01,
                        hws_pk_fid = 10, snr_pk_fid = 0.5,
                        gap_FID = 1){

  # merge all
  # if(typeof(...)[1] == "all"){
  #   fmr <- which(sapply(ls(), function(x) class(get(x))[1]) == "gcead")
  #   list_gcead <- lapply(ls()[fmr],get)
  # }else{
     list_gcead <- list(...) # list_gcead <- list(a04, a05, a07, a08, a09, a10, a11, a12, a14, a15)
  # }


  # check
  sapply(list_gcead, function(X) if (class(X) != "gcead") stop("variables must be a gcead S4 object"))
  if (!is.logical(gcead_names)) if (!is.character(gcead_names))  stop("variables must be a gcead S4 object")
  if (!is.numeric(width_smooth)) stop("'width_smooth' must be a numeric")
  if (!is.numeric(wm_bl_fid)) stop("'wm_bl_fid' must be a numeric")
  if (!is.numeric(ws_bl_fid)) stop("'ws_bl_fid' must be a numeric")
  if (!is.numeric(th_pk_fid)) stop("'th_pk_fid' must be a numeric")
  if (!is.numeric(hws_pk_fid)) stop("'hws_pk_fid' must be a numeric")
  if (!is.numeric(snr_pk_fid)) stop("'snr_pk_fid' must be a numeric")
  if (!is.numeric(gap_FID)) stop("'gap_FID' must be a single value")
  if (length(gap_FID) != 1) stop("'gap_FID' must be a single value")

  if (!is.null(RT_limits)) if(length(RT_limits)!=2) if(!is.numeric(RT_limits)) stop("'RT_limits' must be either a couple of numeric either null.")

  # rename
  if(gcead_names[1] == FALSE){
    nm_ls <- sapply(list_gcead, function(X) X@file)
  } else if(gcead_names[1] == TRUE){
    nm_ls <- sapply(list_gcead, function(X) X@file)
    nm_ls <- str_split(nm_ls, "/",simplify = TRUE)
    nm_ls <- str_split(nm_ls[,ncol(nm_ls)], "_",simplify = TRUE)
    if (class(nm_ls)[1] != "matrix") stop("samples names must be uniform (same number of underscores)")
    fmr <- apply(nm_ls,2, function(X) length(unique(X)))
    nm_ls <- nm_ls[,-which(fmr == 1)] %>% apply(1, str_flatten, collapse = "_")
  } else if(!is.logical(gcead_names)){
    nm_ls <- gcead_names
  }

  ## ALIGNEMENT
  # Align preparation
  list_pk <- lapply(list_gcead, function(X) X@pk_res)
  nb_pk_min <- sapply(list_pk, function(X) nrow(X)) %>% min() # nombre de pic minimum

  ls_pk <- list()
  lg <- length(list_gcead)
  align_pk <- matrix(NA, nb_pk_min, 4, dimnames = list(paste0("pk",1:nb_pk_min),c("sd","detect","num","Tmoy")))

  # Peak by peak
  for (i in 1:nb_pk_min){
    mat_pk <- matrix(NA,lg+3,6)
    rownames(mat_pk) <- c(paste0("pk",i,"_",nm_ls),"mean","median","sd")
    colnames(mat_pk) <- colnames(list_pk[[1]])
    for(j in 1:lg) mat_pk[j,] <- t(list_pk[[j]][i,])
    mat_pk[lg+1,] <- apply(mat_pk[1:lg,],2, mean, na.rm = TRUE)    # calcul des moyennes
    mat_pk[lg+2,] <- apply(mat_pk[1:lg,],2, median, na.rm = TRUE)  # calcul des medianes
    mat_pk[lg+3,] <- apply(mat_pk[1:lg,],2, sd, na.rm = TRUE)      # calcul des écart-types
    ls_pk[[i]] <- as.data.frame(mat_pk)
    align_pk[i,] <- c(mat_pk["sd",5],mat_pk["mean",6],i,mat_pk["mean",5]) # un resumé pour après
  }

  # detection du pic pivot pour l'alignement
  align_pk <- align_pk[which(align_pk[,2]==1),]  # suppresion des pics non detectés partout
  num_pk_align <- align_pk[which.min(align_pk[,1]),"num"] # le pic le mieux aligné

  # detection des index pour chaque echantillon
  ind_pk <- ls_pk[[num_pk_align]][1:lg,1]
  ind_start <- ind_pk - min(ind_pk) + 1
  max_end <- sapply(list_gcead, function(X,nm) nrow(X@GC_EAD)-X@pk_res[nm,1],nm = num_pk_align) %>% min()
  ind_end <- ind_pk + max_end
  coor_mat <- apply(rbind(ind_start,ind_end),2, function(X) X[1]:X[2])

  # mise en forme des matrices EAD, FID, time
  EAD <- matrix(NA,lg,nrow(coor_mat))
  rownames(EAD) <- paste0("raw_EAD_",nm_ls)
  FID <- EAD
  rownames(FID) <- paste0("FID_",nm_ls)
  mTime <- EAD
  rownames(mTime) <- paste0("Time_",nm_ls)
  maxFID <- 1:lg

  for(i in 1:lg){
    EAD[i,] <- list_gcead[[i]]@GC_EAD$EAD[coor_mat[,i]] # centré sur le pic pivot
    FID[i,] <- list_gcead[[i]]@GC_EAD$FID[coor_mat[,i]] # centré sur le pic pivot
    mTime[i,] <- list_gcead[[i]]@GC_EAD$time[coor_mat[,i]] # centré sur le pic pivot
    maxFID[i] <- min(list_pk[[i]][,6], na.rm = TRUE)*max(list_pk[[i]][,2],na.rm = TRUE) # ne prend que les FID qui n'ont pas saturé
  }
  if(sum(maxFID)==0) maxFID <- lapply(list_pk, function(x) max(x[,4])) %>% unlist()

  # gestion du temps
  mTime <- colMeans(mTime) %>% round(5) # un seul temps

  if(!is.null(RT_limits)){
    RT_limits <- sapply(RT_limits, fit.exact.value, mTime)
    rt <- sapply(RT_limits, match, mTime)
    mTime <- mTime[rt[1]:rt[2]]
    EAD <- EAD[,rt[1]:rt[2]]
    FID <- FID[,rt[1]:rt[2]]
  }

  # préapration de l'objet S4
  pk_mat <- new(Class = "mead",
                EAD = as.data.frame(EAD),
                time = mTime,
                FID = FID[which.max(maxFID),]*gap_FID, # un seul FID
                delay = sapply(list_gcead, function(X) return(X@delay)),
                file = sapply(list_gcead, function(X) return(X@file)),
                num = sapply(list_gcead, function(X) return(X@num)),
                wd = list_gcead[[1]]@wd,
                name = nm_ls)

  pk_mat <- pretreatment.mead(pk_mat, wsm = width_smooth, snr_pk_FID = snr_pk_fid,
                              wm_bl_FID = wm_bl_fid, ws_bl_FID = ws_bl_fid,
                              th_pk_FID = th_pk_fid, hws_pk_FID = hws_pk_fid)

  pk_mat@param <- c(pk_mat@param, gap_FID)
  names(pk_mat@param)[length(pk_mat@param)] <- "gap_FID"

  return(pk_mat)
}

#' performs several pre-processings and peaks detection on mead S4 object
#'
#' @param pk_mat a mead object with spectra and metadata
#' @param wsm width smotthing for "fl" procma::savgol function
#' @param wm_bl_FID Width of local window for minimization/maximization for baseline::baseline.rollingBall function
#' @param ws_bl_FID Width of local window for smoothing for baseline::baseline.rollingBall function
#' @param th_pk_FID threshold of FID peak intensity
#' @param hws_pk_FID numeric, half window size for MALDIquant::detectPeaks
#' @param snr_pk_FID single numeric value for MALDIquant::detectPeaks. SNR is an abbreviation for signal-to-noise-ratio. A local maximum has to be higher than SNR*noise to be recognize as peak.
#'
#' @return
#' @export
#'
#' @examples
#' # list_test <- pretreatment.mead()
pretreatment.mead <- function(pk_mat = list_ead, wsm = 125, wm_bl_FID = 400,
                              ws_bl_FID = 10, th_pk_FID = 0.01,
                              hws_pk_FID = 10, snr_pk_FID = 0.5){
  # check ####
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")
  if (!is.numeric(wsm)) stop("'wsm' must be a numeric")
  if (!is.numeric(wm_bl_FID)) stop("'wm_bl_FID' must be a numeric")
  if (!is.numeric(ws_bl_FID)) stop("'ws_bl_FID' must be a numeric")
  if (!is.numeric(th_pk_FID)) stop("'th_pk_FID' must be a numeric")
  if (!is.numeric(hws_pk_FID)) stop("'hws_pk_FID' must be a numeric")
  if (!is.numeric(snr_pk_FID)) stop("'snr_pk_FID' must be a numeric")

  # smoothing ####
  spc <- data.frame(FID = pk_mat@FID) %>% cbind(t(pk_mat@EAD)) %>%
    sapply(savgol, fl = wsm) # wsm = width smotthing

  pk_mat@EAD <- as.data.frame(spc[,-1])

  # baseline FID ####
  spc <- baseline(t(spc[,1]), wm = wm_bl_FID, ws = ws_bl_FID, method = "rollingBall")@corrected
  spc[which(spc < 0)] <- 0

  # detect peak FID ####
  # detection des pics avec MALDIQuant
  MSobj <- createMassSpectrum(pk_mat@time, as.numeric(spc))
  MSobj <- detectPeaks(MSobj, halfWindowSize = hws_pk_FID,
                       method=c("MAD", "SuperSmoother"), SNR = snr_pk_FID)

  # detection des bornes, minimum local ou inferieur au seuil
  fmr <- which(MSobj@intensity > th_pk_FID)
  det <- matrix(NA,3,length(fmr),dimnames = list(c("index","bm_inf","bm_sup"),c(MSobj@mass[fmr])))
  det[1,] <- match(MSobj@mass[fmr],pk_mat@time)
  fmr <- which((det["index",] < 201)|(det["index",] > (ncol(spc)-202)))
  if(length(fmr) > 0) det <- det[,-fmr]

  th_bm_FID <- th_pk_FID/10 # seuil des bornes egal a un dixieme du seuil des pics

  # pour le premier pic
  fmr <- (det[1,1]-199):det[1,1]
  ind_bm <- c(which.min(spc[fmr]), which(spc[fmr] < th_bm_FID)) %>% sort()
  det[2,1] <- det[1,1] + max(ind_bm) - 200
  # pour le dernier pic
  fmr <- (1+det[1,ncol(det)]):(201+det[1,ncol(det)])
  ind_bm <- c(which.min(spc[fmr]), which(spc[fmr] < th_bm_FID)) %>% sort()
  det[3,ncol(det)] <- det[1,ncol(det)] + min(ind_bm)
  # pour les autres pics
  for(i in 1:(ncol(det)-1)){ # i=1
    fmr <- (1+det[1,i]):det[1,i+1]
    ind_bm <- c(which.min(spc[fmr]), which(spc[fmr] < th_bm_FID)) %>% sort()
    det[3,i] <- det[1,i]+min(ind_bm)
    det[2,i+1] <- det[1,i]+max(ind_bm)
  }

  # creation de la fonction creneau, fenetre de chaque pic
  pkz <- rep(0,length(spc))
  for (i in 1:ncol(det)) pkz[det[2,i]:det[3,i]] <- spc[det[1,i]]

  FID <- data.frame(time = pk_mat@time, FID = as.numeric(spc), pk = pkz)
  FID_pk <- data.frame(time =   FID$time[det[1,]], peak =    det[1,],
                       bm_inf = FID$time[det[2,]], ind_inf = det[2,],
                       bm_sup = FID$time[det[3,]], ind_sup = det[3,],
                       intensity = spc[det[1,]])
  FID_pk$tooltip <- apply(FID_pk,1, function(X) paste(X[3],"to",X[5]))

  # Graphe FID + pk windoms
  dyG  <- dygraph(FID, main = pk_mat@name) %>%
    dySeries("pk", label = "FID peak") %>%
    dyAxis("x", label = "Time and mass") %>%
    dyAxis("y", label = "Intensity (a.u.)") %>%
    dyOptions(axisLineColor = "navy", gridLineColor = "lightblue",
              useDataTimezone = FALSE) %>%
    dyRangeSelector() %>% dyUnzoom()

  tooltips <- list()
  for(i in 1:nrow(FID_pk)){
    fmr <- list(list(x = FID_pk[i,1], text = FID_pk[i,1], series = "FID",
                     tooltip = FID_pk[i,8], width=50, tickHeight = 0.1))
    tooltips <- c(tooltips, fmr)
  }

  annotator <- function(x,y){
    d = do.call(dyAnnotation,modifyList(list(dygraph=x),y))
    return(d)
  }

  dyG <- Reduce( annotator, tooltips, init=dyG )

  dyG$x$css = ".dygraphDefaultAnnotation {
                background-color: transparent;
                border: none;
                tickColor: transparent !important;
                color: black !important;
                width: initial !important;
                font-size: 60% !important;
              }"
  print(dyG)

  # mise en forme
  fmr <- c(wsm, wm_bl_FID, ws_bl_FID, th_pk_FID, hws_pk_FID, snr_pk_FID)

  names(fmr) <- c("width_smooth", "wm_bl_fid", "ws_bl_fid", "th_pk_fid", "hws_pk_fid",
                  "snr_pk_fid")

  pk_mat@FID <- as.numeric(spc)
  pk_mat@pk_windows <- pkz
  pk_mat@peak <- as.data.frame(det)
  pk_mat@param <- fmr

  return(pk_mat)
}

#' normalize multi-EAD signal
#'
#' @param shift for shift EAD signals
#' @param amplitude modulated the amplitude of signals
#' @param overlay overlay EAD
#' @param pk_mat a mead object
#'
#' @return a mead object
#' @export
#'
#' @examples
#' # list_ead <- m.EAD.norm(shift = -1, amplitude = 2)
m.EAD.norm <- function(shift = -1, amplitude = NA, overlay = TRUE, pk_mat = list_ead){

  # check ####
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")
  if (!is.logical(overlay)) stop("'overlay' must be logical")
  if (!is.numeric(shift)) stop("'shift' must be numeric")
  if (!length(shift) == 1) stop("The length of 'shift' must be egal at one. Use 'overlay'")
  if (!is.numeric(amplitude)) if (!is.na(amplitude)) stop("'amplitude' must be numeric")
  if (!length(amplitude) %in% c(1,ncol(pk_mat@EAD))) stop("The length of 'amplitude' must be egal at one or at the number of samples")

  # delete the original shift
  EAD <- sapply(pk_mat@EAD, function(X) X - min(X))

  # normalized
  if((is.na(amplitude[1]) == FALSE) & (length(amplitude) == 1)) amplitude <- rep(amplitude,ncol(EAD))

  if(is.na(amplitude[1]) == FALSE){
    for (i in 1:ncol(EAD)){
      EAD[,i] <- max(EAD[,i]) %>%
                  divide_by(EAD[,i],.) %>%
                   multiply_by(amplitude[i]) %>%
                    subtract(amplitude[i])
    }
  }

  # add shift
  shift <- rep(shift,ncol(EAD))
  if(overlay == TRUE) shift <- shift*(1:ncol(EAD))
  pk_mat@EAD <- t(EAD) %>% add(shift) %>% t() %>% as.data.frame()

  #maj param
  pk_mat@param <- c(pk_mat@param, shift[1])
  names(pk_mat@param)[length(pk_mat@param)] <- "shift"

  pk_mat@param <- c(pk_mat@param, amplitude[1])
  names(pk_mat@param)[length(pk_mat@param)] <- "amplitude"

  return(pk_mat)
}

#' delete multi-EAD baseline
#'
#' @param cr_shift for shift the EAD signals
#' @param double_baseline make a second subtraction of baseline
#' @param pk_mat a mead object
#'
#' @return a mead object
#' @export
#'
#' @examples
#' # list_ead <- m.EAD.baseline(shift = -0.1)
m.EAD.baseline <- function(cr_shift = -0.1, double_baseline = FALSE, pk_mat = list_ead){
  # check ####
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")
  if (!is.numeric(cr_shift)) stop("'cr_shift' must be a numeric")

  # create new empty matrix
  pk_mat@bl_EAD <- pk_mat@EAD
  colnames(pk_mat@bl_EAD) <- str_replace_all(colnames(pk_mat@bl_EAD),"raw_","bl_")
  pk_mat@cr_EAD <- pk_mat@EAD
  colnames(pk_mat@cr_EAD) <- str_remove_all(colnames(pk_mat@EAD),"raw_")

  # first baseline
  for (j in 1:ncol(pk_mat@EAD)){ # j = 1
    # create a baseline
    for (i in 1:ncol(pk_mat@peak)){ # i = 10
      xbm <- pk_mat@peak[2:3,i]
      ybm <- pk_mat@EAD[xbm,j]
      tbm <- pk_mat@time[xbm]

      fmr <- (ybm[2]-ybm[1])/(xbm[2]-xbm[1])
      fmr[2] <- ybm[1]-fmr[1]*xbm[1]
      abs <- 1:length(pk_mat@time)
      aff <- abs*fmr[1]+fmr[2]
      pk_mat@bl_EAD[xbm[1]:xbm[2],j] <- aff[xbm[1]:xbm[2]]
    }

    # bl <- savgol(pk_mat@bl_EAD[,j], fl = 701)
    #
    # df_raw <- data.frame(abs = pk_mat@time ,EAD = pk_mat@EAD[,j], bl_raw = pk_mat@bl_EAD[,j], bl = bl)
    #
    # fig <- plot_ly(df_raw, x = ~abs, y = ~EAD, name = "cr", type = 'scatter', mode = 'lines')
    # fig <- fig %>% add_trace(y= ~bl_raw, name = 'bl', mode = 'lines')
    # fig <- fig %>% add_trace(y= ~bl, name = 'bl', mode = 'lines')
    # fig

    # create the corrected electrogramme
    pk_mat@cr_EAD[,j] <- pk_mat@EAD[,j] - pk_mat@bl_EAD[,j] + cr_shift
  }

  # double baseline
  if(double_baseline == TRUE){
    cr_EAD_1 <- pk_mat@cr_EAD

    bl_EAD_2 <- pk_mat@cr_EAD#*0 + cr_shift

    for (j in 1:ncol(pk_mat@EAD)){ # j = 1
      for (i in 1:ncol(pk_mat@peak)){ # i = 6
        # limit :
        xbm <- pk_mat@peak[2:3,i] # pk_mat@time[xbm]
        # peak center :
        ybm <- pk_mat@FID[xbm]
        fmr <- (ybm[2]-ybm[1])/(xbm[2]-xbm[1])
        fmr[2] <- ybm[1]-fmr[1]*xbm[1]
        abs <- 1:length(pk_mat@time)
        aff <- abs*fmr[1]+fmr[2]

        fmr <- pk_mat@FID[xbm[1]:xbm[2]] - aff[xbm[1]:xbm[2]]
        appex <- xbm[1] + which.max(fmr) - 1 # pk_mat@time[appex]

        # detect min inf
        xbm[1] <- xbm[1]+which.max(pk_mat@cr_EAD[xbm[1]:appex,j])-1
        # detect min sup
        xbm[2] <- appex+which.max(pk_mat@cr_EAD[appex:xbm[2],j])-1
                                                      ## pk_mat@time[xbm]

        # calcul affine fonction
        ybm <- pk_mat@cr_EAD[xbm,j]
        tbm <- pk_mat@time[xbm]

        fmr <- (ybm[2]-ybm[1])/(xbm[2]-xbm[1])
        fmr[2] <- ybm[1]-fmr[1]*xbm[1]
        abs <- 1:length(pk_mat@time)
        aff <- abs*fmr[1]+fmr[2]
        bl_EAD_2[xbm[1]:xbm[2],j] <- aff[xbm[1]:xbm[2]]
      }
      # create the corrected electrogramme
      pk_mat@cr_EAD[,j] <- pk_mat@cr_EAD[,j] - bl_EAD_2[,j] + cr_shift
    }

    cr_EAD_2 <- pk_mat@cr_EAD
    df_cr <- data.frame(abs = pk_mat@time,
                        cr_1 = apply(cr_EAD_1,1,mean),
                        cr_2 = apply(cr_EAD_2,1,mean),
                        fid = pk_mat@FID,
                        fid_wd = pk_mat@pk_windows)

    fig <- plot_ly(df_cr, x = ~abs, y = ~cr_1, name = "cr_1", type = 'scatter', mode = 'lines')
    fig <- fig %>% add_trace(y= ~cr_2, name = 'cr_2', mode = 'lines')
    fig <- fig %>% add_trace(y= ~fid, name = 'fid', mode = 'lines')
    fig <- fig %>% add_trace(y= ~fid_wd, name = 'fid_wd', mode = 'lines')
    print(fig)

  } # double_baseline = TRUE

  fmr <- sapply(c("_ave","_med","_sd"), grep, names(pk_mat@cr_EAD)) %>% unlist()
  if(length(fmr)==0) fmr = ncol(pk_mat@EAD) + 1  # pour ne supprimer aucune colonne
  pk_mat@EAD$raw_EAD_sd <- apply(pk_mat@EAD[,-fmr], 1, sd)
  pk_mat@bl_EAD$bl_EAD_sd <- apply(pk_mat@bl_EAD[,-fmr], 1, sd)
  pk_mat@cr_EAD$EAD_sd <- apply(pk_mat@cr_EAD[,-fmr], 1, sd)

  pk_mat@param <- c(pk_mat@param, cr_shift)
  names(pk_mat@param)[length(pk_mat@param)] <- "cr_shift"

  return(pk_mat)
}

#' calculate an average for the multi-EAD
#'
#' @param pk_mat a mead object
#' @param treatment the type of treatment. Choose "mean" or "median"
#' @param shift for shift EAD signals
#' @param overlay overlay EAD
#'
#' @return a mead object
#' @export
#'
#' @examples
#' # list_ead <- m.EAD.average(list_ead,"mean")
m.EAD.average <- function(pk_mat = list_ead, treatment = c("mean","median"),
                          shift = -1, overlay = TRUE){

  # check
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")
  if (class(treatment)[1] != "character") stop("'treatment' must be a character")
  if (ncol(pk_mat@cr_EAD) == 0) stop("calculates the baseline before")
  if (!is.logical(overlay)) stop("'overlay' must be logical")
  if (!is.numeric(shift)) stop("'shift' must be numeric")
  if (!length(shift) == 1) stop("The length of 'shift' must be egal at one. Use 'overlay'")

  # suffixe
  suf <- c("_ave","_med","_sd")
  suffixe <- suf[grep(treatment[1], c("mean","median"))]
  fmr <- grep(suffixe, colnames(pk_mat@EAD))
  if(length(fmr) > 0){
    pk_mat@EAD <- pk_mat@EAD[-fmr]
    pk_mat@bl_EAD <- pk_mat@bl_EAD[-fmr]
    pk_mat@cr_EAD <- pk_mat@cr_EAD[-fmr]
  }

  # find columns of samples
  indS <- 1:ncol(pk_mat@EAD)
  indT <- sapply(suf,grep, x = colnames(pk_mat@EAD)) %>% unlist()
  if(length(indT) >0) indS <- indS[-indT]

  # neutral matrix (whithout shift)
  mat <- pk_mat@EAD[,indS]
  shift <- rep(shift,ncol(mat))
  if(overlay == TRUE) shift <- shift*(1:ncol(mat))
  mat <- t(t(mat) - shift) %>% as.data.frame()
  shift <- mean(shift)

  # mean
  if(treatment[1] == "mean"){
    pk_mat@EAD <- apply(mat,1,mean) %>% add(shift) %>% cbind(.,pk_mat@EAD)
    pk_mat@cr_EAD <- apply(pk_mat@cr_EAD[,indS],1,mean) %>% cbind(.,pk_mat@cr_EAD)

    SD <- pk_mat@cr_EAD$EAD_sd + 10^-9
    fmr <- pk_mat@cr_EAD[which(pk_mat@pk_windows == 0)[1],1]
    pk_mat@cr_EAD$EAD_snr_ave <- -(pk_mat@cr_EAD[,1]-fmr)/SD
  }
  # median
  if(treatment[1] == "median"){
    pk_mat@EAD <- apply(mat,1,median) %>% add(shift) %>% cbind(.,pk_mat@EAD)
    pk_mat@cr_EAD <- apply(pk_mat@cr_EAD[,indS],1,median) %>% cbind(.,pk_mat@cr_EAD)

    SD <- pk_mat@cr_EAD$EAD_sd + 10^-9
    fmr <- pk_mat@cr_EAD[which(pk_mat@pk_windows == 0)[1],1]
    pk_mat@cr_EAD$EAD_snr_med <- -(pk_mat@cr_EAD[,1]-fmr)/SD
  }

  # suffixe
  colnames(pk_mat@EAD)[1] <- paste0("raw_EAD",suffixe)
  colnames(pk_mat@cr_EAD)[1] <- paste0("EAD",suffixe)

  # recalcule baseline
  pk_mat@bl_EAD <- cbind(pk_mat@EAD[,1],pk_mat@bl_EAD)
  colnames(pk_mat@bl_EAD)[1] <- paste0("bl_EAD",suffixe)

  for (i in 1:ncol(pk_mat@peak)){ # i = 10
    xbm <- pk_mat@peak[2:3,i]
    ybm <- pk_mat@EAD[xbm,1]
    tbm <- pk_mat@time[xbm]

    fmr <- (ybm[2]-ybm[1])/(xbm[2]-xbm[1])
    fmr[2] <- ybm[1]-fmr[1]*xbm[1]
    abs <- 1:length(pk_mat@time)
    aff <- abs*fmr[1]+fmr[2]
    pk_mat@bl_EAD[xbm[1]:xbm[2],1] <- aff[xbm[1]:xbm[2]]
  }

  return(pk_mat)
}

# 06 Mesurement ####

#' return information on peaks, amplitude, position, time of depolarisation, SNR
#'
#' @param pk_mat
#'
#' @return a matrix
#' @export
#'
#' @examples
#' # depol_param <- pk.measure(list_ead)
pk.measure <- function(pk_mat = list_ead){ # pk_mat = list_test

  # check #
  if (class(pk_mat)[1] != "mead") stop("'pk_mat' must be a mead S4 object")

  # smart baseline
  bl <- pk_mat@cr_EAD*0 + pk_mat@param["cr_shift"]

  depol_param <- list()

  for (j in 1:ncol(pk_mat@cr_EAD)){ # j = 1 " names(pk_mat@cr_EAD)
    depol_param[[j]] <- matrix(NA,7,ncol(pk_mat@peak),
                               dimnames = list(c("ind_inf","ind_sup",
                                                 "Tinf","Tsup",
                                                 "Tdepol_min","Tdepol_sec",
                                                 "Idepol"),
                                               names(pk_mat@peak)))
    for (i in 1:ncol(pk_mat@peak)){ # i = 10

      # limit :
      xbm <- pk_mat@peak[2:3,i]
      # peak center :
      ybm <- pk_mat@FID[xbm]
      fmr <- (ybm[2]-ybm[1])/(xbm[2]-xbm[1])
      fmr[2] <- ybm[1]-fmr[1]*xbm[1]
      abs <- 1:length(pk_mat@time)
      aff <- abs*fmr[1]+fmr[2]

      fmr <- pk_mat@FID[xbm[1]:xbm[2]] - aff[xbm[1]:xbm[2]]
      appex <- xbm[1] + which.max(fmr) - 1

      # detect max inf
      xbm[1] <- xbm[1]+which.max(pk_mat@cr_EAD[xbm[1]:appex,j])-1

      # detect min sup
      xbm[2] <- appex+which.max(pk_mat@cr_EAD[appex:xbm[2],j])-1

      # calcul affine fonction
      ybm <- pk_mat@cr_EAD[xbm,j]
      tbm <- pk_mat@time[xbm]

      fmr <- (ybm[2]-ybm[1])/(xbm[2]-xbm[1])
      fmr[2] <- ybm[1]-fmr[1]*xbm[1]
      abs <- 1:length(pk_mat@time)
      aff <- abs*fmr[1]+fmr[2]
      bl[xbm[1]:xbm[2],j] <- pk_mat@cr_EAD[xbm[1]:xbm[2],j] - aff[xbm[1]:xbm[2]] + pk_mat@param["cr_shift"]

      fmr <- c(pk_mat@time[xbm[1]],pk_mat@time[xbm[2]])
      depol_param[[j]][,i] <- c(xbm,fmr,diff(fmr),diff(fmr)*60,min(bl[xbm[1]:xbm[2],j]))
    }
  }

  names(depol_param) <- names(pk_mat@cr_EAD) %>% str_remove_all("EAD_")
  if(length(grep("snr",names(depol_param))) > 0) depol_param <- depol_param[-grep("snr",names(depol_param))]

  if((length(which(names(depol_param) == "ave")))&(length(which(names(pk_mat@cr_EAD) == "EAD_snr_ave")))){
    dp <- rbind(depol_param[[which(names(depol_param) == "ave")]],NA)
    rownames(dp)[8] <- "SNR"
    for (i in 1:ncol(pk_mat@peak)) dp[8,i] <- max(pk_mat@cr_EAD$EAD_snr_ave[dp[1,i]:dp[2,i]])
    depol_param[[which(names(depol_param) == "ave")]] <- dp
  }

  if((length(which(names(depol_param) == "med")))&(length(which(names(pk_mat@cr_EAD) == "EAD_snr_med")))){
    dp <- rbind(depol_param[[which(names(depol_param) == "med")]],NA)
    rownames(dp)[8] <- "SNR"
    for (i in 1:ncol(pk_mat@peak)) dp[8,i] <- max(pk_mat@cr_EAD$EAD_snr_ave[dp[1,i]:dp[2,i]])
    depol_param[[which(names(depol_param) == "med")]] <- dp
  }

  depol_param <- lapply(depol_param,round,digits = 3)

  # df <- data.frame(abs = pk_mat@time ,
  #                  cr = pk_mat@cr_EAD[,1],
  #                  bl = bl[,1],
  #                  fid = pk_mat@FID*0.01,
  #                  wdw = pk_mat@pk_windows*0.01)
  # fig <- plot_ly(df, x = ~abs, y = ~cr, name = "cr", type = 'scatter', mode = 'lines')
  # fig <- fig %>% add_trace(y = ~bl, name = 'bl', mode = 'lines')
  # fig <- fig %>% add_trace(y = ~fid, name = 'fid', mode = 'lines')
  # fig <- fig %>% add_trace(y = ~wdw, name = 'wdw', mode = 'lines')
  # fig

  return(depol_param)
}

# 40 S4 section ####

#' object used for rEAD package
#'
#' @slot GC_EAD data.frame.
#' @slot peak data.frame.
#' @slot pk_res data.frame.
#' @slot delay numeric.
#' @slot file character.
#' @slot num numeric.
#' @slot wd character.
#' @slot type character.
#' @slot param numeric.
#'
#' @return a gcead object
#' @export
#'
#' @examples
#' #class(gcead)
setClass("gcead",representation(GC_EAD = "data.frame",
                                peak = "data.frame",
                                pk_res = "data.frame",
                                delay = "numeric",
                                file = "character",
                                num = "numeric",
                                wd = "character",
                                type = "character",
                                param = "numeric"))

#' object used for rEAD package. It's specific for multi EAD.
#'
#' @slot EAD numeric.
#' @slot time numeric.
#' @slot FID numeric.
#' @slot pk_windows numeric.
#' @slot bl_EAD numeric.
#' @slot cr_EAD numeric.
#' @slot peak data.frame.
#' @slot delay numeric.
#' @slot file character.
#' @slot num numeric.
#' @slot wd character.
#' @slot name character.
#' @slot param numeric.
#'
#' @return a mead object
#' @export
#'
#' @examples
#' #class(mead)
setClass("mead",representation(EAD = "data.frame",
                               time = "numeric",
                               FID = "numeric",
                               pk_windows = "numeric",
                               bl_EAD = "data.frame",
                               cr_EAD = "data.frame",
                               peak = "data.frame",
                               delay = "numeric",
                               file = "character",
                               num = "numeric",
                               wd = "character",
                               name = "character",
                               param = "numeric"))

# 98 others functions ####

#' brn.pk
#'
#' @param pki the peak number i
#' @param fdi the FID spectra
#' @param bm half windows maximum for the width of FID peak
#' @param sh half windows minimum for the width of FID peak
#'
#' @return a matriw with positions, limits and intensity of peaks
#' @noRd
brn.pk <- function(pki,fdi=fd, bm = bm_pk_FID, sh = sh_min_FID){ # pki = pk_det[100]

  nv <- fdi[(pki - bm):(pki-sh)] %>% scale()
  fmr <- which(nv > 0.1)
  b1 <- fmr[length(fmr)] %>% add(pki) %>% subtract(bm-sh)

  nv <- fdi[(pki+sh):(pki + bm)] %>% scale()
  b2 <- which(nv > 0.1)[1] %>% add(pki) %>% add(sh)

  return(c(pki,b1,b2))
}

#' fit.exact.value
#'
#' @param app the approx time
#' @param exact the vector with exact values
#'
#' @return the exact time
#' @noRd
fit.exact.value <- function(app = 11.12, exact = pk_mat@GC_EAD$time){
  ind_ex <- subtract(exact,app) %>% abs() %>% which.min()
  return(exact[ind_ex])
}


#' print.FID.peak
#'
#' @param pk_mat a gcead object
#'
#' @return a graph
#' @noRd
print.FID.peak <- function(pk_mat = gc_ead){
  title <- paste(pk_mat@file,"n#",pk_mat@num)

  dyG  <- dygraph(pk_mat@GC_EAD[,1:3], main = title) %>%
    dySeries("pk", label = "FID peak") %>%
    dyAxis("x", label = "Time and mass") %>%
    dyAxis("y", label = "Intensity (a.u.)") %>%
    dyOptions(axisLineColor = "navy", gridLineColor = "lightblue",
              useDataTimezone = FALSE) %>%
    dyRangeSelector() %>% dyUnzoom()

  tooltips <- list()
  for(i in 1:nrow(pk_mat@peak)){
    fmr <- list(list(x = pk_mat@peak[i,1], text = pk_mat@peak[i,1], series = "FID",
                     tooltip = pk_mat@peak[i,8], width=50, tickHeight = 0.1))
    tooltips <- c(tooltips, fmr)
  }

  annotator <- function(x,y){
    d = do.call(dyAnnotation,modifyList(list(dygraph=x),y))
    return(d)
  }

  dyG <- Reduce(annotator, tooltips, init=dyG )

  dyG$x$css = ".dygraphDefaultAnnotation {
                background-color: transparent;
                border: none;
                tickColor: transparent !important;
                color: black !important;
                width: initial !important;
                font-size: 60% !important;
              }"
  print(dyG)
}

# 99 obsolete ####


# END OF CODE ####
