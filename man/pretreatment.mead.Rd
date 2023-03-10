% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rEAD_fct.R
\name{pretreatment.mead}
\alias{pretreatment.mead}
\title{performs several pre-processings and peaks detection on mead S4 object}
\usage{
pretreatment.mead(
  pk_mat = list_ead,
  wsm = 125,
  wm_bl_FID = 400,
  ws_bl_FID = 10,
  th_pk_FID = 0.01,
  hws_pk_FID = 10,
  snr_pk_FID = 0.5
)
}
\arguments{
\item{pk_mat}{a mead object with spectra and metadata}

\item{wsm}{width smotthing for "fl" procma::savgol function}

\item{wm_bl_FID}{Width of local window for minimization/maximization for baseline::baseline.rollingBall function}

\item{ws_bl_FID}{Width of local window for smoothing for baseline::baseline.rollingBall function}

\item{th_pk_FID}{threshold of FID peak intensity}

\item{hws_pk_FID}{numeric, half window size for MALDIquant::detectPeaks}

\item{snr_pk_FID}{single numeric value for MALDIquant::detectPeaks. SNR is an abbreviation for signal-to-noise-ratio. A local maximum has to be higher than SNR*noise to be recognize as peak.}
}
\value{

}
\description{
performs several pre-processings and peaks detection on mead S4 object
}
\examples{
# list_test <- pretreatment.mead()
}
