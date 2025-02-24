% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_UTCI.r
\name{calc_UTCI}
\alias{calc_UTCI}
\title{Universal Climate Thermal Index (UTCI)}
\usage{
calc_UTCI(indat)
}
\arguments{
\item{indat}{= dataframe containing the input variables}
}
\value{
A vector containing the daily UTCI
}
\description{
This function takes a dataframe object as input and computes UTCI.
}
\examples{
calc_UTCI(indat)

Where the dataframe "indat" contains the following variables:

tmin  = Temperature min    (degrees Celsius)
tmax  = Temperature max    (degrees Celsius)
rh =  Relative humidity  (\%)
ws  = Wind speed at 10m  (m/s)
rs  = radiation          (MJ/m^2)

Also rquired is the latitude and longitude
lat = Latitude  (degrees)
lon = Longitude (degrees)
Elevation (needed) is retrieved from the UTCIr internal file "data/elev_dat.rda" which contains data from:
http://www.ecad.eu/download/ensembles/data/Grid_0.1deg_reg_ensemble/elev_ens_0.1deg_reg_v17.0e.nc
Missing values should be converted to NA

}
