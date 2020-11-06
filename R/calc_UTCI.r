#' Universal Climate Thermal Index (UTCI)
#' @description This function takes a dataframe object as input and computes UTCI.
#'
#' @param indat = dataframe containing the input variables
#'
#' @return A vector containing the daily UTCI
#' @examples
#' calc_UTCI(indat)
#'
#' Where the dataframe "indat" contains the following variables:
#'
#' tmin  = Temperature min    (degrees Celsius)
#' tmax  = Temperature max    (degrees Celsius)
#' rh =  Relative humidity  (%)
#' ws  = Wind speed at 10m  (m/s)
#' rs  = radiation          (MJ/m^2)
#' 
#' Also rquired is the latitude and longitude
#' lat = Latitude  (degrees)
#' lon = Longitude (degrees)
#' Elevation (needed) is retrieved from the PETr internal file "data/elev_dat.rda" which contains data from:
#' http://www.ecad.eu/download/ensembles/data/Grid_0.1deg_reg_ensemble/elev_ens_0.1deg_reg_v17.0e.nc
#' Missing values should be converted to NA
#'
#' @export

calc_UTCI <- function(indat) {
  dyn.load(system.file("libs/UTCIr.so", package="UTCIr"))
  data("elev_dat", package="UTCIr")
  tmax <- indat$tmax
  tmin <- indat$tmin
  rh <- indat$rh
  rs <- indat$rs
  ws <- indat$ws
  ndat <- length(tmax)
  
  lat <- indat$lat
  lon <- indat$lon
  # lat index
  lat_index <- which.min(abs(lat - elev_lat))
  # lon index
  lon_index <- which.min(abs(lon - elev_lon))
  
  elev <- elev_dat[lon_index, lat_index]
  if(is.na(elev))
  {
    utci <- array(dim=ndat, NA)
    return(utci)
  }
  # calc day count
  if(is.null(indat$jd))
  {
    jd <- as.POSIXlt(indat$tmax.dates)$yday + 1
  } else {
    jd <- indat$jd
  }
 
  # check whether all required input variables are present
  if (is.null(tmax) || is.null(tmin) || is.null(rh)|| is.null(rs) || is.null(ws)) {

    stop("Terminating - missing imput variables")
  }
  
  # need to change any NAs in input stream to -999.9
  tmax[is.na(tmax)] <- -999.9
  tmin[is.na(tmin)] <- -999.9
  rh[is.na(rh)] <- -999.9
  rs[is.na(rs)] <- -999.9
  ws[is.na(ws)] <- -999.9
  utci <- array(dim=ndat, -999.9)
  result <- .Fortran("calc_UTCI_vector",as.integer(jd), as.double(lat), as.double(tmax), as.double(tmin), as.double(rh), as.double(rs), as.double(ws),
                     as.integer(ndat), as.double(utci), PACKAGE = "UTCIr")
  utci <- result[[9]]
  
  # change -999.9 missing values back to NA
  utci[utci < -999.0] <- NA
  
  return(utci)
}


