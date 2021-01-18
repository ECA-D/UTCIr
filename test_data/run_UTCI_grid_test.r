library(gridclimind)
library(ncdf4)

in.dir <- "input/"
author.data <- list(ERA ="5")
out.dir <- "UTCI_output/"

############ UTCI ############

UTCI.file <- "UTCI_testfile.nc"
out.file <- sprintf("%s%s", out.dir, UTCI.file)
if(file.exists(out.file)){
  file.remove(out.file)
}
input.files <- c(paste0(in.dir,"tn_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"tx_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"rh_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"rs_0.25deg_regular_1979-2018_sub.nc"),
                 paste0(in.dir,"ws_0.25deg_regular_1979-2018_sub.nc"))

create.utci.from.files(input.files, out.file, input.files[1], author.data,  parallel=FALSE)


##### test output #######

in.nc_1= nc_open("UTCI_output/UTCI_testfile_master.nc")
lats <- ncvar_get( in.nc_1, "latitude")   # coordinate variable
nlat <- length(lats)
lons <- ncvar_get( in.nc_1, "longitude")   # coordinate variable
nlon <- length(lons)
tm <- ncvar_get( in.nc_1, "time")/24
nt <- length(tm)
indat1 <- as.matrix(ncvar_get( in.nc_1, "utci", start=c(1,1,1), count=c(nlon,nlat,nt)) )

in.nc_2= nc_open("UTCI_output/UTCI_testfile.nc")
indat2 <- as.matrix(ncvar_get( in.nc_2, "utci", start=c(1,1,1), count=c(nlon,nlat,nt)) )

if(all.equal.numeric(indat1, indat2, tolerance=0.0001)==TRUE)
{
  print("UTCI output matches master file")
}else{
  print("Warning: UTCI output differs from master file")
}
nc_close(in.nc_1)
nc_close(in.nc_2)
