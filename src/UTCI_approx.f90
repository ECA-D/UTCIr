

subroutine calc_UTCI_vector(jd, lat, tmax, tmin, rh, rs, ws, ndat, utci)
!   this index requires several input vectors with daily resolution:
!   jd          yearday count 1-365/366
!   Tmax data     degree Celsius
!   Tmin data     degree Celsius

!   RH data     RH                                  units ????????
!   RS data     radiation - downwelling short wave  units ????????
!   WS data     wind 10m m/s
!   lat         latitude

    implicit none
    real*8  lat, tmax(1:ndat), tmin(1:ndat), rh(1:ndat), rs(1:ndat), ws(1:ndat)
    real*8  alpha, ehPa, Tmrt, utci(1:ndat), UTCI_approx
    integer jd(1:ndat), ndat
    integer i
    alpha = 0.85

 
    do i=1,ndat
        ! check for missing values
        if((tmax(i).lt.-999.0).or.(tmin(i).lt.-999.0).or.(rh(i).lt.-999.0).or.(rs(i).lt.-999.0).or.(ws(i).lt.-999.0)) then
            utci(i) = -999.9
        else

            call mrt(jd(i), lat, tmax(i), alpha, rs(i), Tmrt)
            call vapourpressure(rh(i), tmax(i), tmin(i), ehPa)

        ! vapour pressure is required in hPa (given in kPa)
            ehPa = 10.0*ehPa

            utci(i) = UTCI_approx(tmax(i), ehPa, Tmrt, ws(i))
        endif
    enddo

return
end



   !~ UTCI, Version a 0.002, October 2009
    !~ Copyright (C) 2009  Peter Broede
    
    !~ Program for calculating UTCI Temperature (UTCI)
    !~ released for public use after termination of COST Action 730
    
    !~ replaces Version a 0.001, from September 2009

!~ **********************************************
      DOUBLE precision function UTCI_approx(Ta,ehPa,Tmrt,va)
!~ **********************************************
 !~ DOUBLE PRECISION Function value is the UTCI in degree Celsius
 !~ computed by a 6th order approximating polynomial from the 4 Input paramters 
 !~ 
 !~ Input parameters (all of type DOUBLE PRECISION)
 !~ - Ta       : air temperature, degree Celsius
 !~ - ehPa    : water vapour presure, hPa=hecto Pascal
 !~ - Tmrt   : mean radiant temperature, degree Celsius
 !~ - va10m  : wind speed 10 m above ground level in m/s
 !~ 
 !~  UTCI_approx, Version a 0.002, October 2009
 !~  Copyright (C) 2009  Peter Broede

      implicit none
        !~ type of input of the argument list
       DOUBLE PRECISION Ta,va,Tmrt,ehPa,Pa,D_Tmrt;
          D_TMRT=Tmrt-Ta
       	  PA = ehPa/10.0; !~ use vapour pressure in kPa
        !~ calculate 6th order polynomial as approximation     
      UTCI_approx=Ta+&
		( 6.07562052D-01 )   + &
		( -2.27712343D-02 ) * Ta + &
		( 8.06470249D-04 ) * Ta*Ta + &
		( -1.54271372D-04 ) * Ta*Ta*Ta + &
		( -3.24651735D-06 ) * Ta*Ta*Ta*Ta + &
		( 7.32602852D-08 ) * Ta*Ta*Ta*Ta*Ta + &
		( 1.35959073D-09 ) * Ta*Ta*Ta*Ta*Ta*Ta + &
		( -2.25836520D+00 ) * va + &
		( 8.80326035D-02 ) * Ta*va + &
		( 2.16844454D-03 ) * Ta*Ta*va + &
		( -1.53347087D-05 ) * Ta*Ta*Ta*va + &
		( -5.72983704D-07 ) * Ta*Ta*Ta*Ta*va + &
		( -2.55090145D-09 ) * Ta*Ta*Ta*Ta*Ta*va + &
		( -7.51269505D-01 ) * va*va + &
		( -4.08350271D-03 ) * Ta*va*va + &
		( -5.21670675D-05 ) * Ta*Ta*va*va + &
		( 1.94544667D-06 ) * Ta*Ta*Ta*va*va + &
		( 1.14099531D-08 ) * Ta*Ta*Ta*Ta*va*va + &
		( 1.58137256D-01 ) * va*va*va + &
		( -6.57263143D-05 ) * Ta*va*va*va + &
		( 2.22697524D-07 ) * Ta*Ta*va*va*va + &
		( -4.16117031D-08 ) * Ta*Ta*Ta*va*va*va + &
		( -1.27762753D-02 ) * va*va*va*va + &
		( 9.66891875D-06 ) * Ta*va*va*va*va + &
		( 2.52785852D-09 ) * Ta*Ta*va*va*va*va + &
		( 4.56306672D-04 ) * va*va*va*va*va + &
		( -1.74202546D-07 ) * Ta*va*va*va*va*va + &
		( -5.91491269D-06 ) * va*va*va*va*va*va + &
		( 3.98374029D-01 ) * D_Tmrt + &
		( 1.83945314D-04 ) * Ta*D_Tmrt + &
		( -1.73754510D-04 ) * Ta*Ta*D_Tmrt + &
		( -7.60781159D-07 ) * Ta*Ta*Ta*D_Tmrt + &
		( 3.77830287D-08 ) * Ta*Ta*Ta*Ta*D_Tmrt + &
		( 5.43079673D-10 ) * Ta*Ta*Ta*Ta*Ta*D_Tmrt + &
		( -2.00518269D-02 ) * va*D_Tmrt + &
		( 8.92859837D-04 ) * Ta*va*D_Tmrt + &
		( 3.45433048D-06 ) * Ta*Ta*va*D_Tmrt + &
		( -3.77925774D-07 ) * Ta*Ta*Ta*va*D_Tmrt + &
		( -1.69699377D-09 ) * Ta*Ta*Ta*Ta*va*D_Tmrt + &
		( 1.69992415D-04 ) * va*va*D_Tmrt + &
		( -4.99204314D-05 ) * Ta*va*va*D_Tmrt + &
		( 2.47417178D-07 ) * Ta*Ta*va*va*D_Tmrt + &
		( 1.07596466D-08 ) * Ta*Ta*Ta*va*va*D_Tmrt + &
		( 8.49242932D-05 ) * va*va*va*D_Tmrt + &
		( 1.35191328D-06 ) * Ta*va*va*va*D_Tmrt + &
		( -6.21531254D-09 ) * Ta*Ta*va*va*va*D_Tmrt + &
		( -4.99410301D-06 ) * va*va*va*va*D_Tmrt + &
		( -1.89489258D-08 ) * Ta*va*va*va*va*D_Tmrt + &
		( 8.15300114D-08 ) * va*va*va*va*va*D_Tmrt + &
		( 7.55043090D-04 ) * D_Tmrt*D_Tmrt + &
		( -5.65095215D-05 ) * Ta*D_Tmrt*D_Tmrt + &
		( -4.52166564D-07 ) * Ta*Ta*D_Tmrt*D_Tmrt + &
		( 2.46688878D-08 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt + &
		( 2.42674348D-10 ) * Ta*Ta*Ta*Ta*D_Tmrt*D_Tmrt + &
		( 1.54547250D-04 ) * va*D_Tmrt*D_Tmrt + &
		( 5.24110970D-06 ) * Ta*va*D_Tmrt*D_Tmrt + &
		( -8.75874982D-08 ) * Ta*Ta*va*D_Tmrt*D_Tmrt + &
		( -1.50743064D-09 ) * Ta*Ta*Ta*va*D_Tmrt*D_Tmrt + &
		( -1.56236307D-05 ) * va*va*D_Tmrt*D_Tmrt + &
		( -1.33895614D-07 ) * Ta*va*va*D_Tmrt*D_Tmrt + &
		( 2.49709824D-09 ) * Ta*Ta*va*va*D_Tmrt*D_Tmrt + &
		( 6.51711721D-07 ) * va*va*va*D_Tmrt*D_Tmrt + &
		( 1.94960053D-09 ) * Ta*va*va*va*D_Tmrt*D_Tmrt + &
		( -1.00361113D-08 ) * va*va*va*va*D_Tmrt*D_Tmrt + &
		( -1.21206673D-05 ) * D_Tmrt*D_Tmrt*D_Tmrt + &
		( -2.18203660D-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 7.51269482D-09 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 9.79063848D-11 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 1.25006734D-06 ) * va*D_Tmrt*D_Tmrt*D_Tmrt + &
		( -1.81584736D-09 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt + &
		( -3.52197671D-10 ) * Ta*Ta*va*D_Tmrt*D_Tmrt*D_Tmrt + &
		( -3.36514630D-08 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 1.35908359D-10 ) * Ta*va*va*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 4.17032620D-10 ) * va*va*va*D_Tmrt*D_Tmrt*D_Tmrt + &
		( -1.30369025D-09 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 4.13908461D-10 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 9.22652254D-12 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( -5.08220384D-09 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( -2.24730961D-11 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 1.17139133D-10 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 6.62154879D-10 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 4.03863260D-13 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 1.95087203D-12 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( -4.73602469D-12 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt + &
		( 5.12733497D+00 ) * Pa + &
		( -3.12788561D-01 ) * Ta*Pa + &
		( -1.96701861D-02 ) * Ta*Ta*Pa + &
		( 9.99690870D-04 ) * Ta*Ta*Ta*Pa + &
		( 9.51738512D-06 ) * Ta*Ta*Ta*Ta*Pa + &
		( -4.66426341D-07 ) * Ta*Ta*Ta*Ta*Ta*Pa + &
		( 5.48050612D-01 ) * va*Pa + &
		( -3.30552823D-03 ) * Ta*va*Pa + &
		( -1.64119440D-03 ) * Ta*Ta*va*Pa + &
		( -5.16670694D-06 ) * Ta*Ta*Ta*va*Pa + &
		( 9.52692432D-07 ) * Ta*Ta*Ta*Ta*va*Pa + &
		( -4.29223622D-02 ) * va*va*Pa + &
		( 5.00845667D-03 ) * Ta*va*va*Pa + &
		( 1.00601257D-06 ) * Ta*Ta*va*va*Pa + &
		( -1.81748644D-06 ) * Ta*Ta*Ta*va*va*Pa + &
		( -1.25813502D-03 ) * va*va*va*Pa + &
		( -1.79330391D-04 ) * Ta*va*va*va*Pa + &
		( 2.34994441D-06 ) * Ta*Ta*va*va*va*Pa + &
		( 1.29735808D-04 ) * va*va*va*va*Pa + &
		( 1.29064870D-06 ) * Ta*va*va*va*va*Pa + &
		( -2.28558686D-06 ) * va*va*va*va*va*Pa + &
		( -3.69476348D-02 ) * D_Tmrt*Pa + &
		( 1.62325322D-03 ) * Ta*D_Tmrt*Pa + &
		( -3.14279680D-05 ) * Ta*Ta*D_Tmrt*Pa + &
		( 2.59835559D-06 ) * Ta*Ta*Ta*D_Tmrt*Pa + &
		( -4.77136523D-08 ) * Ta*Ta*Ta*Ta*D_Tmrt*Pa + &
		( 8.64203390D-03 ) * va*D_Tmrt*Pa + &
		( -6.87405181D-04 ) * Ta*va*D_Tmrt*Pa + &
		( -9.13863872D-06 ) * Ta*Ta*va*D_Tmrt*Pa + &
		( 5.15916806D-07 ) * Ta*Ta*Ta*va*D_Tmrt*Pa + &
		( -3.59217476D-05 ) * va*va*D_Tmrt*Pa + &
		( 3.28696511D-05 ) * Ta*va*va*D_Tmrt*Pa + &
		( -7.10542454D-07 ) * Ta*Ta*va*va*D_Tmrt*Pa + &
		( -1.24382300D-05 ) * va*va*va*D_Tmrt*Pa + &
		( -7.38584400D-09 ) * Ta*va*va*va*D_Tmrt*Pa + &
		( 2.20609296D-07 ) * va*va*va*va*D_Tmrt*Pa + &
		( -7.32469180D-04 ) * D_Tmrt*D_Tmrt*Pa + &
		( -1.87381964D-05 ) * Ta*D_Tmrt*D_Tmrt*Pa + &
		( 4.80925239D-06 ) * Ta*Ta*D_Tmrt*D_Tmrt*Pa + &
		( -8.75492040D-08 ) * Ta*Ta*Ta*D_Tmrt*D_Tmrt*Pa + &
		( 2.77862930D-05 ) * va*D_Tmrt*D_Tmrt*Pa + &
		( -5.06004592D-06 ) * Ta*va*D_Tmrt*D_Tmrt*Pa + &
		( 1.14325367D-07 ) * Ta*Ta*va*D_Tmrt*D_Tmrt*Pa + &
		( 2.53016723D-06 ) * va*va*D_Tmrt*D_Tmrt*Pa + &
		( -1.72857035D-08 ) * Ta*va*va*D_Tmrt*D_Tmrt*Pa + &
		( -3.95079398D-08 ) * va*va*va*D_Tmrt*D_Tmrt*Pa + &
		( -3.59413173D-07 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( 7.04388046D-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( -1.89309167D-08 ) * Ta*Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( -4.79768731D-07 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( 7.96079978D-09 ) * Ta*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( 1.62897058D-09 ) * va*va*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( 3.94367674D-08 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( -1.18566247D-09 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( 3.34678041D-10 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( -1.15606447D-10 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa + &
		( -2.80626406D+00 ) * Pa*Pa + &
		( 5.48712484D-01 ) * Ta*Pa*Pa + &
		( -3.99428410D-03 ) * Ta*Ta*Pa*Pa + &
		( -9.54009191D-04 ) * Ta*Ta*Ta*Pa*Pa + &
		( 1.93090978D-05 ) * Ta*Ta*Ta*Ta*Pa*Pa + &
		( -3.08806365D-01 ) * va*Pa*Pa + &
		( 1.16952364D-02 ) * Ta*va*Pa*Pa + &
		( 4.95271903D-04 ) * Ta*Ta*va*Pa*Pa + &
		( -1.90710882D-05 ) * Ta*Ta*Ta*va*Pa*Pa + &
		( 2.10787756D-03 ) * va*va*Pa*Pa + &
		( -6.98445738D-04 ) * Ta*va*va*Pa*Pa + &
		( 2.30109073D-05 ) * Ta*Ta*va*va*Pa*Pa + &
		( 4.17856590D-04 ) * va*va*va*Pa*Pa + &
		( -1.27043871D-05 ) * Ta*va*va*va*Pa*Pa + &
		( -3.04620472D-06 ) * va*va*va*va*Pa*Pa + &
		( 5.14507424D-02 ) * D_Tmrt*Pa*Pa + &
		( -4.32510997D-03 ) * Ta*D_Tmrt*Pa*Pa + &
		( 8.99281156D-05 ) * Ta*Ta*D_Tmrt*Pa*Pa + &
		( -7.14663943D-07 ) * Ta*Ta*Ta*D_Tmrt*Pa*Pa + &
		( -2.66016305D-04 ) * va*D_Tmrt*Pa*Pa + &
		( 2.63789586D-04 ) * Ta*va*D_Tmrt*Pa*Pa + &
		( -7.01199003D-06 ) * Ta*Ta*va*D_Tmrt*Pa*Pa + &
		( -1.06823306D-04 ) * va*va*D_Tmrt*Pa*Pa + &
		( 3.61341136D-06 ) * Ta*va*va*D_Tmrt*Pa*Pa + &
		( 2.29748967D-07 ) * va*va*va*D_Tmrt*Pa*Pa + &
		( 3.04788893D-04 ) * D_Tmrt*D_Tmrt*Pa*Pa + &
		( -6.42070836D-05 ) * Ta*D_Tmrt*D_Tmrt*Pa*Pa + &
		( 1.16257971D-06 ) * Ta*Ta*D_Tmrt*D_Tmrt*Pa*Pa + &
		( 7.68023384D-06 ) * va*D_Tmrt*D_Tmrt*Pa*Pa + &
		( -5.47446896D-07 ) * Ta*va*D_Tmrt*D_Tmrt*Pa*Pa + &
		( -3.59937910D-08 ) * va*va*D_Tmrt*D_Tmrt*Pa*Pa + &
		( -4.36497725D-06 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
		( 1.68737969D-07 ) * Ta*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
		( 2.67489271D-08 ) * va*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
		( 3.23926897D-09 ) * D_Tmrt*D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa + &
		( -3.53874123D-02 ) * Pa*Pa*Pa + &
		( -2.21201190D-01 ) * Ta*Pa*Pa*Pa + &
		( 1.55126038D-02 ) * Ta*Ta*Pa*Pa*Pa + &
		( -2.63917279D-04 ) * Ta*Ta*Ta*Pa*Pa*Pa + &
		( 4.53433455D-02 ) * va*Pa*Pa*Pa + &
		( -4.32943862D-03 ) * Ta*va*Pa*Pa*Pa + &
		( 1.45389826D-04 ) * Ta*Ta*va*Pa*Pa*Pa + &
		( 2.17508610D-04 ) * va*va*Pa*Pa*Pa + &
		( -6.66724702D-05 ) * Ta*va*va*Pa*Pa*Pa + &
		( 3.33217140D-05 ) * va*va*va*Pa*Pa*Pa + &
		( -2.26921615D-03 ) * D_Tmrt*Pa*Pa*Pa + &
		( 3.80261982D-04 ) * Ta*D_Tmrt*Pa*Pa*Pa + &
		( -5.45314314D-09 ) * Ta*Ta*D_Tmrt*Pa*Pa*Pa + &
		( -7.96355448D-04 ) * va*D_Tmrt*Pa*Pa*Pa + &
		( 2.53458034D-05 ) * Ta*va*D_Tmrt*Pa*Pa*Pa + &
		( -6.31223658D-06 ) * va*va*D_Tmrt*Pa*Pa*Pa + &
		( 3.02122035D-04 ) * D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
		( -4.77403547D-06 ) * Ta*D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
		( 1.73825715D-06 ) * va*D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
		( -4.09087898D-07 ) * D_Tmrt*D_Tmrt*D_Tmrt*Pa*Pa*Pa + &
		( 6.14155345D-01 ) * Pa*Pa*Pa*Pa + &
		( -6.16755931D-02 ) * Ta*Pa*Pa*Pa*Pa + &
		( 1.33374846D-03 ) * Ta*Ta*Pa*Pa*Pa*Pa + &
		( 3.55375387D-03 ) * va*Pa*Pa*Pa*Pa + &
		( -5.13027851D-04 ) * Ta*va*Pa*Pa*Pa*Pa + &
		( 1.02449757D-04 ) * va*va*Pa*Pa*Pa*Pa + &
		( -1.48526421D-03 ) * D_Tmrt*Pa*Pa*Pa*Pa + &
		( -4.11469183D-05 ) * Ta*D_Tmrt*Pa*Pa*Pa*Pa + &
		( -6.80434415D-06 ) * va*D_Tmrt*Pa*Pa*Pa*Pa + &
		( -9.77675906D-06 ) * D_Tmrt*D_Tmrt*Pa*Pa*Pa*Pa + &
		( 8.82773108D-02 ) * Pa*Pa*Pa*Pa*Pa + &
		( -3.01859306D-03 ) * Ta*Pa*Pa*Pa*Pa*Pa + &
		( 1.04452989D-03 ) * va*Pa*Pa*Pa*Pa*Pa + &
		( 2.47090539D-04 ) * D_Tmrt*Pa*Pa*Pa*Pa*Pa + &
		( 1.48348065D-03 ) * Pa*Pa*Pa*Pa*Pa*Pa 
      return
      END


!~ **********************************************
      real FUNCTION es(ta)
!~ **********************************************
!~ calculates saturation vapour pressure over water in hPa for input air temperature (ta) in celsius according to:
!~ Hardy, R.; ITS-90 Formulations for Vapor Pressure, Frostpoint Temperature, Dewpoint Temperature and Enhancement Factors in the Range -100 to 100 °C; 
!~ Proceedings of Third International Symposium on Humidity and Moisture; edited by National Physical Laboratory (NPL), London, 1998, pp. 214-221
!~ http://www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf (retrieved 2008-10-01)
      
      implicit none
!
      real ta, tk
      INTEGER I
      REAL :: g(0:7)=(/&
				-2.8365744E3,&
				-6.028076559E3,&
				1.954263612E1,&
				-2.737830188E-2,&
				1.6261698E-5,&
				7.0229056E-10,&
				-1.8680009E-13,&
				2.7150305 /)
!      
      tk=ta+273.15 		! air temp in K
      es=g(7)*log(tk)
      do i=0,6
        es=es+g(i)*tk**(i-2)  
      end do
      es=exp(es)*0.01	! *0.01: convert Pa to hPa
!
      return
      END

!############### end of Broede routines ########################################

subroutine mrt(yearday,phi,Temp,alpha,qir,Tmrt)
! documentation to this subroutine is in rep110610
! or in the UTCI documentation & P.O. Fanger (1970)
!
! INPUT
! yearday: day in the year [1-365/366]
! phi:     latitude of the station [deg]
! Tumrt:   unirradiated mean radiant temperature [deg C]
! alpha:   absorptance of the outer surface of the person [1]
! qir:     irradiation which strikes the person [MJ m^-2 d^-1]
! Temp:    Tmax (tmax)

! OUTPUT
! Tmrt:    mean radiant temperature [deg C]
      implicit none

      integer      yearday
      real*8       xgamma,arg,phi,Tumrt,alpha,qir,Tmrt,fp
      real*8       C2K,K2C,xdum,temp

      Tumrt = C2K(Temp)

! note: the constant here is different than used by Fanger.
! his units are kcal m^-2 hr^-1, here we adhere to the units used
! in the PET routines: MJ m^-2 d^-1. The factor between these units is: 9.968
! which (approximately) explains the difference.
      arg=Tumrt**4 + 0.210*(10**9)*fp(yearday,phi)*alpha*qir

      Tmrt=sqrt(sqrt(arg))

      Tmrt=K2C(Tmrt)

!     write(6,*) 'Tumrt: ',K2C(Tumrt)
!     write(6,*) 'fp*alpha*qir: ',xdum*9.968
!     write(6,*) 'Tmrt: ',Tmrt

return
end
!
!------------------------------------------------------------
double precision function fp(yearday,phi)
! documentation to this subroutine is in rep110610

! input
!     yearday: day of the year [1]
!     phi: latitude of the station
! output
!     fp: projected area factor [1]

      implicit none

      integer      yearday
      real*8       xgamma,arg,phi,gamma_sun_elev,pi
      pi = 3.14158255558979

      xgamma = gamma_sun_elev(yearday,phi)

      arg = xgamma*(1.0d0 - xgamma*xgamma/48402.0d0)
      fp = 0.308d0*cos(2*pi*arg/360.0)

return
end

!-------------------------------------------------------------------
!
double precision function gamma_sun_elev(yearday,phi)
! documentation to this subroutine is in rep110610
! the assumption is made (see rep110610) that the hour angle = 0
!
! input
!     yearday: day of the year [1]
!     phi: latitude of the station [deg]
! output
!     gamma: elevation of the sun [rad]
      implicit none

      integer      yearday
      real*8       delta,phi,singamma,pi,latr
      pi = 3.14158255558979

      latr = phi*(pi/180.0)

! delta = sun declination
      delta=-23.44*cos((pi/180.0)*(yearday+10)*360.0/365.0)

      singamma=cos(2*pi*delta/360.0)*cos(latr) + sin(2*pi*delta/360.0)*sin(latr)

      gamma_sun_elev=acos(singamma)

! convert from radians -> degrees
      gamma_sun_elev = 180.0*gamma_sun_elev/pi

return
end
!
!------------------------------------------------------------
!
subroutine vapourpressure(rh,tmax,tmin,vap)
! calculates vapour pressure, following Allen et al 1994 
! eq. 1.16
! vapourpressure vap is in [kPa]
!
      implicit none

      real*8           rh,tmax,tmin,fac,vap, sattmin, sattmax

      call satvapourpressure(tmax, sattmax)
      call satvapourpressure(tmin, sattmin)

      fac = (50.0/sattmax + 50.0/sattmin)
      vap = rh/fac

return
end
!
!-------------------------------------------------------------------
!
subroutine satvapourpressure(tm,satvap)
! calculates saturation vapour pressure, following Tetens (1930)
! see document of Hess, Cranfield University, form. 2
! see also  Allen et al 1994 eq. 1.10
!
! tm is in [dC]
! satvap is in [kPa]
      implicit none

      real*8           tm,fac,satvap

      fac = 17.27*tm/(tm + 273.3)
!      satvap = 0.61078d0*exp(fac)
      satvap = (2504.0/4099.0)*exp(fac)

return
end
!
!-------------------------------------------------------------------
!
double precision function C2K(temp)
! converts Celsius to Kelvin
      implicit none
!
      real*8       temp

      C2K = temp + 273.15d0

return
end
!
!------------------------------------------------------------
!
double precision function K2C(temp)
! converts Kelvin to Celsius
      implicit none
!
      real*8       temp

      K2C = temp - 273.15d0

return
end
!
!------------------------------------------------------------
!

