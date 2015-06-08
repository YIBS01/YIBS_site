      module flux_com

      implicit none
      save

      integer ncid, start
      integer lonid, latid, timeid, vlonid, vlatid, vtimid
      integer gcid, ciid, coszid, nppid, gppid, gpp0id, roid, cflxid
      integer resprid, resplid, respwid, ctotid, csoilid
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      integer fo3id, dfo3id
#endif
#ifdef ACTIVE_GROWTH
      integer laipid, htpid
#endif
      integer laiid, fvid, o3id
      integer, parameter :: nrecmax = 500000
      integer nhour(nrecmax)

      end module flux_com


      subroutine open_output_netcdf(site, lon, lat, year)

      use flux_com
      implicit none
      include 'netcdf.inc'

      character*(*), intent(in) :: site
      real*8, intent(in)       :: lon
      real*8, intent(in)       :: lat
      integer, intent(in)      :: year
      integer nstat
      integer lenc
      character*120 fsite
      character*8 cyear

      fsite = site(1:lenc(site))//'_flux.nc'
      write(cyear,'(I4.4)') year
      start = 1

      nstat = NF_CREATE(fsite(1:lenc(fsite)),NF_CLOBBER,ncid)
      IF(nstat .NE. NF_NOERR) Call handle_err(nstat)

CCCCCC   DEFINE DIMENSIONS    CCCCCCCC

      nstat = NF_DEF_DIM(ncid, 'lon', 1, lonid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_DIM(ncid, 'lat', 1, latid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_DIM(ncid, 'time',NF_UNLIMITED, timeid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'title',26,
     &       'Standalone vegetation flux')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'lon',NF_DOUBLE,1,lonid,vlonid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_VAR(ncid,'lat',NF_DOUBLE,1,latid,vlatid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_VAR(ncid,'time',NF_INT,1,timeid,vtimid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_PUT_ATT_TEXT(ncid,vlonid,'long_name',9,
     &       'Longitude')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid,vlonid,'units',9,
     &       'degrees_e')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid,vlatid,'long_name',8,
     &       'Latitude')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid,vlatid,'units',9,
     &       'degrees_n')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid,vtimid,'long_name',4,
     &       'Time')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid,vtimid,'units',30,
     &       'hours since '//cyear(1:lenc(cyear))//'-1-1 00:00:0.0')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid,vtimid,'calendar',7,
     &       '365_day')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'gc',NF_FLOAT,1,timeid,gcid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, gcid,'long_name',33,
     &       'Canopy conductance of water vapor')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, gcid,'units',5,
     &       'm s-1')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'ci',NF_FLOAT,1,timeid,ciid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, ciid,'long_name',34,
     &       'Internal foliage CO2 concentration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, ciid,'units',7,
     &       'mol m-3')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'cosz',NF_FLOAT,1,timeid,coszid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, coszid,'long_name',28,
     &       'cosine of solar zenith angle')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, coszid,'units',4,
     &       'none')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'NPP',NF_FLOAT,1,timeid,nppid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, nppid,'long_name',24,
     &       'Net primary productivity')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, nppid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'GPP',NF_FLOAT,1,timeid,gppid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, gppid,'long_name',26,
     &       'Gross primary productivity')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, gppid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'GPP0',NF_FLOAT,1,timeid,gpp0id)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, gpp0id,'long_name',34,
     &       'Gross primary productivity offline')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, gpp0id,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'R_auto',NF_FLOAT,1,timeid,roid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, roid,'long_name',23,
     &       'Autotrophic respiration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, roid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'CO2flux',NF_FLOAT,1,timeid,cflxid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, cflxid,'long_name',15,
     &       'Net CO2 flux up')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, cflxid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'resp_r',NF_FLOAT,1,timeid, resprid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, resprid,'long_name',16,
     &       'Root respiration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, resprid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'resp_l',NF_FLOAT,1,timeid, resplid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, resplid,'long_name',16,
     &       'Leaf respiration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, resplid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'resp_w',NF_FLOAT,1,timeid, respwid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, respwid,'long_name',16,
     &       'Wood respiration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, respwid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'Carbon_total',NF_FLOAT,1,timeid, ctotid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, ctotid,'long_name',25,
     &       'Total Land Carbon Storage')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, ctotid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'Carbon_soil',NF_FLOAT,1,timeid, csoilid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, csoilid,'long_name',25,
     &       'Total Soil Carbon Storage')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, csoilid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)

      nstat = NF_DEF_VAR(ncid,'FO3',NF_FLOAT,1,timeid,fo3id)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, fo3id,'long_name',21,
     &       'ozone flux to stomata')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, fo3id,'units',15,
     &       'nmol O3 m-2 s-1')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'dFO3',NF_FLOAT,1,timeid,dfo3id)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, dfo3id,'long_name',28,
     &       'excess ozone flux to stomata')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, dfo3id,'units',15,
     &       'nmol O3 m-2 s-1')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#endif

#ifdef ACTIVE_GROWTH

      nstat = NF_DEF_VAR(ncid,'lai_p',NF_FLOAT,1,timeid,laipid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, laipid,'long_name',14,
     &       'prognostic LAI')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, laipid,'units',6,
     &       'm2 m-2')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'ht_p',NF_FLOAT,1,timeid,htpid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, htpid,'long_name',17,
     &       'prognostic height')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, htpid,'units',1,
     &       'm')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#endif

      nstat = NF_DEF_VAR(ncid,'O3conc',NF_FLOAT,1,timeid,o3id)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, o3id,'long_name',16,
     &       'O3 concentration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, o3id,'units',3,
     &       'ppb')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'lai',NF_FLOAT,1,timeid,laiid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, laiid,'long_name',15,
     &       'Leaf area index')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, laiid,'units',6,
     &       'm2 m-2')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncid,'vfrac',NF_FLOAT,1,timeid,fvid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, fvid,'long_name',18,
     &       'Total veg fraction')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncid, fvid,'units',8,
     &       'fraction')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_ENDDEF(ncid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_PUT_VAR_DOUBLE(ncid,vlonid,lon)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VAR_DOUBLE(ncid,vlatid,lat)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      end subroutine open_output_netcdf


      subroutine close_output_netcdf

      use flux_com,   only: ncid
      implicit none
      include 'netcdf.inc'

      integer nstat

      nstat=NF_CLOSE(ncid)
      IF(nstat .NE. NF_NOERR) Call handle_err(nstat)

      end subroutine close_output_netcdf


      subroutine write_output_netcdf(GCANOPY,Ci,NPP,GPP,GPP0,CO2flux,
     &       resp_r, resp_l, resp_w, c_tot, c_soil,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &       FO3, dFO3,
#endif
#ifdef ACTIVE_GROWTH
     &       lai_p, ht_p,
#endif
     &       R_auto, coszen, o3s, lai0, fv0, ny, nm, nd, nh, year)

      use flux_com
      implicit none
      include 'netcdf.inc'

      real*8, intent(in)   :: GCANOPY, Ci, NPP, GPP, GPP0, CO2flux
      real*8, intent(in)   :: resp_r, resp_l, resp_w
      real*8, intent(in)   :: c_tot,  c_soil          
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8, intent(in)   :: FO3, dFO3
#endif
#ifdef ACTIVE_GROWTH
      real*8, intent(in)   :: lai_p,  ht_p
#endif
      real*8, intent(in)   :: R_auto, coszen, o3s, lai0, fv0
      integer, intent(in)  :: ny, nm, nd, nh, year

      integer, parameter :: mdays(12) = (/
     &    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      integer totdays, n, nstat

      If (start .ge. nrecmax) then
         print*, 'Record length exceed limit, increase nrecmax!'
         stop
      Endif

      totdays = (ny-year)*365
      do n = 1, nm
         totdays = totdays + mdays(n)
      enddo
      totdays = totdays - mdays(nm) + nd
      nhour(start) = (totdays-1)*24 + nh

      nstat = NF_PUT_VARA_REAL(ncid,gcid,start,1,real(gcanopy))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,ciid,start,1,real(ci))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,nppid,start,1,real(NPP))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,gppid,start,1,real(GPP))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,gpp0id,start,1,real(GPP0))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,cflxid,start,1,real(CO2flux))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,resprid,start,1,real(resp_r))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,resplid,start,1,real(resp_l))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,respwid,start,1,real(resp_w))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,ctotid,start,1,real(c_tot))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,csoilid,start,1,real(c_soil))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      nstat = NF_PUT_VARA_REAL(ncid,fo3id,start,1,real(FO3))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,dfo3id,start,1,real(dFO3))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#endif
#ifdef ACTIVE_GROWTH
      nstat = NF_PUT_VARA_REAL(ncid,laipid,start,1,real(lai_p))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,htpid,start,1,real(ht_p))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#endif
      nstat = NF_PUT_VARA_REAL(ncid,o3id,start,1,real(o3s))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,roid,start,1,real(R_auto))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,coszid,start,1,real(coszen))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,laiid,start,1,real(lai0))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncid,fvid,start,1,real(fv0))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VAR_INT(ncid,vtimid,nhour(1:start))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      start = start + 1

      end subroutine write_output_netcdf




