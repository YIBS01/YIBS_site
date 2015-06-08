      module site_com

      implicit none
      save

      real*4 t1s, t1c, q1s, p1s, w1s, pr1s, pf1s, co21, o31, ch1, fw1
      real*4 st1(6), sm1(6), si1(6), sp1(6)
      integer, parameter :: nsmax = 300          ! maximum number of sites
      character*10 :: sites(nsmax)
      integer :: idall(nsmax), nyst(nsmax), nyed(nsmax)
      real*8  :: lons(nsmax), lats(nsmax)
      real*8  :: lons_old, lats_old
      integer :: nsite, irun(nsmax)
      integer :: start, start2(2),ncut2(2), ncid
      integer :: dimlen, jdmin, jdmax
      integer :: jdayst, jdayed, jdnow, jhnow
      integer :: jh_old  = -1
      integer :: iflag   = -1
      logical :: ifirst  = .true.
      logical :: recycle = .false.

      end module site_com


      subroutine read_site(site, lon, lat, nh, nd, nm, ny,
     &                     jday, year, jday2, year2,
     &                     tas, tcans, qs, ps, ws, 
     &                     prs, pfs, co2s, o3s, chs, fws, 
     &                     soilt, soilm, soili, soilp, iflag0)

      use site_com
   
      implicit none
      include 'netcdf.inc'
 
      character*(*), intent(in) :: site
      real*8, intent(in)  :: lon, lat
      integer, intent(in) :: nh, nd, nm, ny
      integer, intent(in) :: jday, year, jday2, year2
      real*8, intent(out) :: tas, tcans, qs, ps, ws, prs, pfs
      real*8, intent(out) :: co2s, o3s, chs, fws
      real*8, intent(out) :: soilt(6), soilm(6), soili(6), soilp(6)
      integer, intent(out):: iflag0
 
      character*120 dirs, fsite
      integer nstat, varid, dimid
      integer nd1, nm1, ny1
      real*4  nh1, jdays, jhours
      integer, parameter :: ncount=1
      integer, parameter :: mds(12) = (/
     &    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      integer lenc, n

      jdays = (ny-1990)*365+nd-mds(nm)
      do n = 1, nm
         jdays = jdays + mds(n)
      enddo
      jhours = (jdays-1)*24+nh
      if (jhours .ge. jh_old) then
         jh_old = jhours
      else
         recycle = .true.
         jh_old  = -1
      endif

      If (.not. ifirst) then
      If (abs(lons_old-lon) .gt. 1.0d-3 .or.
     &    abs(lats_old-lat) .gt. 1.0d-3) then 
          ifirst = .true.
          recycle = .false.
          iflag = -1
      Endif
      If (iflag .ne. -3 .and. recycle) then       ! end of site forcing file
         nstat=NF_CLOSE(ncid)
         IF(nstat .NE. NF_NOERR) Call handle_err(nstat)
      Endif
      If (iflag .eq. -3 .and. .not. recycle) then
         iflag0 = iflag
         return     ! end of site forcing file
      Endif
      Endif

      If (iflag .eq. -2) then 
          iflag0 = iflag
          return
      Endif
      
      If (ifirst .or. recycle) then
         lons_old = lon
         lats_old = lat
         jdmin = (year-1990)*365+jday              ! days refer to 1990/01/01
         jdmax = (year2-1990)*365+jday2             ! days refer to 1990/01/01
         ifirst     = .false.
         recycle    = .false.

         dirs  = '/home/YIBS_site/Input/DRIVER_DATA/'
         fsite = dirs(1:lenc(dirs))//site(1:lenc(site))//
     &          '.forcing.nc'
         print*, 'Read Site OK : '//site(1:lenc(site))//
     &          '.forcing.nc'
         nstat    = NF_OPEN(fsite(1:lenc(fsite)),NF_NOWRITE,ncid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)

         nstat    = NF_INQ_DIMID(ncid, 'time', dimid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_DIMLEN(ncid, dimid, dimlen)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)

         nstat    = NF_INQ_VARID(ncid,'year',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,1,ncount,ny1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'month',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,1,ncount,nm1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'dom',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,1,ncount,nd1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'start_hr',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,1,ncount,nh1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         jdayst   = (ny1-1990)*365+nd1-mds(nm1)
         do n = 1, nm1
            jdayst = jdayst + mds(n)
         enddo
         
         nstat    = NF_INQ_VARID(ncid,'year',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,dimlen,ncount,ny1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'month',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,dimlen,ncount,nm1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'dom',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,dimlen,ncount,nd1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         jdayed   = (ny1-1990)*365+nd1-mds(nm1)
         do n = 1, nm1
            jdayed = jdayed + mds(n)
         enddo
         print*, ny1, nm1, nd1, nh1
         print*, jdayst, jdayed, jdmin, jdmax

         if (jdayst .gt. jdmax .or. jdayed .lt. jdmin) then
            print*, 'No appropriate site forcing data found!'
            iflag = -2
            iflag0 = iflag
            nstat=NF_CLOSE(ncid)
            IF(nstat .NE. NF_NOERR) Call handle_err(nstat)
            return 
         endif

         if (jdayst .ge. jdmin) then
            start = 1
            jdnow = jdayst
         else
            start = 1
            jdnow = jdayst
            do while (jdnow .lt. jdmin) 
            start = start + 1
            nstat    = NF_INQ_VARID(ncid,'year',varid)
            If(nstat .NE. NF_NOERR) Call handle_err(nstat)
            nstat    = NF_GET_VARA_INT(ncid,varid,start,ncount,ny1)
            If(nstat .NE. NF_NOERR) Call handle_err(nstat)
            nstat    = NF_INQ_VARID(ncid,'month',varid)
            If(nstat .NE. NF_NOERR) Call handle_err(nstat)
            nstat    = NF_GET_VARA_INT(ncid,varid,start,ncount,nm1)
            If(nstat .NE. NF_NOERR) Call handle_err(nstat)
            nstat    = NF_INQ_VARID(ncid,'dom',varid)
            If(nstat .NE. NF_NOERR) Call handle_err(nstat)
            nstat    = NF_GET_VARA_INT(ncid,varid,start,ncount,nd1)
            If(nstat .NE. NF_NOERR) Call handle_err(nstat)
            jdnow   = (ny1-1990)*365+nd1-mds(nm1)
            do n = 1, nm1
               jdnow = jdnow + mds(n)
            enddo
            enddo
            start = start - 1
         endif
         jhnow = (jdnow-1)*24 - 1

      Endif

      if (jhnow .gt. jhours) then
         iflag = -1
         iflag0 = iflag
         return
      endif
      do while (jhnow .lt. jhours .and. start .le. dimlen)
         nstat    = NF_INQ_VARID(ncid,'year',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,start,ncount,ny1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'month',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,start,ncount,nm1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'dom',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,start,ncount,nd1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'start_hr',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,nh1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         jdnow   = (ny1-1990)*365+nd1-mds(nm1)
         do n = 1, nm1
            jdnow = jdnow + mds(n)
         enddo
         jhnow = (jdnow-1)*24+nh1
         start = start + 1
      enddo
      If (jhnow .eq. jhours) then
         iflag = 1
         iflag0 = iflag
         start  = max(start - 1, 1)
         start2(1) = 1
         start2(2) = start
         ncut2(1)  = 6
         ncut2(2)  = 1
         nstat    = NF_INQ_VARID(ncid,'year',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,start,ncount,ny1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'month',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,start,ncount,nm1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'dom',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_INT(ncid,varid,start,ncount,nd1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'start_hr',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,nh1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(ncid,'Tair',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,t1s)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         tas = t1s 
         nstat    = NF_INQ_VARID(ncid,'Tcanopy',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,t1c)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         tcans = t1c 
         nstat    = NF_INQ_VARID(ncid,'Qcanopy',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,q1s)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         qs = q1s
         nstat    = NF_INQ_VARID(ncid,'Psurf',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,p1s)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         ps = p1s
         nstat    = NF_INQ_VARID(ncid,'Wind',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,w1s)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         ws = w1s
         nstat    = NF_INQ_VARID(ncid,'Pardif',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,pf1s)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         pfs = pf1s
         nstat    = NF_INQ_VARID(ncid,'Pardir',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,pr1s)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         prs = pr1s
         nstat    = NF_INQ_VARID(ncid,'CO2',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,co21)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         co2s = co21
         nstat    = NF_INQ_VARID(ncid,'O3',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,o31)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         o3s = o31
         nstat    = NF_INQ_VARID(ncid,'Ch',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,ch1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         chs = ch1
         nstat    = NF_INQ_VARID(ncid,'fw',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start,ncount,fw1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         fws = fw1
         nstat    = NF_INQ_VARID(ncid,'Soiltemp',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start2,ncut2,st1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         soilt = st1
         nstat    = NF_INQ_VARID(ncid,'Soilmoist',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start2,ncut2,sm1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         soilm = sm1
         nstat    = NF_INQ_VARID(ncid,'Soilmp',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start2,ncut2,sp1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         soilp = sp1
         nstat    = NF_INQ_VARID(ncid,'fice',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VARA_REAL(ncid,varid,start2,ncut2,si1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         soili = si1
         If (tas .lt. 0) iflag = 0
         If (co2s .lt. 0) iflag = 0
         If (start .eq. dimlen) then       ! end of site forcing file
            nstat=NF_CLOSE(ncid)
            IF(nstat .NE. NF_NOERR) Call handle_err(nstat)
            iflag = -3
         Endif
      else
         iflag = -1
         return
      Endif
      iflag0 =iflag
     
      return
      End subroutine read_site
      
