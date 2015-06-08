      
      module vegfor_com
      
      implicit none
      save

      integer, parameter :: nlon=144
      integer, parameter :: nlat=90
      integer idx, jdy
      integer, parameter :: mdays(12) = (/
     &    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      integer, PARAMETER :: JMPERY = 12
      integer :: JDmidOfM(0:JMPERY+1) = (
     *     /-15,16,45,75,106,136,167,197,228,259,289,320,350,381/)

      real*8, dimension(nlon,nlat,16,12)  :: vlai
      real*8, dimension(nlon,nlat,16)     :: vfrac
      real*8, dimension(nlon,nlat,16)     :: vheight
      real*8, dimension(nlon,nlat,16)     :: vlaix
      real*8, dimension(nlon,nlat)        :: vcrop
      real*8, dimension(nlon,nlat,2)      :: crop_cal

      end module vegfor_com


      subroutine getlaix_m16(laimax,hdata)

      use vegfor_com
      implicit none

      real*8, intent(out) :: laimax(18), hdata(18)
      integer i
       
      laimax(:) = 0.0d0
      hdata(:)  = 0.0d0
      do i = 1, 16
         laimax(i) = vlaix(idx, jdy, i)
         hdata(i)  = vheight(idx, jdy, i)
      enddo

      return
      end subroutine getlaix_m16

  
      subroutine getlai_m16(jday, lai)

      use vegfor_com
      implicit none

      integer, intent(in)  :: jday
      real*8, intent(out)  :: lai(18)
      integer, dimension(0:JMperY+1) :: startday
      integer :: jmon, offset, itd, totdays, k, nm1, nm2
      real*8  :: alpha, beta

      startday(0:JMperY+1)=jdMIDofM(0:JMperY+1)-1

      jmon  = 1
      totdays = 0
      do while (totdays < jday)
         totdays = totdays + mdays(jmon)
         jmon = jmon + 1
      end do
      jmon = jmon - 1

      if(jday < startday(jmon))then
        offset=-1
        nm1 = mod(jmon+10, 12) + 1 
        nm2 = jmon
      else
        offset=0
        nm1 = jmon
        nm2 = mod(jmon,12) + 1
      endif
      itd = startday(jmon+1+offset) - startday(jmon+offset)

      beta  = dble(jday-startday(jmon+offset)) / dble(itd)
      alpha = 1.d0 - beta
      lai   = 0.0d0

      do k = 1, 16
         lai(k) = alpha*vlai(idx, jdy, k, nm1)
     &          + beta*vlai(idx, jdy, k, nm2)
      enddo

      end subroutine getlai_m16


      subroutine readlai_m16(year,lon,lat)
!@sum READDLAI read in leaf area indicies from selected file
      use vegfor_com
      use filemanager, only: openunit,closeunit

      implicit none
      
      integer, intent(in) :: year
      real*8, intent(in)  :: lon, lat
      integer :: k,iunit, n, i, j
      real*4  :: vlai1(nlon,nlat)
      real*8 lons(nlon), lats(nlat), lon0, dlon, dlat
      character*80 :: fbin, head
      character*2  :: cmonth

! read first month's LAI's:
      vlai(:,:,:,:) = 0.0d0
      do n = 1,12 
      write(cmonth, '(I2.2)') n
      fbin='LAIM16_'//cmonth
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,16
        read(iunit) head, vlai1
        vlai(:,:,k,n) = vlai1
      enddo
      call closeunit(iunit)
      enddo

      vfrac(:,:,:) = 0.0d0
      call openunit('VEG16',iunit,.true.,.true.)
      do k=1,16
        read(iunit) head, vlai1
        vfrac(:,:,k) = vlai1
      enddo
      call closeunit(iunit)
       
      vheight(:,:,:) = 0.0d0
      call openunit('VHT16',iunit,.true.,.true.)
      do k=1,16
        read(iunit) head, vlai1
        vheight(:,:,k) = vlai1
      enddo
      call closeunit(iunit)
       
      vlaix(:,:,:) = 0.0d0
      call openunit('VLAIX16',iunit,.true.,.true.)
      do k=1,16
        read(iunit) head, vlai1
        vlaix(:,:,k) = vlai1
      enddo
      call closeunit(iunit)

      vcrop(:,:) = 0.0d0
      call readcrop_m16(year,vcrop)
 
      do i = 1, nlon
         lons(i) = -178.75d0 + dble(i-1)*2.5d0
      enddo
      do j = 1, nlat
         lats(j) = -89.0d0 + dble(j-1)*2.0d0
      enddo
      lon0 = lon
      If (lon .gt. 180.0d0) lon0 = lon - 360.0d0
      dlon = lons(10) - lons(9)
      dlat = lats(10) - lats(9)
      idx  = nint((lon0-lons(1))/dlon)+1
      jdy  = nint((lat-lats(2))/dlat)+2

      return
      end subroutine readlai_m16


      subroutine readcrop_m16(year,cropdata)
!@sum Read in crop cover fractions from file.
!@+   Calculates crop fraction for given year.
      use vegfor_com
      use FILEMANAGER, only : openunit,closeunit,nameunit
      implicit none
      include 'netcdf.inc'

      integer, intent(in) :: year
      real*8, intent(out) :: cropdata(nlon,nlat)
      integer i, nstat, ncid, id1, id2
      !----------
      integer :: iu_CROPS
      integer :: year1, year2
      real*4 crop4(nlon,nlat)
      real*8 wt, crop1(nlon,nlat), crop2(nlon,nlat)
      real*4 date1(nlon,nlat), date2(nlon,nlat)
      character*80 title
      
      !* Calculate fraction for given gcmtime:  interpolate between years*/
        
      year1 = -32768 ; crop1(:,:) = 0.d0
      year2 = -32767 ; crop2(:,:) = 0.d0
      wt = 1.d0
          
      call openunit("VCROPS",iu_CROPS,.true.,.true.)
      do while( year2 < year )
        year1 = year2
        crop1(:,:) = crop2(:,:)
        read (iu_CROPS,end=10) title , crop4
        read(title,*) year2 !Read year integer out of character array title
        crop2 = crop4
      enddo
      wt = (year-year1)/(real(year2-year1,kind=8))
 10   continue
      call closeunit(iu_CROPS)

      cropdata(:,:) = max(0.d0, crop1(:,:)
     &     + wt * (crop2(:,:) - crop1(:,:))) 

      nstat=NF_OPEN('CROPS_CAL',NCNOWRIT,ncid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_INQ_VARID(ncid,'crop_plant_day',id1)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_GET_VAR_REAL(ncid,id1,date1)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_INQ_VARID(ncid,'crop_harvest_day',id2)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_GET_VAR_REAL(ncid,id2,date2)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_CLOSE(ncid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      crop_cal(:,:,1)=date1
      crop_cal(:,:,2)=date2

      return
      end subroutine readcrop_m16

      

      subroutine modify_lai(jday,npft,lai)

      implicit none

      integer, intent(in)  :: jday, npft
      real*8, intent(inout):: lai

      if (npft .le. 14) return   

! narrow the growth period for crop and grass

      if (jday .le. 90 .or. jday .gt. 330) lai = 0.05d0*lai
      if (jday .gt. 90 .and. jday .le. 150) lai = 0.1d0*lai
      if (jday .gt. 270 .and. jday .le. 330) lai = 0.1d0*lai
      if (jday .gt. 150 .and. jday .le. 180)
     &    lai = (0.1d0+0.03d0*dble(jday-150))*lai
      if (jday .gt. 240 .and. jday .le. 270)
     &    lai = (0.1d0+0.03d0*dble(270-jday))*lai

      return
  
      end subroutine modify_lai

