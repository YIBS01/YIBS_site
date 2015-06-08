#include "rundeck_opts.h"

      module yibs_prescribed_drv

      !*********************************************************************
!@sum !*    SUBROUTINES TO READ IN prescribed VEGETATION DATA SETS 
!@+   !+      for initialization and updating prescribed cover like crops.
!@+   !*    Array data only, no ycells or patches info.
!@+   !*    Interfaces with yibs_prescr_veg for Matthews pft-level calculations
!@+   !*      or with yibs_prescribed_drv_geo for geographic initialization.
      !*********************************************************************

      use yibs_const
      use yibs_pfts
      use yibs_prescr_veg

      implicit none
      private
      save

      public 
     &     yibs_init_params,
     &     init_canopy_physical,
     &     prescr_vegdata,
     &     prescr_veg_albedodata

      public prescr_get_soilpools
      public prescr_get_laidata, prescr_get_soil_C_total
      public prescr_get_hdata
      public prescr_get_pft_vars

      contains

!***************************************************************************
      subroutine yibs_init_params()
!@sum Initialize some YIBS parameters.
      use yibs_prescr_veg, only : init_params

      !call prescr_calcconst !renamed init_params
      call init_params()

      end subroutine yibs_init_params
!***************************************************************************
      subroutine init_canopy_physical(
     & I0,I1,J0,J1,Ci_ini, CNC_ini, Tcan_ini, Qf_ini)
!@sum For old Friend & Kiang (2005) biophysics. Initialize LSM outputs.
      integer,intent(in) :: I0,I1,J0,J1
      real*8, DIMENSION(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini

      Ci_ini(:,:) = 0.0127d0
      CNC_ini(:,:) = 0.d0
      Tcan_ini(:,:) = 0.d0        !Should be a forcing from land surface model.
      Qf_ini(:,:) = 0.d0          !Should be a forcing from land surface model.

      end subroutine init_canopy_physical
      
!***************************************************************************
      subroutine prescr_get_soilpools(I0,I1,J0,J1, lon, lat,
     &     soil_C_total, Tpool_ini)
!@sum Prescribe initial partitioning of soil carbon pools given total.
      !* For global runs, this routine gets soil_C_total from subroutine 
      !get_soil_C_total,which reads in total soil pool amounts (measured). 
      !Then individual soil pool fractions (modeled pft-dependent values 
      !from spinup runs by PK),are used to prescribe individual amounts.
      !**all carbon amounts should be in g/m2** -PK 12/07, NK 7/11
      !
      !* For site runs, this routine reads in all soil carbon fractions.

      use FILEMANAGER, only : openunit,closeunit
      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(in) :: lon, lat
      real*8,intent(in) ::
     &     soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      real*8,intent(out) :: 
     &      Tpool_ini(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,  
     &                I0:I1,J0:J1)!prescribed soil pools, g/m2
      !-----Local------
!      first 3 for eventually reading in globally gridded dataset, e.g. ISRIC-WISE
      integer :: iu_SOILCARB, iunit
      integer :: n,p,nn,k,i,j,ii,jj
      real*8, dimension(N_PFT,NPOOLS-NLIVE,N_CASA_LAYERS) :: Cpool_fracs
      real*4 :: cpool1(360,181)
      real*8 :: cpool_pft(360,181,N_PFT), lons(360), lats(181)
      real*8 :: dlon, dlat, lon0
      character*80 :: head

      Tpool_ini(:,:,:,:,:,:) = 0.d0  !initialize all pools to zero (g-C/m^2)

#ifdef PFT_MODEL_YIBS
!YK - temp. values, modified from 8 GISS pfts below
!NK - later these arrays should be moved to yibs_pfts
        Cpool_fracs(1,:,1) = (/ !ever_ES_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(2,:,1) = (/ !ever_LS_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(3,:,1) = (/ !ever_ES_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(4,:,1) = (/ !ever_LS_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(5,:,1) = (/ !cold_ES_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(6,:,1) = (/ !cold_LS_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(7,:,1) = (/ !drought_broad
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(8,:,1) = (/ !decid_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(9,:,1) = (/ !shrub_cold
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(10,:,1) = (/ !shrub_arid
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(11,:,1) = (/ !c3grass
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(12,:,1) = (/ !c4grass
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(13,:,1) = (/ !c3grass_ann
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(14,:,1) = (/ !c3grass_arctic
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(15,:,1) = (/ !cropsc4
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(16,:,1) = (/ !cropstree
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
#else
!***  for now define for 8 GISS pfts, one-layer only -PK 1/23/08***
        Cpool_fracs(1,:,1) = (/ !tundra (for now=C3 grass)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(2,:,1) = (/ !C3 grass (Vaira)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(3,:,1) = (/ !shrub (for now=savanna)
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(4,:,1) = (/ !savanna (Tonzi)
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(5,:,1) = (/ !decid broadl (MMSF)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(6,:,1) = (/ !evergr needl (for now=decid broadl)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(7,:,1) = (/ !trop rainf (for now=decid broadl)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(8,:,1) = (/ !crops (for now=C3 grass)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
#endif

#ifdef RESTART_CPOOLS
      print*, "Restart Soil Carbon Pools"
      call openunit("SOILCARB_restart",iunit,.true.,.true.)
      do k=1,N_PFT
        read(iunit) head, cpool1
        cpool_pft(:,:,k) = cpool1
      enddo
      call closeunit(iunit)
      
      do i = 1, 360
         lons(i) = -180.0d0 + dble(i-1)*1.d0 + 1.d0/2.d0
      enddo
      do j = 1, 181
         lats(j) = -90.0d0 + dble(j-1)*1.0d0
      enddo

      dlon = lons(10) - lons(9)
      dlat = lats(10) - lats(9)
      lon0 = lon
      If (lon .gt. 180.0d0) lon0 = lon - 360.0d0
      II   = nint((lon0-lons(1))/dlon)+1
      JJ   = nint((lat-lats(2))/dlat)+2

#endif

      !assign Tpool_ini values (pft-specific)
      do p=1,N_PFT
        do n=1,N_CASA_LAYERS 
          do nn=NLIVE+1,NPOOLS
#ifdef RESTART_CPOOLS
            Tpool_ini(p,CARBON,nn-NLIVE,n,I0:I1,J0:J1) =
     &           Cpool_fracs(p,nn-NLIVE,n)
     &           * cpool_pft(II,JJ,p)*1.d3
#else
            Tpool_ini(p,CARBON,nn-NLIVE,n,I0:I1,J0:J1) =
     &           Cpool_fracs(p,nn-NLIVE,n)
     &           * soil_C_total(n,I0:I1,J0:J1)*1.d3
#endif
          end do
        end do
      end do
ccc#endif

      end subroutine prescr_get_soilpools


      subroutine prescr_get_soil_C_total(I0,I1,J0,J1,lon,lat,
     &     soil_C_total)
!@sum (gC/m2) Read map of total soil carbon from file.
      use FILEMANAGER, only : openunit,closeunit
      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(in)  :: lon, lat
      real*8,intent(out) ::
     &     soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      !---
      real*4 :: buf(N_CASA_LAYERS,144,90)
      character*80 :: title
      integer :: iu_SOILCARB
      real*8 :: lons(144), lats(90)
      real*8 :: dlon, dlat, lon0
      integer k,i,j, ii, jj

      call openunit("SOILCARB_global",iu_SOILCARB,.true.,.true.)
      read (iu_SOILCARB) title, buf !data in kg/m2 (converted to g/m2 below)
      call closeunit(iu_SOILCARB)

      do i = 1, 144
         lons(i) = -178.75d0 + dble(i-1)*2.5d0
      enddo
      do j = 1, 90
         lats(j) = -89.0d0 + dble(j-1)*2.0d0
      enddo

      dlon = lons(10) - lons(9)
      dlat = lats(10) - lats(9)
      lon0 = lon
      If (lon .gt. 180.0d0) lon0 = lon - 360.0d0
      II   = nint((lon0-lons(1))/dlon)+1
      JJ   = nint((lat-lats(2))/dlat)+2
      do k=1,N_CASA_LAYERS
          soil_C_total(k,I0:I1,J0:J1) = buf(k,ii,jj)
      enddo

      end subroutine prescr_get_soil_C_total


      subroutine read_soilcarbon_site(I0,I1,J0,J1,Tpool_ini)
!@sum (gC/m2 by soil layer) Read site values for soil carbon pools from file.
!External file should be named as below and should be organized as follows:
!(1) there should be 1 or 2 columns (corresponding to each soil bgc layer);
!(2) first non-header row should have total site-measured pool (in g/m2);
!(3) 9 subsequent rows correspond to modeled 9 soil pool fractions
      use FILEMANAGER, only : openunit,closeunit,nameunit

      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(out) :: Tpool_ini(N_PFT,PTRACE,NPOOLS-NLIVE
     &     ,N_CASA_LAYERS, I0:I1,J0:J1)!prescribed soil pools, g/m2
      !---Local------------------------
      !Site total measured soil C_org:
      integer :: iu_SOILCARB  !File ID
      real*8, dimension(N_CASA_LAYERS) :: total_Cpool  !g-C/m^2
      real*8, dimension(NPOOLS-NLIVE,N_CASA_LAYERS) :: Cpool_fracs_in !fraction

      !Variables for calculating soil carbon pools
      real*8, dimension(N_PFT,NPOOLS-NLIVE,N_CASA_LAYERS) :: Cpool_fracs !fraction
      integer :: nn,p,n

      call openunit("SOILCARB_site",iu_SOILCARB,.false.,.true.)  !csv dataset
      read(iu_SOILCARB,*)  !skip optional header row(s)
      read(iu_SOILCARB,*) total_Cpool(:)
      do nn=1,NPOOLS-NLIVE
        read(iu_SOILCARB,*) Cpool_fracs_in(nn,:)
      end do

      do p=1,N_PFT      
       do n=1,N_CASA_LAYERS 
        do nn=NLIVE+1,NPOOLS
         Cpool_fracs(p,nn-NLIVE,n) = Cpool_fracs_in(nn-NLIVE,n)
         Tpool_ini(p,CARBON,nn-NLIVE,n,I0:I1,J0:J1) =
     &          Cpool_fracs(p,nn-NLIVE,n)*total_Cpool(n)
        end do
       end do
      end do
      end subroutine read_soilcarbon_site

!***************************************************************************
      ! This module is not used for GCM runs, but only for standalone runs.
      subroutine prescr_vegdata(jday, year, IM,JM,I0,I1,J0,J1,
     &     lon, lat, vegdata,albedodata,laidata,hdata,
     &     rootprofdata,soil_color,soil_texture,
     &     Tpooldata, vegdata2, cropdata2,
     &     do_soilinit,do_phenology_activegrowth)
!@sum prescr_vegdata - File reading and Matthews prescribed calculations
!@+   of vegetation structure. Off-line (standalone) runs only.
!     This is a general driver routine that calls routines
!     in this module only to read or calculate vegetation and soil structure
!     and carbon pools; the module subroutines call routines in yibs_prescr_veg,
!     which do the explicit calculations and are not available to outside
!     drivers.  
!     This routine may be imitated by drivers used for off-line runs
!     or coupled runs to GCMs.
      implicit none
      integer,intent(in) :: jday, year
      integer,intent(in) :: IM,JM,I0,I1,J0,J1 !long/lat grid number range
      real*8,intent(in)  :: lon, lat
      real*8,intent(in)  :: vegdata2(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(in)  :: cropdata2
      real*8,intent(out) :: vegdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: albedodata(N_BANDS,N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: laidata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: hdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH)
      integer,intent(out) :: soil_color(N_COVERTYPES)
      real*8,intent(out) :: soil_texture(N_SOIL_TEXTURES,I0:I1,J0:J1)
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpooldata !in g/m2 -PK
      logical,intent(in) :: do_soilinit
      logical,intent(in) :: do_phenology_activegrowth

      !-----Local------
      integer :: i,j,k, jeq, p
      integer hemi(I0:I1,J0:J1)
      REAL*8 :: soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      real*8 :: laimaxdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8 :: laimaxdata1(N_COVERTYPES,I0:I1,J0:J1)
      real*8 :: laidata1(N_COVERTYPES,I0:I1,J0:J1)
      real*8 :: hdata1(N_COVERTYPES,I0:I1,J0:J1)

      jeq = JM/2
      do j=J0,J1
        hemi(:,j) = 1

        if (j <= jeq) hemi(:,j) = -1
      enddo

      call init_vfraction(IM,JM,I0,I1,J0,J1,vegdata,vegdata2)   !veg fractions
      call prescr_update_vegcrops(year,IM,JM,I0,I1,J0,J1
     &     ,vegdata,cropdata2)

         laimaxdata = 0.d0
         hdata = 0.d0
         laidata = 0.d0
         call getlaix_m16(laimaxdata1,hdata1)
         call getlai_m16(jday, laidata1)
         call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata) !lai
         call prescr_get_hdata(I0,I1,J0,J1,hdata)
         call prescr_get_laimaxdata(I0,I1,J0,J1,laimaxdata)
         do j=J0,J1
           do i=I0,I1
             do k=1,16
               if (laidata1(k,i,j) .gt. 0.d0) then
                  laidata(k,i,j) = laidata1(k,i,j)
               endif
               if (hdata1(k,i,j) .gt. 0.d0) then
                  hdata(k,i,j) = hdata1(k,i,j)
               endif
               if (laimaxdata1(k,i,j) .gt. 0.d0) then
                  laimaxdata(k,i,j) = laimaxdata1(k,i,j)
               endif
             enddo
           enddo
         enddo

      call prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)

      call prescr_get_pft_vars(rootprofdata,soil_color)

      call prescr_get_soiltexture(I0,I1,J0,J1,lon,lat,soil_texture)
      Tpooldata(:,:,:,:,:,:) = 0.d0
      if ( do_soilinit ) then
#ifdef SOILCARB_SITE
         !Site soil carbon pools
           print *,"Getting site soil carbon"
           call read_soilcarbon_site(I0,I1,J0,J1,Tpooldata)
#else
         !Global soil carbon pools
          print *,'Reading global soil carbon data'
          call prescr_get_soil_C_total(I0,I1,J0,J1,lon,lat,soil_C_total)
          call prescr_get_soilpools(I0,I1,J0,J1,lon,lat,
     &                              soil_C_total,Tpooldata)
#endif
      endif

      end subroutine prescr_vegdata


!***************************************************************************
      subroutine init_vfraction(im,jm,I0f,I1f,J0f,J1f,vfraction,vfrac)
!@sum Read in vegetation and soil cover from file (mosaicked fractions).
      use FILEMANAGER, only : openunit,closeunit,nameunit
      integer, intent(in) :: im,jm,I0f,I1f,J0f,J1f
      real*8, intent(in)  :: vfrac(N_COVERTYPES,I0f:I1f,J0f:J1f) 
      real*8, intent(out) :: vfraction(N_COVERTYPES,I0f:I1f,J0f:J1f) 
      !------Local---------------------
      !1    2    3    4    5    6    7    8    9   10   11    12
      !BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE GRAC4
      character*80 :: title
      real*4 :: buf(im,jm)
      integer :: iu_VEG
      integer :: k,i,j
      real*8 :: s

      ! Make sure that unused fractions are set to 0
      vfraction(:,:,:) = 0.d0
      do i=I0f,I1f
         do j=J0f,J1f
            vfraction(:,i,j) = vfrac(:,i,j)
         enddo
      enddo

      !NK - This looks like it could be setting water to SAND??
      ! make sure that veg fractions are reasonable
      do j=J0f,J1f
         do i=I0f,I1f
          do k=1,N_COVERTYPES
            ! get rid of unreasonably small fractions
            if ( vfraction(k,i,j) < 1.d-4 ) vfraction(k,i,j) = 0.d0
          enddo
          s = sum( vfraction(:,i,j) )
          if ( s > .3d0 ) then !Keep if at least 90% specified, scale to 1.
             vfraction(:,i,j) = vfraction(:,i,j)/s
          else if ( s < .1d0 ) then
             print *, "missing veg data at ",i,j,"assume bare soil",s
             vfraction(:,i,j) = 0.d0
             vfraction(COVER_SAND,i,j) = 1.d0
          else
c             call stop_model("Incorrect data in VEG file",255)
          endif
        enddo
      enddo

      print *,'End of init_vfraction'
      end subroutine init_vfraction

!**************************************************************************

      subroutine prescr_update_vegcrops(year,IM,JM,I0,I1,J0,J1,
     &     vegdata, vcrop)
!@sum Rescales natural vegetation cover fractions given new crop cover.
      integer,intent(in) :: year
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      real*8, intent(in) :: vcrop
      real*8, intent(inout) :: vegdata(N_COVERTYPES,I0:I1,J0:J1)
      !--------
      real*8,ALLOCATABLE,dimension(:,:) :: cropdata !grid array
      integer :: i,j
      real*8 crops_old

      ALLOCATE(cropdata(I0:I1,J0:J1))

      !* Loop *!
      print *,"Getting crop cover."
      cropdata(I0:I1,J0:J1) = vcrop

      !* If cropdata was prepared somewhere else, then cover is as simple as
      !* modifying the vegetation fractions.  Need to update cohort and
      !* patch summary variables.
      do j=J0,J1
        do i=I0,I1
          if ( cropdata(i,j) == 1.d0 ) then
            vegdata(:,i,j) = 0.d0
            vegdata(CROPS+COVEROFFSET,i,j) = 1.d0
          else
            crops_old = vegdata(CROPS+COVEROFFSET,i,j)
c            if ( crops_old == 1.d0 ) then
c              call stop_model("incompatible crops: old=1, new<1",255)
c            endif
            vegdata(:,i,j) = vegdata(:,i,j)
     $           * (1.d0-cropdata(i,j))/(1.d0-crops_old)
            vegdata(CROPS+COVEROFFSET,i,j) = cropdata(i,j)
          endif
        end do
        !write(*,*) vegdata(CROPS+COVEROFFSET,:,j)
      end do

      DEALLOCATE(cropdata)
 
      end subroutine prescr_update_vegcrops

!**************************************************************************

      subroutine prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata)
!@sum Returns prescr GCM leaf area index for entire grid and given jday.
!@+   Based on seasonal sinusoidally varying LAI by veg type from R&A 1997.
      use yibs_const,only : N_COVERTYPES
      integer,intent(in) :: jday
      integer, intent(in) :: I0,I1,J0,J1
      real*8 :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      integer :: hemi(I0:I1,J0:J1) !@var hemi =1 in N. hemisphere, =-1 in South
      !----------
      integer :: n !@var cover type
      integer i,j,jeq

      !jeq = JM/2

      do j=J0,J1
        !hemi = 1
        !if (j <= jeq) hemi = -1
        do i=I0,I1
          do n=1,N_COVERTYPES
            laidata(n,i,j) = prescr_calc_lai(n,jday,hemi(i,j))
          enddo
        enddo
      enddo

      !* Return lai for each vegetation type.
      end subroutine prescr_get_laidata

!**************************************************************************

      subroutine prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)
!@sum Assign seasonal prescribed albedo by cover type (Matthews, 1983).
      integer,intent(in) :: jday
      integer, intent(in) :: I0,I1,J0,J1
      real*8 :: albedodata(N_BANDS,N_COVERTYPES,I0:I1,J0:J1)
      integer :: hemi(I0:I1,J0:J1)
      !----------
      !integer :: pft
      integer :: ncov
      integer i,j,jeq
      
      !jeq = JM/2

      do j=J0,J1
        !hemi = 1
        !if (j <= jeq) hemi = -1
        do i=I0,I1
          do ncov = 1, N_COVERTYPES
            call prescr_veg_albedo(hemi(i,j),ncov,jday,
     &           albedodata(:,ncov,i,j))
          end do
        enddo
      enddo

      end subroutine prescr_veg_albedodata

!**************************************************************************


      subroutine prescr_get_laimaxdata(I0,I1,J0,J1
     &     ,laimaxdata)
!@sum Make arrays of prescribed LAI: alamax, alamin from Matthews (1983).
      !* This should access alamax and alamin via yibs_prescr_veg, but
      !* so many levels is annoying for a trivial assignment. ...NK
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(out) :: laimaxdata(N_COVERTYPES,I0:I1,J0:J1) 
      !real*8,intent(out) :: laimindata(N_COVERTYPES,I0:I1,J0:J1) 
      !--- Local ----
      integer :: i,j
      do i=I0,I1
         do j=J0,J1
            laimaxdata(:,i,j) = alamax(:)
            !laimindata(:,i,j) = alamin(:)
         enddo
      enddo
      end subroutine prescr_get_laimaxdata

!**************************************************************************
      subroutine prescr_get_hdata(I0,I1,J0,J1,hdata) !height
!@sum Return array of prescribed canopy heights (m) by vegetation type
!@+   (Matthews, 1983).
      use yibs_const,only : N_COVERTYPES
      use yibs_prescr_veg, only : prescr_calc_hdata
      integer, intent(in) :: I0,I1,J0,J1
      real*8,intent(out) :: hdata(N_COVERTYPES,I0:I1,J0:J1)
      !---Local----
      integer :: i,j

      !* Read file of height data - TBA

      !* Calculate prescribed model tree heights
      do j=J0,J1
         do i=I0,I1
            call prescr_calc_hdata(hdata(:,i,j))
         enddo
      enddo
      
      end subroutine prescr_get_hdata

!**************************************************************************

      subroutine prescr_get_pft_vars(rootprofdata,soil_color)
!@sum Routine to assign pft-dependent misc. arrays that are not gridded.
      use yibs_prescr_veg, only : prescr_calc_initnm,
     &     prescr_calc_rootprof_all,prescr_calc_soilcolor

      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH)
      integer,intent(out) :: soil_color(N_COVERTYPES)

      call prescr_calc_rootprof_all(rootprofdata)
      call prescr_calc_soilcolor(soil_color)
      
      end subroutine prescr_get_pft_vars

!*************************************************************************

      subroutine prescr_get_soiltexture(I0,I1,J0,J1,lon,lat,
     &                                  soil_texture)
!@sum Return arrays of GISS soil color and texture.
      use FILEMANAGER, only : openunit,closeunit
      integer,intent(in)  :: I0,I1,J0,J1
      real*8, intent(in)  :: lon, lat
      real*8, intent(out) :: soil_texture(N_SOIL_TEXTURES,I0:I1,J0:J1)
      !------
      real*8 :: buf(144,90,N_SOIL_TEXTURES)
      real*8 :: lons(144), lats(90)
      real*8 :: dlon, dlat, lon0
      integer :: iu_SOIL
      integer k,i,j, ii, jj

      call openunit("soil_textures",iu_SOIL,.true.,.true.)
      !print *,IM,JM,N_COVERTYPES !#DEBUG
      read(iu_SOIL) buf
      call closeunit(iu_SOIL)
      !print *,"soil fractions:",buf(I0,J0,:)!#DEBUG
      do i = 1, 144
         lons(i) = -178.75d0 + dble(i-1)*2.5d0
      enddo
      do j = 1, 90
         lats(j) = -89.0d0 + dble(j-1)*2.0d0
      enddo

      dlon = lons(10) - lons(9)
      dlat = lats(10) - lats(9)
      lon0 = lon
      If (lon .gt. 180.0d0) lon0 = lon - 360.0d0
      II   = nint((lon0-lons(1))/dlon)+1
      JJ   = nint((lat-lats(2))/dlat)+2
      do k=1,N_SOIL_TEXTURES
          soil_texture(k,I0:I1,J0:J1) = buf(ii,jj,k)
      enddo

      end subroutine prescr_get_soiltexture
!*************************************************************************
      end module yibs_prescribed_drv

