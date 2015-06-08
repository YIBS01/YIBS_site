      module yibs_prog_mod
!@sum Utilities for off-line run of YIBS model.

      use yibs_mod
      use FILEMANAGER
      implicit none

      contains

!************************************************************************

      subroutine yibs_init_vegstruct( cells,
     &     IM, JM, I0, I1, J0, J1, jday, year, lat2d, lon, lat, id_veg,
     &     do_soilinit,do_phenology_activegrowth)
      use yibs_prescribed_drv, only:init_canopy_physical,prescr_vegdata
     &     ,yibs_init_params
      !use yibs_prescr_veg, only : prescr_calcconst
      use yibs_const
      use vegfor_com, only: vfrac, vcrop, crop_cal, idx, jdy

      implicit none
      type(ycelltype_public), intent(inout) :: cells(I0:I1,J0:J1)
      integer, intent(in) :: IM, JM, I0, I1, J0, J1, jday, year
      real*8, dimension(I0:I1,J0:J1) :: lat2d
      real*8, intent(in)  :: lon, lat
      integer, intent(in) :: id_veg
      logical, intent(in) :: do_soilinit, do_phenology_activegrowth
      !---Local variables-----
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata2 !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      integer, dimension(N_COVERTYPES) :: soildata ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(I0:I1,J0:J1) :: plant_date,harvest_date
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &                  I0:I1,J0:J1):: Tpooldata  !in g/m2 
      real*8 :: cropdata2
      !------
      integer :: iu,i,j
      
      call yibs_init_params()

 
#ifdef PFT_MODEL_YIBS
      write(*,*) "PFT_MODEL_YIBS defined"
#endif

         vegdata2 = 0.0d0
         do i=I0,I1
         do j=J0,J1
            vegdata2(1:16,i,j) = vfrac(idx,jdy,:)
            plant_date(i,j)    = crop_cal(idx,jdy,1)
            harvest_date(i,j)  = crop_cal(idx,jdy,2)
         enddo
         enddo
         cropdata2 = vcrop(idx,jdy)

         call prescr_vegdata(jday, year, IM,JM,I0,I1,J0,J1,
     &     lon, lat, vegdata,albedodata,laidata,hdata,
     &     rootprofdata, soildata,soil_texture,
     &     Tpooldata, vegdata2, cropdata2,
     &     do_soilinit,do_phenology_activegrowth)

           vegdata = 0.d0
           vegdata(id_veg,I0:I1,J0:J1) = 1.0d0
           laidata = 0.d0
           laidata(id_veg,I0:I1,J0:J1) = 2.0d0     ! default, just for initial

         !Translate gridded data to YIBSdata structure
         call init_canopy_physical(I0, I1, J0, J1,
     &        Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

         call yibs_cell_set(cells,jday,lat2d,vegdata, laidata, hdata,
     &        rootprofdata, soildata, albedodata, soil_texture,
     &        Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpooldata,
     &        plant_date, harvest_date,
     &        reinitialize=.true.)  

      end subroutine yibs_init_vegstruct

!************************************************************************

      end module yibs_prog_mod
!************************************************************************


!************************************************************************
      !program yibs_prog
      subroutine YIBS_1D
      ! this driver just passes the bounds of the region, so that all
      ! arrays can be allocated on a stack
      use site_com

      implicit none
      integer IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer jday, year, jday2, year2
      real*8 dt !seconds
      real*8  lon, lat
      logical :: do_soilinit,do_soilresp
      logical :: do_phenology_activegrowth, do_structuralgrowth
      logical :: do_frost_hardiness
      logical :: do_patchdynamics
      integer :: id_veg, i, lenc
      character*10 site
      character*4 cyear1, cyear2 

      print *, "Started."

      !* Default configuration
      dt = 1800.d0              !second. Default value may be over-ridden by yibs_input
      do_soilinit = .true.
      do_soilresp = .true.
      do_phenology_activegrowth = .false.
      do_structuralgrowth = .false.
      do_frost_hardiness = .true.
      do_patchdynamics = .false.

      !* Set world bounds (should correspond to format of input files)
      IM = 72; JM = 46

      !* Default, entire grid.
      i0 = 1; i1 = 72
      j0 = 1; j1 = 46

      !* dims to default forcings file
      i0f = 1; i1f = 72
      j0f = 1; j1f = 46

      !* Default date start for GISS GCM 10-day forcings
      jday = 152                !June 1
      year = 1980
      jday2 = jday + 10
      year2 = -1 !Initialize

      write(*,*) 'Got here before read_input_parameters'
      call read_input_parameters(
     &     I0, I1, J0, J1, i0f, i1f, j0f, j1f,
     &     jday, year, jday2, year2,dt, lon, lat, id_veg, site)
      if (year2.eq.-1) year2=year  !If year2 was not specified, default same as year.

      if (.not.do_phenology_activegrowth .and.
     &   do_structuralgrowth) then
         print*,"impossible combinations of input parameters"
         stop
      endif

      print *,"yibs_input: "
     &     , jday, year, jday2, year2,dt
     &     , i0f, i1f, j0f, j1f
     &     , 'site location: ', lon, lat


125   format(1x,a10, 6x, i8, i8, i8, f13.4, f13.4, i3)

      print *,"starting program"

      Write(cyear1, '(i4)') year
      Write(cyear2, '(i4)') year2
      print*, "Simulating site "//site(1:lenc(site))//
     $        " from "//cyear1//" to "//cyear2

      call run_offline(IM, JM, I0, I1, J0, J1, jday, year,jday2,year2,dt
     &     ,i0f, i1f, j0f, j1f, lon, lat
     &     ,do_soilinit, do_soilresp, do_phenology_activegrowth
     &     ,do_structuralgrowth, do_frost_hardiness 
     &     ,do_patchdynamics, id_veg, site)

      print *,"YIBS run completed."

      end subroutine YIBS_1D

!************************************************************************
      subroutine read_input_parameters(
     &      I0, I1, J0, J1, i0f, i1f, j0f, j1f,
     &      jday, year, jday2, year2,dt, lon, lat, id_veg, site)

      use filemanager
      implicit none
      integer I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer jday, year,jday2, year2
      real*8 dt
      real*8 lon, lat
      !----
      integer iu_yibs_input, id_veg, i,j
      real*8 long(72), latg(46), dlon, dlat, lon0
      character*80, parameter :: file_yibs_input="yibs_input"
      character*10 site
      namelist /input_parameters/ 
     &      jday, year, jday2, year2,dt, lon, lat,
     &      id_veg, site

      print *,"reading input parameters from ", file_yibs_input
      call openunit(trim(file_yibs_input),iu_yibs_input,.false.,.true.)
      read(iu_yibs_input, NML=input_parameters, ERR=10)
      call closeunit(iu_yibs_input)

      do i = 1, 72
         long(i) = -177.5 + dble(i-1)*5.
      enddo
      do j = 1, 46
         latg(j) = -90.0 + dble(j-1)*4.
      enddo
      lon0 = lon
      If (lon .gt. 180.0d0) lon0 = lon - 360.0d0
      dlon = long(10) - long(9)
      dlat = latg(10) - latg(9)
      I0   = nint((lon0-long(1))/dlon)+1
      J0   = nint((lat-latg(2))/dlat)+2
      I1   = I0
      I0f  = I0
      I1f  = I0
      J1   = J0
      J0f  = J0
      J1f  = J0

      return
 10   continue
      print *,"error reading namelist file:", file_yibs_input
      stop 255
      end subroutine read_input_parameters


!************************************************************************
      subroutine run_offline(IM, JM, I0, I1, J0, J1
     &     ,jday, year, jday2, year2,dt
     &     ,i0f, i1f, j0f, j1f, lon ,lat
     &     ,do_soilinit
     &     ,do_soilresp, do_phenology_activegrowth
     &     ,do_structuralgrowth, do_frost_hardiness
     &     ,do_patchdynamics, id_veg, site)
      !****************************************************************
      !* Example program to run YIBS coupled to a GCM.
      !* This version assumes an explicit scheme for calculation of
      !* canopy conductance, photosynthesis, and temperature.
      !* - For parallelization, the GCM provides the grid bounds to YIBS.
      !* - YIBS initializes vegetation structure parameters within these bounds.
      !* - GCM provides initial surface meteorological state variables.
      !* - YIBS initializes GCANOPY, Ci, Qf, zeroes GPP
      !* - Simulation loop:
      !*   a. GCM provides meteorological drivers and canopy temperature
      !*   b. YIBS updates.
      !*   c. Program gets from YIBS: GCANOPY, Ci, Qf, GPP, TRANS_SW, albedo
      !*   d. GCM updates Qf, Tcanopy, given GCANOPY, TRANS_SW, albedo
      !*      and runs tracers on C and N.
      !****************************************************************

      use yibs_mod
      use yibs_prog_mod
      use filemanager
      use yibs_const 
      
      implicit none
      integer, intent(in) :: IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer, intent(in) :: jday, year,jday2,year2
      real*8,  intent(in) :: lon, lat
      logical, intent(in) :: do_soilinit
      logical, intent(in) :: do_soilresp,do_frost_hardiness
      logical, intent(in) :: do_phenology_activegrowth
      logical, intent(in) :: do_structuralgrowth
      logical, intent(in) :: do_patchdynamics

      integer, intent(in) :: id_veg
      character*(*),intent(in) :: site
      real*8, intent(in) :: dt
      !---Local----
      integer :: jdaycount  
      integer :: nh, nd, nm, ny, totdays, nymet, nm_old=-1
      integer, parameter :: mdays(12) = (/
     &    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      
      logical :: do_rewind

      real*8 :: max_time !In seconds
      real*8,parameter :: max_days = 365 !In days
      real*8, parameter :: save_interval=86400.d0*10d0 !save every 10 days

      !Coupling variables
      !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
      real*8 :: tas, tcans, qs, ps, ws, prs, pfs, co2s, o3s, chs, fws
      real*8 :: td, cosz1, cosz2, cosz
      real*8 :: rtmp
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) :: soilt
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) :: soilm
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) :: soilp
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) :: soili
      real*8, dimension(I0:I1,J0:J1) ::
     &     lat2d                !2-d latitude - for phenology
     &     ,TairC               !Air temperature (Celsius) !KIM - for phenology
     &     ,TcanopyC            !Canopy temperature (Celsius)
     &     ,Qf                  !Foliage surface specif humidity (kg vapor/ kg air)
     &     ,P_mbar              !Atmospheric pressure (mb)
     &     ,Ca                  !@Atmos CO2 conc at surface height (mol/m3).
     &     ,Ch                  !Ground to surface heat transfer coefficient 
     &     ,U                   !Surface layer wind speed (m s-1)
     &     ,IPARdif             !Incident diffuse PAR (vis) (W m-2)
     &     ,IPARdir             !Incident direct PAR (vis) (W m-2)
     &     ,IPARfrac            ! IPARdif/(IPARdif+IPARdir)
     &     ,CosZen              !cos of solar zenith angle
     &     ,O3                  !Surface [O3] (nmol/m3)
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) ::  !Changed to GCM layering -NK
     &      Soiltemp           !soil temp 
     &     ,Soilmoist          !soil volum moist 
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) ::
     &     Soilmp              !Soil matric potential
     &     ,fice                !Fraction of soil water that is ice.

      !Coupling and diagnostic variables specific to YIBS (output)
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidatam
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: fv0
      real*8, dimension(I0:I1,J0:J1) ::
     &     GCANOPY              !Canopy conductance of water vapor (m s-1). 
     &     ,Ci                  !Internal foliage CO2 concentration (mol/m3)
     &     ,NPP                 !Net primary productivity (kg[C]/m2/s).
     &     ,GPP                 !Gross primary productivity (kg[C]/m2/s).
     &     ,GPPD                !Gross primary productivity (kg[C]/m2/s).
     &     ,FO3                 !ozone flux to stomata (nmol O3 m-2 s-1).
     &     ,dFO3                !excess ozone flux to stomata (nmol O3 m-2 s-1).
     &     ,TRANS_SW            !Transmittance of shortwave through canopy to soil 
     &     ,z0, CO2flux         !Variables from YIBS to GCM
     &     ,R_auto              !Variables from YIBS to GCM
     &     ,lai0                !Variables from YIBS to GCM
     &     ,lai_p               !Prognostic LAI
     &     ,ht_p                !Prognostic HT
     &     ,resp_r              !Root respiration
     &     ,resp_l              !Leaf respiration
     &     ,resp_w              !Wood respiration
     &     ,carb_tot
     &     ,carb_soil
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) ::
     &     betadl
      real*8, dimension(N_BANDS,I0:I1,J0:J1) ::
     &     albedo               !Variables from YIBS to GCM
      !---------------------------------------------------------------------

      type(ycelltype_public) :: cells(I0:I1,J0:J1)
      integer lenc
      real*8 time, time_since_last_save
      integer hemi(I0:I1,J0:J1) !hemisphere flags = 1 for N., =-1 for S.
      logical :: update_day    !For prescribed phenology, litter
      real*8 :: fw(I0:I1,J0:J1) ! fraction of wet canopy
      real*8 fv_tot

      !---Other local vars
      integer :: i,j,iu_ycells,iflag, iflagm !For printing ycells - NK

      print *,"started run_offline with:", IM, JM, I0, I1, J0, J1

      fw(:,:) = 0.d0   ! all canopy is dry


      !* Now YIBS should be initialized before any other calls to yibs_*
      call yibs_init_config(
     &     do_soilresp=do_soilresp
     &     ,do_phenology_activegrowth=do_phenology_activegrowth
     &     ,do_structuralgrowth=do_structuralgrowth
     &     ,do_frost_hardiness=do_frost_hardiness
     &     ,do_patchdynamics=do_patchdynamics)

      !* Set hemisphere flags.
      if ( J0<=JM/2 )   hemi(:,J0:min(JM/2,J1))   = -1    ! S.
      if ( J1>=JM/2+1 ) hemi(:,max(JM/2+1,J0):J1) =  1    ! N.
      print *,"set hemi ok"
      
      !* Initialize yibs cells (makes head ycell).
      call yibs_cell_construct( cells )
      print *,"yibs_cell_construct passed ok" 


      !* Initialize vegetation cover.  
      call readlai_m16(year,lon,lat)
      lat2d = lat
      call yibs_init_vegstruct( cells, IM, JM, I0, I1, J0, J1, 
     &        jday, year, lat2d, lon, lat, id_veg,
     &        do_soilinit,do_phenology_activegrowth)
      print *,"yibs_init_vegstruct passed ok"
    

      max_time = 86400.*(365*(year2-year) + jday2)
      ny   = year
      print *,"max_time: ",max_time
      time = 86400.*(jday-1) !Allows starting run in middle of year at jday1.
      time_since_last_save = 0.d0
      jdaycount = jday
      update_day = .true. !Initialize

      call open_output_netcdf(site, lon, lat, year)
      print *,"---------------------------------------------------"
      print *,"started time step, time=", time, "dt=", dt

      do while( time < max_time )

        totdays = 0
        nm  = 1
        do while (totdays < jdaycount) 
           totdays = totdays + mdays(nm)
           nm = nm + 1
        end do
        nm = nm - 1
        totdays = totdays - mdays(nm)
        nd = jdaycount - totdays
        nh = mod(time/3600,24.)
        nymet = max(ny, year)
        print *,"hour, day, month, year =", nh, nd, nm, ny

        td = jdaycount+(dble(nh))/24.0d0        
        call calc_solarzen(td,lat,cosz1)
        td = jdaycount+(dble(nh)+0.5)/24.0d0        
        call calc_solarzen(td,lat,cosz2)
        cosz = (cosz1+cosz2)/2.0d0

        call  read_site(site, lon, lat, nh, nd, nm, nymet,
     &                  jday, year, jday2, year2,
     &                  tas, tcans, qs, ps, ws,  
     &                  prs, pfs, co2s, o3s, chs, fws,
     &                  soilt, soilm, soili, soilp, iflag)
        If (iflag .ge. 0) then 
           If (tas .ge. -99.) then 
              TairC = tas
              TcanopyC = tcans
              Qf    = qs
              P_mbar= ps
              Ca    = co2s*(1.0D-06)*P_mbar*100.0/gasc/(TcanopyC+tfrz)
              O3    = o3s*P_mbar*100.0/gasc/(TcanopyC+tfrz)
              Ch    = chs
              U     = ws
              IPARdir = prs
              IPARdif = pfs
              CosZen  = cosz
              fw = fws
              soiltemp  = soilt
              soilmoist = soilm
              soilmp    = soilp
              fice      = soili
            Endif
        Endif

        !* Set forcings

        call yibs_set_forcings( cells,
     &       air_temperature=TairC,
     &       canopy_temperature=TcanopyC,
     &       canopy_air_humidity=Qf,
     &       surf_pressure=P_mbar,
     &       surf_CO2=Ca,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &       surf_O3=O3,
#endif
     &       heat_transfer_coef=Ch,
     &       wind_speed=U,
     &       total_visible_rad=IPARdir+IPARdif,
     &       direct_visible_rad=IPARdir,
     &       cos_solar_zenith_angle=CosZen,
     &       canopy_wet_fraction=fw,
     &       soil_temp=Soiltemp,  
     &       soil_moist=Soilmoist,
     &       soil_matric_pot=Soilmp,
     &       soil_ice_fraction=fice
     &       )

      !* NEW STREAMLINED CONTROL *!
      if (update_day) then
          laidata = 0.d0
          call getlai_m16(jdaycount, laidata)
          call yibs_prescribe_vegupdate(cells,hemi,jdaycount,year,
     &         do_giss_phenology=.false.,
     &         do_giss_albedo=.true.,
     &         do_giss_lai=.false.,
     &         laidata=laidata,
     &         update_crops=.false.,init=.false.)
          print *,'yibs_prescribe_vegupdate with calculated lai, ht'
      endif                 


      call yibs_run(cells,dt,update_day) 
 
 
        call yibs_get_exports( cells,
     &       canopy_conductance=GCANOPY,
     &       beta_soil_layers=betadl,
     &       shortwave_transmit=TRANS_SW,
     &       leafinternal_CO2=Ci,
     &       foliage_humidity=Qf,
     &       canopy_npp=NPP,
     &       canopy_gpp=GPP,
     &       canopy_gpp0=GPPD,
     &       canopy_resp_r=Resp_r,
     &       canopy_resp_l=Resp_l,
     &       canopy_resp_w=Resp_w,
     &       canopy_carb_tot=Carb_tot,
     &       canopy_carb_soil=Carb_soil,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &       ozone_flux=FO3,
     &       excess_ozone_flux=dFO3,
#endif
#ifdef ACTIVE_GROWTH
     &       LAI_prognostic=lai_p,
     &       height_prognostic=ht_p,
#endif
     &       roughness_length=z0,
     &       flux_CO2=CO2flux,
     &       R_auto=R_auto,
     &       albedo=albedo
     &     ,  leaf_area_index=lai0
     &     ,  vegetation_fractions=fv0
     &     )
  
        fv_tot = sum(fv0)
        if (ny .ge. year) then
        call write_output_netcdf(GCANOPY, Ci, NPP, GPP, GPPD, CO2flux, 
     &       resp_r, resp_l, resp_w, carb_tot, carb_soil,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &       FO3, dFO3,
#endif
#ifdef ACTIVE_GROWTH
     &       lai_p, ht_p,
#endif
     &       R_auto, coszen, O3s, lai0, fv_tot, ny, nm, nd, nh, year)
        endif

        time = time + dt  !Moved time update to before check of update_day-NK
        update_day=(dmod(time,86400.d0) .eq. 0.d0).and.(time.ne.0.d0)

        if (update_day) then
          jdaycount = jdaycount + 1
          if (jdaycount>365) then !numdays=365 for 1-year data, etc.
            jdaycount = 1       !reset
            ny = ny + 1
          endif
        endif
        
      enddo

      call close_output_netcdf

      end subroutine run_offline

!***************************************************************************

