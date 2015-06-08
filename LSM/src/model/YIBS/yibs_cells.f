      module ycells
!@sum ycells contains routines to perform summaries at the yibs grid cell
!@+   level.


      !YIBS MODULES TO USE
      use yibs_const
      use yibs_types
      use yibs_pfts
      use patches
      use cohorts

      implicit none
      private
      save

      public zero_ycell, summarize_ycell
      public assign_ycell, assign_ycell_soilcarbon
      public init_simple_ycell, ycell_construct, ycell_destruct
      public ycell_extract_pfts, ycell_carbon

      contains
!**************************************************************************

      subroutine zero_ycell_patchsum(ecp)
!@sum Zeros ycell patch summary variables.
      implicit none
      type(ycelltype) :: ecp

      !--Cohort-level properties--!
      ecp%LAI = 0.d0
      ecp%LAIpft(:) = 0.d0
      ecp%HTpft(:) = 0.d0
      ecp%Cdead(:) = 0.d0
      ecp%Clive(:) = 0.d0
      ecp%Phenfpft(:) = 0.d0
      ecp%h = 0.d0
      ecp%lai_p = 0.d0
      ecp%ht_p  = 0.d0
      ecp%fracroot(:) = 0.d0
      ecp%Ci = 0.0127D0         !Internal foliage CO2 (mol/m3) non-zero
      ecp%GCANOPY = 0.d0
      ecp%GPP = 0.d0
      ecp%GPP0= 0.d0
      ecp%GPPpft(:) = 0.d0
      ecp%IPPpft(:) = 0.d0
      ecp%MTPpft(:) = 0.d0
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      ecp%FO3 = 0.d0
      ecp%dFO3= 0.d0
#endif
      ecp%IPP = 0.d0
      ecp%MTP = 0.d0
      ecp%NPP = 0.d0
      ecp%resp_r = 0.d0
      ecp%resp_l = 0.d0
      ecp%resp_w = 0.d0
      ecp%resp_p = 0.d0
      ecp%R_auto = 0.d0
      ecp%R_root = 0.d0  !PK 5/15/07
      ecp%carb_tot = 0.d0
      ecp%carb_soil = 0.d0

      !--Patch-level properties, to be calculated from patch averages --!
      ecp%z0 = 0.d0
      ecp%albedo(:) = 0.d0
      ecp%betad = 0.d0
      ecp%betadl(:) = 0.d0
      ecp%TRANS_SW = 0.d0 
      ecp%CO2flux = 0.d0
      ecp%Soil_resp = 0.d0
      ecp%Tpool(:,:,:) = 0.d0

      ecp%C_total = 0.d0
      ecp%C_growth = 0.d0

      end subroutine zero_ycell_patchsum
!**************************************************************************

      subroutine zero_ycell(ecp)
!@sum Zeros import/export variables of entdata type
      implicit none
      type(ycelltype) :: ecp

      ecp%area = 0.d0

      call zero_ycell_patchsum(ecp)
      
      !VEGETATION - EXPORT STATE
      ecp%fv = 0.d0
      ecp%heat_capacity = 0.d0
      ecp%fwet_canopy = 0.d0
      !SOIL - IMPORT
      ecp%soil_Phi = 0.0         !Soil porosity (m3/m3)
      ecp%soildepth = 0.0        !Soil depth (m)
      ecp%theta_max = 0.0        !Saturated soil water volume (m/m)
      ecp%k_sat = 0.0            !Saturated hydraulic conductivity
      ecp%root_Phi = 0.0         !Infiltration factor promoted by roots (units?)

      !-----
      !METEOROLOGICAL - IMPORT STATE VARIABLES

      !VEGETATION - PRIVATE - Initial values not zero.
!      ecp%Ci = 0.0127D0         !Internal foliage CO2 (mol/m3) !!Cohort level or patch??
      ecp%Qf = 3.D-6            !Foliage surface vapor mixing ratio (kg/kg)
      
      !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
      ecp%jday = 0
      ecp%lat  = 0.d0
C Nadine
      ecp%plantdate = 0.d0          ! crop plant date
      ecp%harvestdate = 0.d0          ! crop harvest date

      ecp%TairC = 0.0           !Air temperature (Celsius) 
      ecp%TcanopyC = 0.d0       !Canopy temperatue (Celsius)
!      ecp%Qv = 0.0               !Canopy air specif humidity (kg vapor/ kg air)
      ecp%P_mbar = 0.d0         !Atmospheric pressure (mb)
      ecp%Ca = 0.d0             !@Atmos CO2 conc at surface height (mol/m3).
      ecp%Soilmoist(:) = 0.d0   !Soil moisture (volumetric fraction), depth-structured  -PK 6/28/06
      ecp%Soiltemp(:) = 0.d0    !Soil temp (Celsius), depth-structured 
      ecp%fice = 0.d0           !Fraction of soil layer that is ice
      ecp%Ch = 0.d0             !Ground to surface heat transfer coefficient 
      ecp%U = 0.d0              !Surface layer wind speed (m s-1)

      !Radiation - IMPORT STATE VARIABLE
      !may later be broken down into hyperspectral increments in an array
      ecp%IPARdif = 0.0         !Incident diffuse PAR 400-700 nm (W m-2)
      ecp%IPARdir = 0.0         !Incident direct PAR 400-700 nm (W m-2)
      ecp%CosZen = 0.0          !cos of solar zenith angle

      !Phenology/Growth
      ecp%soiltemp_10d = 0.0d0 !0.7d0! for Hyytiala10-day running avergeage of soil temp (degC)
      ecp%airtemp_10d = 0.0d0 !-3.17!for Hyytiala  !10-day running average of air temp (degC) 
      ecp%paw_10d = 0.5d0!10-day running average of soil moisture (-)
      ecp%par_10d = 100.d0
      ecp%gdd = 0.0d0 !growing degree day
      ecp%fdd = 0.0d0 !falling degree day
      ecp%sgdd = 0.d0
      ecp%sfdd = 0.d0
      ecp%ncd = 0.0d0 !number of chilling day
      ecp%nmd = 0.0d0 !number of mature day (based on air temp)
      ecp%nmsd = 0.0d0 !number of mature day (based on soil temp)
      ecp%daylength(:) = 0.0d0  !day length (min)
      ecp%fall = 1 !KIM- now it's integer...true. !KIM - starting in the winter
      ecp%soiltemp_ltm = 0.d0
      ecp%jday_old = -1
      ecp%jday_tot = 0
      ecp%drought_state = 0

      ecp%C_total = 0.d0
      ecp%C_growth = 0.d0
      ecp%Soilmp(:) = 0.d0

      end subroutine zero_ycell
!**************************************************************************
      subroutine summarize_ycell(ecp)
!@sum Summarize patch properties within a grid cell.
      type(ycelltype) :: ecp
      !-----Local variables------------
      type(patch),pointer :: pp
      integer :: ip             !#patches
      integer :: ia, ib         !counter variable
      real*8 :: fa, laifa, laifasum

      ecp%fv = 0.d0
      call zero_ycell_patchsum(ecp)

      ip = 0
      fa = 0.d0
      laifasum = 0.d0
      pp => ecp%oldest
      do while (ASSOCIATED(pp)) 
        ip = ip + 1

        call summarize_patch(pp) ! make sure patch is summarized, or assume

        fa = fa + pp%area       !For doing wtd avgs
        laifa = pp%area * pp%LAI !For doing wtd avgs
        laifasum = laifasum + laifa

        !- - - - Cohort - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ecp%lai_p = ecp%lai_p + pp%lai_p*pp%area_p
        ecp%ht_p  = ecp%ht_p + pp%ht_p*laifa

        do ia=1,N_COVERTYPES
          ecp%LAI = ecp%LAI + pp%LAIpft(ia) * pp%area !wtd avg by area
          ecp%LAIpft(ia) = ecp%LAIpft(ia) + pp%LAIpft(ia) !xyue not wtd
          ecp%HTpft(ia)  = ecp%HTpft(ia)  + pp%HTpft(ia) !xyue not wtd
          ecp%Phenfpft(ia) = ecp%Phenfpft(ia) + pp%Phenfpft(ia) !xyue not wtd
          ecp%Cdead(ia)  = ecp%Cdead(ia) + pp%Cdead(ia)
          ecp%Clive(ia)  = ecp%Clive(ia) + pp%Clive(ia)
        end do

        ecp%h = ecp%h + pp%h * laifa
        do ia=1,N_DEPTH
          ecp%fracroot(ia) = ecp%fracroot(ia) + pp%fracroot(ia)*pp%area
        end do

        !* IMPORT/EXPORT
        ecp%Ci = ecp%Ci + pp%Ci*pp%area
        ecp%GCANOPY = ecp%GCANOPY + pp%GCANOPY*pp%area

        !* DIAGNOSTICS
        ecp%GPP = ecp%GPP + pp%GPP*pp%area
        ecp%GPP0= ecp%GPP0+ pp%GPP0*pp%area
        do ia=1,N_COVERTYPES
          ecp%GPPpft(ia) = ecp%GPPpft(ia) + pp%GPPpft(ia) * pp%area
          ecp%IPPpft(ia) = ecp%IPPpft(ia) + pp%IPPpft(ia) * pp%area 
          ecp%MTPpft(ia) = ecp%MTPpft(ia) + pp%MTPpft(ia) * pp%area 
        enddo

#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
        ecp%FO3 = ecp%FO3 + pp%FO3*pp%area
        ecp%dFO3= ecp%dFO3+ pp%dFO3*pp%area
#endif
        ecp%IPP = ecp%IPP + pp%IPP*pp%area
        ecp%MTP = ecp%MTP + pp%MTP*pp%area
        ecp%NPP = ecp%NPP + pp%NPP*pp%area
        ecp%resp_r = ecp%resp_r + pp%resp_r*pp%area
        ecp%resp_l = ecp%resp_l + pp%resp_l*pp%area
        ecp%resp_w = ecp%resp_w + pp%resp_w*pp%area
        ecp%resp_p = ecp%resp_p + pp%resp_p*pp%area
        ecp%R_auto = ecp%R_auto + pp%R_auto*pp%area
        ecp%R_root = ecp%R_root + pp%R_root*pp%area  !PK 5/15/07
        ecp%carb_tot  = ecp%carb_tot  + pp%carb_tot*pp%area
        ecp%carb_soil = ecp%carb_soil + pp%carb_soil*pp%area

        !- - - - Patch - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !* IMPORT-PRESCRIBED, EXPORT-SIMULATED
        ecp%z0 = ecp%z0 + pp%z0*pp%area !Area-weighted average
        !* EXPORT
        do ia=1,N_BANDS  !wtd avg by area
          ecp%albedo(ia) = ecp%albedo(ia) + pp%albedo(ia)*pp%area
        end do

        ecp%betad = ecp%betad + pp%betad*pp%area
        do ia=1,N_DEPTH
          ecp%betadl(ia) = ecp%betadl(ia) + pp%betadl(ia)*pp%area
        end do
        ecp%TRANS_SW = ecp%TRANS_SW + pp%TRANS_SW*pp%area !Area-weighted average
        ecp%CO2flux = ecp%CO2flux + pp%CO2flux*pp%area
        
        !* DIAGNOSTICS
        ecp%Soil_resp = ecp%Soil_resp + pp%Soil_resp*pp%area   
        ecp%Tpool(:,:,:) = ecp%Tpool(:,:,:) + pp%Tpool(:,:,:)*pp%area   

        !* IMPORT Variables calculated by GCM/EWB - downscaled from grid cell
        ecp%Soilmoist(:) = ecp%Soilmoist(:) + pp%Soilmoist(:)*pp%area !##

        if ( associated( pp%tallest ) ) ecp%fv = ecp%fv + pp%area

        !* daigs and hacks
        ecp%C_total = ecp%C_total + pp%C_total*pp%area
        ecp%C_growth = ecp%C_growth + pp%C_growth*pp%area
       
        pp => pp%younger
      end do
      
!!!!
      !!!CHECK IF ECP%AREA IS ZERO!
      if (ASSOCIATED(ecp%oldest)) then
        !- - - - Cohort - - - - - - - - - - - - - - - - - - - - - - - - - 
        !laifasum weighting
         if (laifasum.gt.0.) then
            ecp%h = ecp%h/laifasum
         else
            ecp%h = 0.d0
         endif

         ecp%LAI = ecp%LAI/fa
         ecp%lai_p = ecp%lai_p/fa
         ecp%ht_p  = ecp%ht_p/laifasum

         do ia=1,N_COVERTYPES
            ecp%GPPpft(ia) = ecp%GPPpft(ia)/fa
            ecp%IPPpft(ia) = ecp%IPPpft(ia)/fa
            ecp%MTPpft(ia) = ecp%MTPpft(ia)/fa
         enddo

         ecp%fracroot = ecp%fracroot/fa

        !* Flux variables for GCM/EWB - patch total wtd averages
        ecp%Ci = ecp%Ci/fa
        ecp%GCANOPY = ecp%GCANOPY/fa
        ecp%GPP = ecp%GPP/fa
        ecp%GPP0= ecp%GPP0/fa
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
        ecp%FO3 = ecp%FO3/fa
        ecp%dFO3= ecp%dFO3/fa
#endif
C NADINE - IS THIS CORRECT?
        ecp%IPP = ecp%IPP/fa
        ecp%MTP = ecp%MTP/fa
        ecp%NPP = ecp%NPP/fa
        ecp%resp_r = ecp%resp_r/fa
        ecp%resp_l = ecp%resp_l/fa
        ecp%resp_w = ecp%resp_w/fa
        ecp%resp_p = ecp%resp_p/fa
        ecp%R_auto = ecp%R_auto/fa
        ecp%R_root = ecp%R_root/fa  !PK 5/15/07
        ecp%carb_tot  = ecp%carb_tot/fa
        ecp%carb_soil = ecp%carb_soil/fa

        !- - - - - -  Patch-level summary values - PHYSICAL ------------------

        ecp%z0 = ecp%z0/fa !Area-weighted average
        do ia=1,N_BANDS         !Area-weighted average
          ecp%albedo(ia) = ecp%albedo(ia)/fa
        end do 
        ecp%betad = ecp%betad/fa
        ecp%TRANS_SW = ecp%TRANS_SW/fa !Area-weighted average
        ecp%CO2flux = ecp%CO2flux/fa
        ecp%Soil_resp = ecp%Soil_resp/fa  
        ecp%Tpool(:,:,:) = ecp%Tpool(:,:,:)/fa       
        
        !* Variables calculated by GCM/EWB - up/downscaled to/from grid cell
        ecp%Soilmoist(:) = ecp%Soilmoist(:)/fa
        do ia=1,N_DEPTH
          ecp%betadl(ia) = ecp%betadl(ia)/fa
        end do

        ecp%area = fa !## NK used in yibs_diag.
      end if

#ifndef MIXED_CANOPY
      call ycell_update_shc_mosaicveg(ecp)
#else
      call ycell_update_shc(ecp)
#endif

      end subroutine summarize_ycell
!**************************************************************************

      subroutine ycell_update_shc( ecp )
!@sum ycell_update_shc. 
!     Mixed canopies calculation of shc, or generic any veg structure.
!     Old GISS GCM version: shc(avg(lai*vfraction)).
!     Correct version: avg(shc(lai)*vfraction).
      use patches, only : shc_patch
!      use yibs_prescr_veg, only : GISS_shc
      implicit none
      type(ycelltype),pointer :: ecp
      !-----Local---------
      type(patch),pointer :: pp
      real*8 :: shc, pfrac

      shc = 0.d0
      pfrac = 0.d0
      pp => ecp%oldest
      do while (associated(pp))
         shc = shc + shc_patch(pp)*pp%area
         pfrac = pfrac + pp%area
         pp=>pp%younger
      enddo

      if (pfrac>EPS) then 
         shc = shc/pfrac        !Correct way.
      else
         !shc = GISS_shc(0.d0)  !Non-zero shc for zero lai from zero patch area.
         shc = 0.d0
      endif

      ecp%heat_capacity = shc
      end subroutine ycell_update_shc
!******************************************************************

      subroutine ycell_update_shc_mosaicveg( ecp )
!@sum Returns GISS GCM specific heat capacity for ycell.
!     This version preserves old, slightly incorrect lai averaging.
!     shc_ycell = shc(mean lai*vfraction)
!Old ycell_update_shc      
      use yibs_const
      use yibs_pfts, only: COVEROFFSET, alamax, alamin
      use allometryfn, only: do_geo,GISS_shc
      type(ycelltype) :: ecp
      !-----Local---------
      real*8 vfraction(N_COVERTYPES) ! needed for a hack to compute canopy
      type(patch),pointer :: pp      
      real*8 lai, fsum  !lai is mean annual ycell lai.
      integer pft, anum

      lai = 0.d0
      fsum = 0.d0
      vfraction(:) = 0.d0

      !Cover-weighted average
      if (.not.do_geo) then
         !Matthews mean annual LAI
         call ycell_extract_pfts( ecp, vfraction )
         do pft=1,N_PFT
            anum = pft+COVEROFFSET
            lai = lai + .5d0*(alamax(anum) 
     &           + alamin(anum))*vfraction(anum)
            fsum = fsum + vfraction(anum)
         enddo
      endif
      if ( fsum > EPS ) then 
         lai = lai/fsum
      else
         lai = 0.d0
      endif
         
      ecp%heat_capacity=GISS_shc(lai)

      end subroutine ycell_update_shc_mosaicveg

!*************************************************************************
      subroutine init_simple_ycell( ecp, jday, lat, 
     i     vegdata,laidata,hdata,
     i     fracrootdata,soildata,albedodata,soil_texture,
     i     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini,
     i     pdate,hdate,
     i     reinitialize)
!@sum Initializes an ycell assuming one cohort per patch.
      use patches, only : summarize_patch
      type(ycelltype) :: ecp
      integer,intent(in) :: jday
      real*8,intent(in) :: lat
      real*8,intent(in) :: vegdata(N_COVERTYPES) !Veg cover fractions.
      real*8,intent(in) :: laidata(N_COVERTYPES) !LAI
      real*8,intent(in) :: hdata(N_COVERTYPES) !Height
      real*8,intent(in) :: fracrootdata(N_COVERTYPES,N_DEPTH) !Root profile.
      integer,intent(in) :: soildata(N_COVERTYPES)
      real*8,intent(in) :: albedodata(N_BANDS,N_COVERTYPES) !patch, NOTE:snow
      real*8,intent(in) :: soil_texture(N_SOIL_TEXTURES) !soil texture fractions.
      real*8 :: Ci_ini, CNC_ini, Tcan_ini, Qf_ini
      real*8 :: pdate,hdate
      real*8,intent(in) ::
     &      Tpool_ini(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS)  !prescribed soil pools, g/m2 -PK
      logical,intent(in) :: reinitialize
      !-----Local---------
      integer :: ncov, pft
      type(patch),pointer :: pp, pp_tmp, pp_ncov
      real*8 :: sandfrac,clayfrac,smpsat,bch,watsat,watdry

      if ( reinitialize ) then
        ! destroy all existing patches since we are going to 
        ! re-initialize the cell
        pp => ecp%oldest      
        do while ( associated(pp) )
          pp_tmp => pp
          pp => pp%younger
          call patch_destruct(pp_tmp)
        enddo
        nullify(ecp%oldest)
        nullify(ecp%youngest)
      else
        ! just set all areas to 0 since we will reset them
        pp => ecp%oldest      
        do while ( associated(pp) )
          pp%area = 0.d0
          pp => pp%younger
        enddo

      endif

      do ncov=1,N_COVERTYPES           !One patch with one cohort per pft or bare
        pft = ncov - COVEROFFSET
        !### Get from GISS GCM ## vfraction of grid cell and area.

        if (vegdata(ncov)>0.0) then

          call get_patch_by_cover(ecp,ncov,pp_ncov)
          if ( associated(pp_ncov) ) then
            ! if patch is present - just resize it
            pp_ncov%area = vegdata(ncov)
            cycle
          endif
          
         !call insert_patch(ecp,GCMgridareas(j)*vegdata(pnum))
          call insert_patch(ecp,vegdata(ncov),soildata(ncov))
          pp => ecp%youngest

          !## Supply also geometry, clumping index

          ! insert cohort only if population density > 0 (i.e. skip bare soil)
          if ( vegdata(ncov) > EPS ) then 
            if ( pft < 1 .or. pft > N_PFT ) then
              print *,"init_simple_ycell: wrong pft:", pft
              call stop_model("init_simple_ycell: wrong pft",255)
            endif
            call assign_patch(pp,Ci_ini, CNC_ini, pft, Tpool_ini)
            call insert_cohort(pp,pft,hdata(ncov), laidata(ncov),
     &           fracrootdata(ncov,:),
     &           Ci_ini, CNC_ini,0.d0,0.d0,0.d0,0.d0,
     &           0.d0, 1.d0,0.d0,1.d0,
     &           1.d0, 0.d0,
     &           1.d0, -999.d0) !KIM-added for phenology/growth
          endif
          call summarize_patch(pp)

          !CALL CALC_ALBEDO HERE
          pp%albedo = albedodata(:,ncov) !##GISS HACK
        end if
      end do

      ! get rid of patches with 0 area
      pp => ecp%oldest      
      do while ( associated(pp) )
        pp_tmp => pp
        pp => pp%younger
        if( pp_tmp%area == 0.d0 ) call delete_patch(ecp, pp_tmp)
      enddo


      if ( reinitialize ) then
        !Initialize canopy met variables.
        ecp%jday = jday
        ecp%lat = lat
        ecp%TcanopyC = Tcan_ini
        ecp%Qf = Qf_ini
C Nadine crop harvest and plant dates
        ecp%plantdate = pdate
        ecp%harvestdate = hdate

        ! soil textures for CASA
        ecp%soil_texture(:) = soil_texture(:)

      !Soil porosity and wilting? hygroscopic? point for water stress2 calculation. From soilbgc.f.
        sandfrac = soil_texture(1)
        clayfrac = soil_texture(3)

        watsat =  0.489d0 - 0.00126d0*sandfrac !porosity, saturated soil fraction
        smpsat = -10.d0*(10.d0**(1.88d0-0.0131d0*sandfrac))
        bch = 2.91d0 + 0.159d0*clayfrac
        watdry = watsat * (-316230.d0/smpsat) ** (-1.d0/bch)
!      watopt = watsat * (-158490.d0/smpsat) ** (-1.d0/bch)
        ecp%soil_Phi = watsat
        ecp%soil_dry = watdry
      endif

#ifdef OFFLINE
      write(*,*) "soil_Phi, soil_dry",ecp%soil_Phi, ecp%soil_dry
#endif

      call summarize_ycell(ecp)

      

      end subroutine init_simple_ycell

!*************************************************************************
      subroutine assign_ycell( ecp, jday, lat,
     i     soil_texture,    
     i     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini, 
     i     pdate,hdate,
     i     reinitialize)
!@sum assign_ycell. Assigns ycell level values as passed in parameters.
!+    NOTE:  soil_type is a patch-level variable by cover type in Matthews for
!+     the purpose of calculating albedo, whereas soil_texture is an 
!+     ycell-level variable. 

      use patches, only : summarize_patch
      implicit none
      type(ycelltype) :: ecp
      integer,intent(in) :: jday
      real*8,intent(in) :: lat
      real*8,intent(in) :: soil_texture(N_SOIL_TEXTURES)!soil texture fractions.

      real*8 :: Ci_ini, CNC_ini, Tcan_ini, Qf_ini
      real*8 :: pdate,hdate
      real*8,intent(in) ::
     &      Tpool_ini(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS)  
      logical,intent(in) :: reinitialize
      !-----Local---------
      integer :: ncov, pft
      type(patch),pointer :: pp, pp_tmp, pp_ncov
      real*8 :: sandfrac,clayfrac,smpsat,bch,watsat,watdry


      !Assign ycell variables.
      ecp%jday = jday
      ecp%lat = lat
      ecp%TcanopyC = Tcan_ini
      ecp%Qf = Qf_ini

      ! soil textures for CASA
      ecp%soil_texture(:) = soil_texture(:)

C Nadine crop harvest and plant dates
      ecp%plantdate = pdate
      ecp%harvestdate = hdate

      !Soil porosity and wilting? hygroscopic? point for water stress2 calculation. From soilbgc.f.
! is it supposed to be
      sandfrac = soil_texture(1)
      clayfrac = soil_texture(3)

      watsat =  0.489d0 - 0.00126d0*sandfrac !porosity, saturated soil fraction
      smpsat = -10.d0*(10.d0**(1.88d0-0.0131d0*sandfrac))
      bch = 2.91d0 + 0.159d0*clayfrac
      watdry = watsat * (-316230.d0/smpsat) ** (-1.d0/bch)
!      watopt = watsat * (-158490.d0/smpsat) ** (-1.d0/bch)

      ecp%soil_Phi = watsat
      ecp%soil_dry = watdry


#ifdef OFFLINE
      write(*,*) "soil_Phi, soil_dry",ecp%soil_Phi, ecp%soil_dry
#endif

      call summarize_ycell(ecp)

      !print *,"leaving assign_ycell"
      !call ycell_print(6,ecp)

      end subroutine assign_ycell

  !*********************************************************************
      subroutine assign_ycell_soilcarbon(ecp,Tpools)
!@sum Distribute ycell grid soil carbon (e.g. from file) to subgrid pathes.
      use patches, only : assign_patch, summarize_patch
      implicit none
      type(ycelltype) :: ecp
      real*8,intent(in) ::
     &      Tpools(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS)  
      !--- Local ---
      type(patch),pointer :: pp

      pp => ecp%oldest      
      do while ( associated(pp) )
         pp%Tpool(NITROGEN,(NLIVE+1):NPOOLS,:) = 0.d0 !## HACK WHILE NO PROGNOSTIC N #
            if (ASSOCIATED(pp%tallest)) then
                call assign_patch_soilcarbon(pp,pp%tallest%pft,Tpools)
                call summarize_patch(pp) 
            endif
          pp => pp%younger
      enddo

      end subroutine assign_ycell_soilcarbon
  !*********************************************************************

      subroutine ycell_construct(ecp)
!@sum ycell_construct  Allocate memory for an ycell and nullify pointers.
      implicit none
      type(ycelltype), pointer :: ecp

      ! allocate memory
      allocate( ecp )
      allocate( ecp%LAIpft(N_COVERTYPES) )
      allocate( ecp%HTpft(N_COVERTYPES) )
      allocate( ecp%GPPpft(N_COVERTYPES) )
      allocate( ecp%IPPpft(N_COVERTYPES) )
      allocate( ecp%MTPpft(N_COVERTYPES) )
      allocate( ecp%Cdead(N_COVERTYPES) )
      allocate( ecp%Clive(N_COVERTYPES) )
      allocate( ecp%Phenfpft(N_COVERTYPES) )
      allocate( ecp%Soilmp(N_DEPTH) )
      allocate( ecp%fice(N_DEPTH) )
      allocate( ecp%fracroot(N_DEPTH) )
      allocate( ecp%betadl(N_DEPTH) )

      ! set pointers
      nullify( ecp%youngest )
      nullify( ecp%oldest   )

      ! for now set all values to zero or defaults
      call zero_ycell(ecp)

      ! maybe the following is not needed, but I guess we wanted
      ! to initialize the cells by default to bare soil (sand)
      call insert_patch(ecp,1.d0, soil_color_prescribed(SAND))

      end subroutine ycell_construct

 !*********************************************************************

      subroutine ycell_destruct(ecp)
!@sum ycell_destruct  Deallocate memory for an ycell.
      implicit none
      type(ycelltype), pointer :: ecp
      !---
      type(patch), pointer :: pp, pp_tmp

      ! destroy other patches
      pp => ecp%oldest      
      do while ( associated(pp) )
        pp_tmp => pp
        pp => pp%younger
        call patch_destruct(pp_tmp)
      enddo

      ! deallocate memory
      deallocate( ecp )
      nullify( ecp )

      end subroutine ycell_destruct

 !*********************************************************************
      

      subroutine ycell_extract_pfts(ecp, vfraction)
!@sum ycell_extract_pfts Extract cover fraction of subgrid patches.
      type(ycelltype) :: ecp
      real*8 :: vfraction(:)
      !---
      type(patch), pointer :: pp
      real*8 :: vfraction_patch(size(vfraction))

      vfraction(:) = 0.d0
      pp => ecp%oldest
      do while( associated(pp) )
        call patch_extract_pfts(pp, vfraction_patch)
        vfraction(:) = vfraction(:) + vfraction_patch(:)*pp%area
        pp => pp%younger
      enddo

      end subroutine ycell_extract_pfts

 !*********************************************************************

      subroutine get_patch_by_cover(ecp, ncov, pp_ncov)
!@sum Point to a patch of a given cover type in the given ycell.
      type(ycelltype) :: ecp
      integer :: ncov
      type(patch), pointer :: pp_ncov
      !---
      type(patch), pointer :: pp

      !write(0,*) "entered get_patch_by_cover", ncov
      nullify( pp_ncov )
      pp => ecp%oldest
      do while( associated(pp) )
        !write(0,*) "inside loop"
        if ( associated(pp%tallest) ) then
          !write(0,*) pp%tallest%pft+COVEROFFSET
          if ( pp%tallest%pft == ncov - COVEROFFSET ) then
            pp_ncov => pp
            return
          endif
        else
          if ( ncov == COVER_SAND .and. pp%soil_type == 1 ) then
            !write(0,*) 1
            pp_ncov => pp
            return
          else if ( ncov == COVER_DIRT .and. pp%soil_type == 2 ) then
            !write(0,*) 10
            pp_ncov => pp
            return
          endif
        endif
        pp => pp%younger
      enddo
      
      !write(0,*) "null"

      end subroutine get_patch_by_cover

!----------------------------------------------------------------------
      real*8 function ycell_carbon(ecp, 
     &     ecp_Cfol,ecp_Cstem,ecp_Croot, ecp_Csoil)
!@sum ycell_carbon (kg-C m-2). Return total carbon in ycell per land
!@+       area (exclude any water area in ycell grid).  
!@+       Optionally return component carbon pools.
      type(ycelltype) :: ecp
      real*8, optional :: ecp_Cfol
      real*8, optional :: ecp_Cstem
      real*8, optional :: ecp_Croot
      real*8, optional :: ecp_Csoil
      !--Local----
      type(patch), pointer :: pp
      real*8 :: ecp_kgCm2
      real*8 :: pp_Cfol, pp_Cstem, pp_Croot, pp_Csoil
      real*8 :: sumarea

      if (present(ecp_Cfol)) ecp_Cfol = 0.d0
      if (present(ecp_Cstem)) ecp_Cstem = 0.d0
      if (present(ecp_Croot)) ecp_Croot = 0.d0
      if (present(ecp_Csoil)) ecp_Csoil = 0.d0

      ecp_kgCm2 = 0.d0
      pp_Cfol = 0.d0
      pp_Cstem = 0.d0
      pp_Croot = 0.d0
      pp_Csoil = 0.d0
      sumarea = 0.d0

      pp => ecp%oldest

      do while (associated(pp))
         sumarea = sumarea + pp%area 
         ecp_kgCm2 = ecp_kgCm2 + pp%area*patch_carbon(pp,
     &        pp_Cfol, pp_Cstem, pp_Croot, pp_Csoil)

         if (present(ecp_Cfol)) ecp_Cfol = ecp_Cfol + pp%area*pp_Cfol
         if (present(ecp_Cstem)) ecp_Cstem = ecp_Cstem +pp%area*pp_Cstem
         if (present(ecp_Croot)) ecp_Croot = ecp_Croot +pp%area*pp_Croot
         if (present(ecp_Csoil)) ecp_Csoil = ecp_Csoil +pp%area*pp_Csoil

         pp => pp%younger
      end do

      if (sumarea.eq.0.d0) then
         ycell_carbon = 0.d0
      else
         ycell_carbon = ecp_kgCm2/sumarea
         if (present(ecp_Cfol)) ecp_Cfol = ecp_Cfol/sumarea
         if (present(ecp_Cstem)) ecp_Cstem = ecp_Cstem/sumarea
         if (present(ecp_Croot)) ecp_Croot = ecp_Croot/sumarea
         if (present(ecp_Csoil)) ecp_Csoil = ecp_Csoil/sumarea
      endif

      end function ycell_carbon
!----------------------------------------------------------------------

      end module ycells
