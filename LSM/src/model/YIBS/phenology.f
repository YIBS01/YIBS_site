      module phenology
!@sum Routines to calculate phenological change in an ycell:
!@sum budburst/leafout, albedo change, senescence
!@auth Y. Kim
#ifdef YIBS_1D_DIAG
#define PHENOLOGY_DIAG
!#define DEBUG
#endif

      use yibs_types
      use yibs_const
      use yibs_pfts
      use allometryfn

      implicit none
!      public veg_init
      public clim_stats
      public pheno_update
      public litter_cohort, litter_patch   !Now called from veg_update

      private pheno_update_coldwoody
      private pheno_update_coldherbaceous
      private pheno_update_drought
      private photosyn_acclim

      !**********************************************************************

      !** GROWTH MODEL CONSTANTS - phenology & carbon allocation 
      !*l_fract: fraction of leaf C retained after leaf fall (unitless) (value from ED) 
      real*8, parameter :: l_fract = 0.50d0 
      real*8, parameter :: C2B = 2.0d0 
      !*temperature constraint for cold-deciduous PFTs (Botta et al. 1997)
      !*airtemp_par !base temperature to calculate the growing degree days (gdd)
      !*gdd_par1/2/3: paramters to estimate the threshold for gdd
      !*gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)    
      real*8, parameter :: airtemp_par  = 5.d0 
      real*8, parameter :: airtemp_fall = 20.d0 
      real*8, parameter :: gdd_par1 = -110.d0 
      real*8, parameter :: gdd_par2 = 550.d0   
      real*8, parameter :: gdd_length = 380.d0 
      real*8, parameter :: fdd_threshold = -140.d0 
      real*8, parameter :: fdd_length = 410.d0 
      real*8, parameter :: gdd_par3 = -0.01d0 
      real*8, parameter :: soiltemp_fall = 10.d0 
      real*8, parameter :: sfdd_threshold = -80.d0 
      real*8, parameter :: sfdd_length = 100.d0 
      real*8, parameter :: soilt_base = -5.d0 !Added from YK
      !*sgdd_threshold & length - tuning parameters (tuned for Barrow)
      real*8, parameter :: sgdd_threshold = 100.d0
      real*8, parameter :: soiltemp_par = 0.d0
      real*8, parameter :: sgdd_length=100.d0   
      !*ld_max/ld_min (minute): light length constraint for cold-deciduous woody PFTs 
      real*8, parameter :: ld_min = 585.d0   
      real*8, parameter :: ld_max = 695.d0 
      ! For crop grow length (days)
      real*8, parameter :: dgrow  = 30.d0 
      ! For crop fall length (days)
      real*8, parameter :: dfall  = 15.d0 
      real*8, parameter :: ddfacu = 1.d0/15.d0
      !_w for woody & _h  (tunned for MMSF) - max/min/resistance parameters for woody
      real*8, parameter :: betad_max_w = 1.0d0    ! modify_04.28: 0.1 --> 0.5
      real*8, parameter :: betad_min_w = 0.4d0 
      real*8, parameter :: betad_res_w = 1.0d0    ! modify_05.26: 0.25 -> 1.0
      ! _h for herbaceous  (tunned for Vaira/Tonzi) - max/min/resistance prameters for woody
      real*8, parameter :: betad_max_h = 0.9d0 
      real*8, parameter :: betad_min_h = 0.3d0
      real*8, parameter :: betad_res_h = 1.0d0      
      !*light-controll phenology model (Kim et al; originally for ED2)
      !*different from the original implementation, as it cannot be directly implemented due to model difference.
      !*e.g.) PAR instead of Rshort & other differences in parameterization requires the model to be tunned for YIBS.
      !*par_turnover_int & par_turnover_slope  - tunning parameters (tunned for TNF; not finalized)
      real*8, parameter :: par_turnover_int = -12.d0 
      real*8, parameter :: par_turnover_slope = 0.18d0  
      !*r_fract: fraction of excess c going to seed reproduction (value from ED)
      real*8, parameter :: r_fract = 0.3d0
      !*c_fract: fraction of excess c going to clonal reproduction - only for herbaceous (value from ED)
      real*8, parameter :: c_fract = 0.7d0
      !*mort_seedling: mortality rate for seedling (value from ED)
      real*8, parameter :: mort_seedling = 0.90d0 

      contains


      !*********************************************************************
      subroutine clim_stats(dtsec, ecp, config,dailyupdate)
!@sum Calculate climate statistics such as 10 day running mean 
!@+   Called by ent.f every physical time step.
      use soilbgc, only : Soillayer_convert_YIBS 
      use yibs_prognostic, only: calc_daylength
      real*8,intent(in) :: dtsec           !dt in seconds
      type(ycelltype) :: ecp      
      type(yibs_config) :: config
      logical, intent(in) :: dailyupdate  
      !-----local--------
      type(patch), pointer :: pp
      type(cohort), pointer :: cop 
      !*local variables for ycell-level envrionment variable
      real*8 :: airtemp        !air temperature degC 
      real*8 :: soiltemp       !soil temperature degC
      real*8 :: par            !photosynthetic active radiation (PAR)
      !*local variables for ycell-level variables, 
      !*updated in this subroutine
      real*8 :: airtemp_10d    !10 day running mean of air temperature
      real*8 :: soiltemp_10d   !10 day running mean of soil temperature
      real*8 :: soiltemp_ltm   !long-term mean soil temperature
      real*8 :: ddt            !number of total days
      real*8 :: par_10d        !10 day running mean of PAR
      real*8 :: gdd            !growing degree days, based on air temperature
      real*8 :: ncd            !number of chilling days, based on air temperature
      real*8 :: sgdd           !growing degree days, based on soil temperature
Cxyue
      real*8 :: fdd            !falling degree days, based on air temperature
      real*8 :: sfdd           !falling degree days, based on soil temperature
Cxyue
      !*PAR-limited phenology parameters
      real*8 :: par_crit       !PAR threshold
      logical :: par_limit     !logical whether PAR-limited phenology parameterization is applied or not for certain PFTs
      real*8 :: turnover0      !turnover amplitude, calculated with the phenology parameterization
      real*8 ::  llspan0       !leaf life span, calculated with the phenology parameterization
      !*soil temperature for CASA layers
      real*8 :: Soiltemp2layer(N_CASA_LAYERS)  
     
  
      airtemp = ecp%TairC
      call Soillayer_convert_YIBS(ecp%Soiltemp(:), SOILDEPTH_m, 
     &     Soiltemp2layer) 
      soiltemp = Soiltemp2layer(1)

      soiltemp_10d = ecp%soiltemp_10d
      airtemp_10d = ecp%airtemp_10d
      par_10d = ecp%par_10d  
      gdd = ecp%gdd
      ncd = ecp%ncd
      sgdd = ecp%sgdd
      fdd  = ecp%fdd
      sfdd = ecp%sfdd
      soiltemp_ltm = ecp%soiltemp_ltm

      !*10-day running mean of Air Temperature
      airtemp_10d = running_mean(dtsec, 10.d0, airtemp, airtemp_10d)

      !*10-day running mean of Soil Temperature
      soiltemp_10d = running_mean(dtsec, 10.d0, soiltemp, soiltemp_10d)

      !*Long-term running mean of Soil Temperature
      If (ecp%jday_old .ne. ecp%jday) then
          ecp%jday_tot = ecp%jday_tot+1
          ecp%jday_old = ecp%jday
          If (mod(ecp%jday_tot,365) .eq. 0) then        ! check drought state annually
              If (soiltemp_ltm .gt. 12.d0) then
                  ecp%drought_state = 1
              Else
                  ecp%drought_state = 0
              Endif
          Endif
      Endif
      ddt = max(ecp%jday_tot*1.0d0, 1.0d0)
      soiltemp_ltm = running_mean(dtsec, ddt, soiltemp, soiltemp_ltm)

      !*10-day running mean of PAR
      par = ecp%IPARdif + ecp%IPARdir !total PAR is the sum of diffused and direct PARs
      par_10d = running_mean(dtsec, 10.d0, par, par_10d)

      !*daylength
Cxyue
c      if (ecp%CosZen > 0.d0) then
c         ecp%daylength(2) = ecp%daylength(2) + dtsec/60.d0
c      end if
      ecp%daylength(2) = calc_daylength(ecp%jday, ecp%lat)
c      write(205,*) ecp%jday, ecp%lat, ecp%daylength(2)
Cxyue

      !*GDD & NCD - Update Once a day 
      if (dailyupdate) then
         !*Calculate Growing degree days
         if (airtemp_10d .ge. airtemp_par) then
            gdd = gdd + ( airtemp_10d - airtemp_par )
         end if
         !*Calculate Growing degree days for soil temperature 
         if (soiltemp_10d .ge. soiltemp_par) then
            sgdd = sgdd + ( soiltemp_10d -soiltemp_par )
         end if
         !*Number of chilling days 
         !number of days below airtemp_par - chilling requirements
         if (airtemp_10d .lt. airtemp_par) then
            ncd = ncd +  1.d0
         end if
Cxyue
         !*Calculate Falling degree days
         if (airtemp_10d .lt. airtemp_fall) then
            fdd = fdd + ( airtemp_10d - airtemp_fall )
         end if
         !*Calculate Falling degree days
         if (soiltemp_10d .lt. soiltemp_fall) then
            sfdd = sfdd + ( soiltemp_10d - soiltemp_fall )
         end if
Cxyue
         !*If the season is fall or not (if fall, it's 1; else, it's 0)
         !1) it is to control the phenological status      
         !2) it is determined according to whether the daylength is decreasing (i.e., fall) or not. 
         if (NInt(ecp%daylength(2)) .lt. NInt(ecp%daylength(1)) ) then
            ecp%fall = 1
         else if (NInt(ecp%daylength(2)).gt.NInt(ecp%daylength(1))) then
            ecp%fall = 0
         end if 

Cxyue for polar regions, daylength is always 24 hours during summer
         If ( ecp%lat .ge. 0 ) then
            If (ecp%jday .ge. 356 .or. ecp%jday .lt. 173) then
               ecp%fall = 0
            else
               ecp%fall = 1
            Endif
         else
            If (ecp%jday .ge. 356 .or. ecp%jday .lt. 173) then
               ecp%fall = 1
            else
               ecp%fall = 0
            Endif
         Endif
            
       end if

      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 

        cop => pp%tallest
        do while(ASSOCIATED(cop))
          
          !*10-day running mean of stressH2O (betad)
          cop%betad_10d = running_mean(dtsec, 10.d0, 
     &                    cop%stressH2O, cop%betad_10d)
     &                     
          !*Daily carbon balance
          !it is used for the carbon allocation
          cop%CB_d =  cop%CB_d + cop%NPP*dtsec/cop%n*1000.d0

          !*********************************************
          !* evergreen broadleaf - PAR limited - not finalized yet!
          !********************************************* 

          !if it is evergreen & broadleaf, radiation-limited phenology is working    
          par_limit = ((pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
     &                (pfpar(cop%pft)%leaftype.eq.BROADLEAF))
          !par_limit = .false. !temp. suppress
          !raidation-limited phenology model 
          if (par_limit) then
             par_crit = - par_turnover_int/par_turnover_slope 

             !calculate the turnover amplitude 
             !(relative ratio of turnover compared to its intrinsic turnover rate)
             !based on PAR
             turnover0 = min(100.d0, max(0.01d0, 
     &          par_turnover_slope*par_10d + par_turnover_int))

             if (par_10d .lt. par_crit) turnover0 = 0.01d0

             !calculate 10 day running mean of turnover amplitude
             cop%turnover_amp = running_mean(dtsec,10.d0, 
     &                          turnover0, cop%turnover_amp)

             !calculate the leaf life span based on turnover amplitude
             !lrage is in year, llspan is in month, and then 12 is used to convert the units
             llspan0 = pfpar(cop%pft)%lrage*12.d0/cop%turnover_amp 

             !calculate 90 day running mean of llspan
             cop%llspan = running_mean(dtsec, 90.d0,llspan0,cop%llspan)

          else

             cop%turnover_amp = 1.d0 
             cop%llspan = undef !-999.d0

          endif
          

          !**************************************************************
          !* Update photosynthetic acclimation factor for evergreen veg
          !**************************************************************
          if (config%do_frost_hardiness) then 
             if (((pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
     &           (pfpar(cop%pft)%leaftype.eq.NEEDLELEAF)).or.
     &           (pfpar(cop%pft)%phenotype.eq.COLDDECID).or.
     &           (pfpar(cop%pft)%phenotype.eq.COLDDROUGHTDECID)) then
                call photosyn_acclim(dtsec,airtemp_10d,cop%Sacclim) 
             else
                cop%Sacclim = 25.d0 !Force no cold hardening, mild temperature.
             endif
	  else
              cop%Sacclim = 25.d0 !Force no cold hardening, mild temperature.
          endif

          cop => cop%shorter  
        end do        
        pp => pp%younger 
      end do 

      ecp%soiltemp_10d = soiltemp_10d
      ecp%soiltemp_ltm = soiltemp_ltm
      ecp%airtemp_10d = airtemp_10d
      ecp%par_10d = par_10d
      ecp%gdd = gdd
      ecp%ncd = ncd
      ecp%sgdd = sgdd
Cxyue
      ecp%fdd  = fdd
      ecp%sfdd = sfdd
Cxyue

      end subroutine clim_stats

      !*********************************************************************   
      subroutine photosyn_acclim(dtsec,Ta,Sacc)
!@sum Calculate Sacclim accumulator state used to calculate 
!@+   acclimation/frost hardiness for boreal coniferous forests.  
!@+   Called by clim_stats during update of phenology climate stats.
!@+   Based on Repo et al (1990), 
!@+   Hanninen & Kramer (2007),and  Makela et al (2006)
      implicit none
      real*8,intent(in) :: dtsec ! time step size [sec]
      real*8,intent(in) :: Ta ! air temperature [deg C]
      real*8,intent(inout) :: Sacc ! state of acclimation [deg C]

      !----Local-----
      real*8,parameter :: tau_inv = 2.22222e-6 
                          ! inverse of time constant of delayed
                          ! response to ambient temperature [sec] = 125 hr 
                          ! Makela et al (2004) for Scots pine

!      Use a first-order Euler scheme
       Sacc = Sacc + dtsec*(tau_inv)*(Ta - Sacc) 

!     Predictor-corrector method requires temperature from next timestep
!       Sacc_old = Sacc
!       Sacc = Sacc_old + dtsec*(1/tau_acclim)*(Ta - Sacc ) 
!       Sacc = Sacc_old + ((1/tau_acclim)*(Ta - Sacc_old)+
!     &                     (1/tau_acclim)*(Ta_next - Sacc))*0.5d*dtsec

      end subroutine photosyn_acclim


      !*********************************************************************   
      subroutine pheno_update(pp)
!@sum Update statstics for phneology_update    
!@sum DAILY time step.
      use yibs_const

      type(patch) :: pp
      !--Local-----
      type(cohort), pointer :: cop
      integer :: pft
      integer :: phenotype
      integer :: dstate
      real*8 :: airtemp_10d    !10 day running mean of air temperature
      real*8 :: soiltemp_10d   !10 day running mean of soil temperature    
      real*8::  betad_10d      !10 day running mean of betad (calculated with stressH2O)  
      real*8 :: ld             !day length in minutes
      real*8 :: gdd            !growing degree days, based on air temperature
      real*8 :: ncd            !number of chilling days, based on air temperature
      real*8 :: sgdd           !growing degree days, based on soil temperature 
Cxyue
      real*8 :: jday           !julian day
      real*8 :: pday           !crop plant day
      real*8 :: hday           !crop harvest day
      real*8 :: dday           !crop growing day from plant to harvest
      real*8 :: fdd            !falling degree days, based on air temperature
      real*8 :: nmd            !number of mature days, based on air temperature
      real*8 :: nmsd           !number of mature days, based on soil temperature
      real*8 :: sfdd           !falling degree days, based on soil temperature
      real*8 :: phenostatus_c
      real*8 :: phenostatus_d
Cxyue
      real*8 :: airt_adj       !adjustment for air temperature threshold (airt_max & airt_min)
      real*8 :: soilt_adj      !adjustment for soil temperature threshold (soilt_max & soilt_min)
      !*phenofactor : phenological elongation factor, ranging 0 (no leaf) to 1 (full leaf)
      !_c for the cold-deciduous & _d for the drought-deciduous
      real*8 :: phenofactor
      real*8 :: phenofactor_c
      real*8 :: phenofactor_d
      !*phenostatus :  define phenological status
      !1 - no leaf (phenofactor, equal to 0; after leaf-off in the  fall, until leaf green-up in the next spring)
      !2 - growing leaf (phenofactor, increasing from 0 to 1 in the spring) 
      !3 - leaf in full growth (phenofactor, equal to 1)
      !4 - leaf senescence (phenofactor, decreasing from 1 to 0 in the fall)
      real*8 :: phenostatus
      !*phenogy type 
      logical :: cold_limit    !.true.=cold deciduous
      logical :: drought_limit  !.true.=drought deciduous
      !*GDD treshold (White et al.)      
      real*8 :: gdd_threshold
      !*whther the season is fall or not (determined in clim_stats)
      logical :: fall         !.true. for fall
      logical :: fall_old     !.true. for fall (the day before)
      !*whether  PFT is wood or not
      logical :: woody
      logical :: iscrop
      !*X day to mature: mature=1+X/1000 
      !1.1 means 100 days to mature.
      !it is devised to prevent abrupt grass green-up 
      !in the middle of winter due to a couple of warm days. 
      real*8 :: mature = 1.1d0 
 
      soiltemp_10d = pp%cellptr%soiltemp_10d
      airtemp_10d = pp%cellptr%airtemp_10d
      gdd = pp%cellptr%gdd
      ncd = pp%cellptr%ncd
      if ( pp%cellptr%fall == 1 ) then
        fall = .true.
      else
        fall = .false.
      endif
      sgdd = pp%cellptr%sgdd
Cxyue
      jday = pp%cellptr%jday+0.0d0
      pday = pp%cellptr%plantdate
      hday = pp%cellptr%harvestdate
      dday = hday - pday
      if (dday .le. 0.d0) dday = dday + 365.d0
      if ( pp%cellptr%fall_old == 1 ) then
        fall_old = .true.
      else
        fall_old = .false.
      endif
      sfdd = pp%cellptr%sfdd
      fdd  = pp%cellptr%fdd
      nmd  = pp%cellptr%nmd
      nmsd = pp%cellptr%nmsd
      dstate = pp%cellptr%drought_state
Cxyue
      ld = pp%cellptr%daylength(2)

      cop => pp%tallest
      do while(ASSOCIATED(cop))                
         phenofactor_c=cop%phenofactor_c
         phenofactor_d=cop%phenofactor_d
         phenofactor=cop%phenofactor
         cop%phen_old=max(phenofactor,0.01d0)         ! xyue for active growth
         phenostatus_c=cop%phenostatus_c
         phenostatus_d=cop%phenostatus_d
         betad_10d=cop%betad_10d
         pft=cop%pft
         phenotype=pfpar(pft)%phenotype
         woody = pfpar(pft)%woody

         !***********************************************
         !*Determine whther PFT cold or drought deciduous
         !***********************************************
         if (phenotype .eq. COLDDECID) then 
            cold_limit = .true.
            drought_limit = .false.
         else if (phenotype .eq. DROUGHTDECID) then 
            cold_limit = .true.
            drought_limit = .true.
         else if (phenotype .eq. EVERGREEN) then
            cold_limit = .false.
            drought_limit = .false.
         else !any of cold and drought deciduous
            cold_limit = .true.
            drought_limit = .true.            
         end if

         !*Set the air/temperature adjustment for specific PFTs
         airt_adj=0.d0
         if (pft.eq.DROUGHTDECIDBROAD)airt_adj=10.d0
         soilt_adj=0.d0
         if (pft .eq. GRASSC3ARCTIC)soilt_adj=-5.d0      
         if (pft .eq. CROPS) then
             iscrop = .true.
         else
             iscrop = .false.
         endif

         !*******************************************
         !*Update the phenology for Cold-deciduous
         !*******************************************
         if (cold_limit)then

           !*Cold-deciduous Woody
           if (woody) then  
             !*Update phenofactor and phenostatus
             call pheno_update_coldwoody(cop%phenostatus, 
     i            fall, fall_old, airtemp_10d, airt_adj, ld,
     o            gdd, ncd, fdd, nmd,
     o            phenofactor_c, phenostatus_c)

          !*Cold-decidous Herbaceous
          else 
             !*Update phenofactor and phenostatus
             call pheno_update_coldherbaceous(cop%phenostatus,
     i            fall, fall_old, soiltemp_10d, soilt_adj, 
     o            sgdd, sfdd, nmsd, ld, 
     o            phenofactor_c, phenostatus_c)
           end if

         else !cold_limit=.false.
           phenofactor_c = 1.d0
           phenostatus_c  = 3.d0
         end if          
                 
         !*******************************************
         !*Update the phenology for Drought-deciduous
         !*******************************************
         if (drought_limit)then

           !*Drought-deciduous Woody
           if (woody) then 
             !*Update phenofactor and phenostatus
             call pheno_update_drought(cop%phenostatus, 
     i            mature, 
     i            betad_10d, betad_min_w, betad_max_w, betad_res_w, 
     o            phenofactor_d, phenostatus_d)

         !*Drought-decidous Herbaceous
          else 
             !*Update phenofactor and phenostatus
             call pheno_update_drought(cop%phenostatus, 
     i            mature, 
     i            betad_10d, betad_min_h, betad_max_h, betad_res_h, 
     o            phenofactor_d, phenostatus_d)
           end if   
    
         else !drought_limit=.false.
           phenofactor_d = 1.d0
           phenostatus_d = 3.d0
         end if  

         !*******************************************
         !*Update the phenology according to PFTs 
         !*******************************************
         if (phenotype .eq. EVERGREEN) then !leaf in full growth
            phenofactor = 1.d0
            phenostatus = 3.d0
         elseif (phenotype .eq. COLDDECID .or. dstate .eq. 0) then
            phenofactor = phenofactor_c
            phenostatus = phenostatus_c
         elseif (woody .and. dstate .eq. 1) then 
            phenofactor = phenofactor_d
            phenostatus = phenostatus_d
         elseif (phenotype .eq. DROUGHTDECIDBROAD .and. 
     &            .not.phenostatus_c.lt.3.d0) then
             phenofactor = phenofactor_c
             phenostatus = phenostatus_c
Cxyue
         else !any of cold and drought deciduous

           phenofactor = min(phenofactor_c, phenofactor_d)
           if(phenofactor .le. 0.01d0) then
              phenostatus = 1.d0
           elseif (phenofactor .ge. 0.95d0) then
              phenostatus = 3.d0
           elseif (fall) then
              phenostatus = 4.d0
           else
              phenostatus = 2.d0
           endif
Cxyue
         end if
         
         if (iscrop .and. dday .gt. dgrow+dfall) then

         call pheno_update_crop(jday, pday, hday, 
     &        phenofactor, phenostatus)
   
         endif
        
         !*increment phenostatus by 0.001 
         !to track how many days after phenostatus has been changed
Cxyue
c         if (aint(cop%phenostatus).eq.aint(phenostatus))then
c            phenostatus = cop%phenostatus + 1.d0/1000.d0  
c         end if
Cxyue
    
#ifdef DEBUG
             write(202,*) pp%cellptr%fall
     &      ,phenofactor
     &      ,phenofactor_c,phenofactor_d
     &      ,phenostatus 
     &      ,phenostatus_c, phenostatus_d
     &      ,betad_10d,soiltemp_10d
     &      ,gdd, sgdd, fdd, sfdd, ncd, nmd, nmsd 
     &      ,dstate, pp%cellptr%soiltemp_ltm
#endif
         
         cop%phenofactor_c=max(phenofactor_c,0.01d0)
         cop%phenofactor_d=max(phenofactor_d,0.01d0)
         cop%phenofactor=max(phenofactor,0.01d0)
         cop%phenostatus=phenostatus
         cop%phenostatus_c=phenostatus_c
         cop%phenostatus_d=phenostatus_d
   
         cop => cop%shorter 
      
      end do   

      pp%cellptr%gdd = gdd
      pp%cellptr%ncd = ncd
      pp%cellptr%sgdd = sgdd
      pp%cellptr%fdd  = fdd
      pp%cellptr%nmd  = nmd
      pp%cellptr%nmsd = nmsd
      pp%cellptr%sfdd = sfdd
      
      end subroutine pheno_update
      !*********************************************************************  
      subroutine pheno_update_coldwoody(phenostatus, 
     i            fall, fall_old, airtemp_10d, airt_adj, ld,
     o            gdd, ncd, fdd, nmd, 
     o            phenofactor_c, phenostatus_c)
!@sum Update phenology for cold-decidous woody PFTs
!@sum Called from pheno_update

      use yibs_const

      !input variables
      real*8, intent(in) :: phenostatus !phenological status (refer pheno_update for details) 
      logical,intent(in) :: fall        !.true. if the season is fall
      logical,intent(in) :: fall_old     !.true. if the season is fall
      real*8, intent(in) :: airtemp_10d !10 day running mean of air temperature
      real*8, intent(in) :: airt_adj    !adjustment for air temperature threshold (airt_max & airt_min)
      real*8, intent(in) :: ld          !day length in minutes
      !in/output variables
      real*8, intent(inout) :: gdd      !growing degree days, based on air temperature
      real*8, intent(inout) :: ncd      !number of chilling days, based on air temperature  
Cxyue
      real*8, intent(inout) :: fdd      !falling degree days, based on air temperature
      real*8, intent(inout) :: nmd      !number of mature days, based on phenostatus  
Cxyue
      !output variables
      real*8, intent(out) :: phenofactor_c !phenological factor for cold deciduous (refer pheno_update for details) 
      real*8, intent(inout) :: phenostatus_c !phenological status for cold deciduous (refer pheno_update for details)
      !local variables
      real*8 :: gdd_threshold   !GDD treshold (White et al.)      
      real*8 :: a

      !*GDD threshold for leaf green-up 
      gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)  

      !*FDD threshold for leaf off
      if (abs(phenostatus_c-3.d0) .le. 0.01d0)then
          nmd = nmd + 1.0d0
      endif

      if (.not. fall .and. fall_old ) then
          if (phenostatus_c.ge.3.d0) then
             gdd = phenofactor_c*gdd_length+gdd_threshold
             fdd = 0.d0
             phenostatus_c = 2.d0
          else
             gdd = 0.d0
             ncd = 0.d0
          endif
      endif
      if (.not. fall) fdd = 0.d0

      !*Leaf-on in the spring, triggered by thermal sum
      !if gdd is larger than its threshold 
      !in the spring (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X)), 
      !determine the phenofactor and corresponding phenostatus.
       if ((phenostatus_c.lt.3.d0).and.(gdd.gt.gdd_threshold))then  

         !determine phenofactor by scaling gdd with gdd_threshold and gdd_length
         phenofactor_c = min (1.d0,(gdd-gdd_threshold)/gdd_length)
         if (phenofactor_c .lt. 1.d0) then
            phenostatus_c = 2.d0 !growing leaf
         else 
            phenostatus_c = 3.d0 !leaf in full grwoth
            nmd = 0.d0           !reinitialize nmd for falling
            fdd = 0.d0           !reinitialize fdd for falling
         end if
         return

      end if

      !*Leaf-off in the fall 
      !*Leaf-off triggered by air temperature
      if ((phenostatus_c.ge.3.d0).and.(fdd.lt.fdd_threshold))then

         phenofactor_c = min(phenofactor_c,max(0.d0,
     &      1.0d0+(fdd-fdd_threshold)/fdd_length))

         if (phenofactor_c .le. 0.01d0) then
            phenostatus_c = 1.d0   !no leaf
            ncd = 0.d0             !zero-out ncd once complete leaf-off occurs.
            gdd = 0.d0             !zero-out gdd once complete leaf-off occurs.
         else  
            phenostatus_c = 4.d0   !leaf senescence
         end if 

      end if
      a = max(0.d0, min(1.0d0+(fdd-fdd_threshold)/fdd_length,1.d0))

      !*Leaf-off triggered by day-length
      !if day length is falling shorter than its maximum 
      if ((phenostatus_c.ge.3.d0).and.(ld.lt.ld_max)) then

         !dtermine phenofactor by scaling the day length with its min and max.
         phenofactor_c = min(phenofactor_c, a*max(0.d0,
     &       (ld - ld_min)/(ld_max-ld_min)))
         if (phenofactor_c .le. 0.01d0) then
            phenostatus_c =1.d0    !no leaf
            ncd = 0.d0             !zero-out ncd once complete leaf-off occurs.
            gdd = 0.d0             !zero-out gdd once complete leaf-off occurs.
         else
            phenostatus_c = 4.d0   !leaf senescence
         end if

      end if    
       
      end subroutine pheno_update_coldwoody

      !*********************************************************************  
      subroutine pheno_update_coldherbaceous(phenostatus, 
     i            fall, fall_old, soiltemp_10d, soilt_adj,
     o            sgdd, sfdd, nmsd, ld, 
     o            phenofactor_c, phenostatus_c)
!@sum Update phenology for cold-decidous herbaceous PFTs
!@sum Called from pheno_update

      use yibs_const

      !input variables
      real*8, intent(in) :: phenostatus !phenological status (refer pheno_update for details) 
      logical,intent(in) :: fall         !.true. if the season is fall
      logical,intent(in) :: fall_old     !.true. if the season is fall
      real*8, intent(in) :: soiltemp_10d !10 day running mean of soil temperature
      real*8, intent(in) :: soilt_adj    !adjustment for soil temperature threshold (soilt_max & soilt_min)
      !in/output variables
      real*8, intent(inout) :: sgdd      !growing degree days, based on soil temperature
      real*8, intent(inout) :: sfdd      !falling degree days, based on soil temperature
      real*8, intent(inout) :: nmsd      !number of mature days, based on phenostatus  
      real*8, intent(in) :: ld           !day length in minutes
      !output variables
      real*8, intent(out) :: phenofactor_c !phenological factor for cold deciduous (refer pheno_update for details) 
      real*8, intent(inout) :: phenostatus_c !phenological status for cold deciduous (refer pheno_update for details)

Cxyue
      !*SFDD threshold for leaf off
      if (abs(phenostatus_c-3.d0) .le. 0.01d0)then
          nmsd = nmsd + 1.0d0
      endif

      if (.not. fall .and. fall_old ) then
          if (phenostatus_c.ge.3.d0) then
             sgdd = phenofactor_c*sgdd_length+sgdd_threshold
             sfdd = 0.d0
             phenostatus_c = 2.d0
          else
             sgdd = 0.d0
          endif
      endif
      if (.not. fall) sfdd = 0.d0
Cxyue

      !*Leaf-on in the spring, triggered by thermal sum, based on soil temperature
      !if sgdd is larger than its threshold 
      !in the spring (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if ((phenostatus_c.lt.3.d0).and.(sgdd.gt.sgdd_threshold)) then

         phenofactor_c  
     &      = min (1.d0,(sgdd-sgdd_threshold)/sgdd_length)

         if (phenofactor_c .lt. 1.d0) then
            phenostatus_c = 2.d0    !growing leaf
         else 
            phenostatus_c = 3.d0    !leaf in full growth
            nmsd = 0.d0             !reinitialize nmsd for falling
            sfdd = 0.d0             !reinitialize sfdd for falling
         end if
         return
      end if

      !*Leaf-off in the fall 
      !*Leaf-off triggered by soil temperature
      if ((phenostatus_c.ge.3.d0).and.(sfdd.lt.sfdd_threshold)) then

         phenofactor_c = min(phenofactor_c,max(0.d0,
     &      1.0d0+(sfdd-sfdd_threshold)/sfdd_length))

         if (phenofactor_c .eq. 0.d0) then
            phenostatus_c = 1.d0    !no leaf
            sgdd = 0.d0             !zero-out sgdd once complete leaf-off occurs.
         else  
            phenostatus_c = 4.d0    !leaf senescence
         end if 

      end if
    
      end subroutine pheno_update_coldherbaceous

      !*********************************************************************  
      subroutine pheno_update_drought(phenostatus, 
     i            mature, betad_10d, betad_min, betad_max, betad_res,
     o            phenofactor_d, phenostatus_d)
!@sum Update phenology for drought-deciduous PFTs
!@sum Called from pheno_update

      use yibs_const

      !input variables
      real*8, intent(in) :: phenostatus !phenological status (refer pheno_update for details) 
      real*8, intent(in) :: mature      !X day to mature: mature=1+X/1000  (refer pheno_update for details) 
      real*8, intent(in) :: betad_10d   !10 day running mean of betad (calculated with stressH2O)  
      real*8, intent(in) :: betad_min   !betad minimum
      real*8, intent(in) :: betad_max   !betad maximum
      real*8, intent(in) :: betad_res   !betad resistance
      !output variables
      real*8, intent(out) :: phenofactor_d !phenological factor for drought deciduous (refer pheno_update for details) 
      real*8, intent(inout) :: phenostatus_d !phenological status for drought deciduous (refer pheno_update for details)

      !*Leaf-on in the spring, triggered by water stress 
      !if betad is larger than its minimum 
      !in the spring (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X), 
      !and after long enough after the leaf-off), 
      !determine the phenofactor and corresponding phenostatus.
      if ((phenostatus_d.lt.3.d0).and. (betad_10d.gt.betad_min))then

         !determine phenofactor by scaling betad with its mim, max and resistance factor
         phenofactor_d = min(1.d0,
     &      ((betad_10d-betad_min)/(betad_max-betad_min))**betad_res)

         if (phenofactor_d .ge. 0.95d0) then
            phenostatus_d = 3.d0    !leaf in full growth
         else
            phenostatus_d = 2.d0    !growing leaf
         end if

      !*Leaf-off in the fall 
      !*Leaf-off triggered by water stress
      !if betad is falling smaller than its maximum 
      !in the fall (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X)), 
      !determine the phenofactor and corresponding phenostatus.
      elseif((phenostatus_d.ge.3.d0).and. (betad_10d.lt.betad_max))then

         !determine phenofactor by scaling betad with its mim, max and resistance factor
         phenofactor_d = max(0.d0,
     &      ((betad_10d-betad_min)/(betad_max-betad_min))**betad_res)
         if (phenofactor_d .le. EPS) then
            phenostatus_d = 1.d0    !no leaf
         else
            phenostatus_d = 4.d0    !leaf senescence
         end if 

      end if    
      
      end subroutine pheno_update_drought

      !*********************************************************************  
      subroutine pheno_update_crop(jday, pday, hday,
     o           phenofactor_c, phenostatus_c)
!@sum Update phenology for crop

      real*8, intent(in) :: jday    ! julian day
      real*8, intent(in) :: pday    ! plant day
      real*8, intent(in) :: hday    ! harvest day
      real*8, intent(out):: phenofactor_c    ! cold phenology factor
      real*8, intent(out):: phenostatus_c    ! cold phenology status
       
      real*8 pd1, pd2, hd1, hd2, p, p2
 
      pd1 = pday
      pd2 = mod(pday+dgrow, 365.d0)
      hd1 = mod(hday+365-dfall, 365.d0)
      hd2 = hday
  
      if (hd2 .gt. pd1) then
         if (jday .le. hd1) then
            p = min(1.d0, max(0.d0, (jday-pd1)/(pd2-pd1)))
         else
            p = min(1.d0, max(0.d0, 1.d0-(jday-hd1)/(hd2-hd1)))
         endif
      elseif (hd1 .gt. hd2) then
         if (jday .ge. hd2 .and. jday .le. hd1) then
            p = min(1.d0, max(0.d0, (jday-pd1)/(pd2-pd1)))
         elseif (jday .gt. hd1) then
            p = min(1.d0, max(0.d0,1.d0-(jday-hd1)/(hd2-hd1+365.d0)))
         else
            p = min(1.d0, max(0.d0,-(jday-hd2)/(hd2-hd1+365.d0)))
         endif
      elseif (pd2 .gt. hd1) then
         if (jday .ge. hd2) then
            p = min(1.d0, max(0.d0, (jday-pd1)/(pd2-pd1)))
         else
            p = min(1.d0, max(0.d0, 1.d0-(jday-hd1)/(hd2-hd1)))
         endif
      else
         if (jday .ge. pd2 .and. jday .le. pd1) then
            p = min(1.d0, max(0.d0, 1.d0-(jday-hd1)/(hd2-hd1)))
         elseif (jday .gt. pd1) then
            p = min(1.d0, max(0.d0,(jday-pd1)/(pd2-pd1+365.d0)))
         else
            p = min(1.d0, max(0.d0,1.d0+(jday-pd2)/(pd2-pd1+365.d0)))
         endif
      endif
  
      If (p .le. 0.01d0) p2 = 1    ! dormancy
      If (p .ge. 0.98d0) p2 = 3    ! mature
      If (pd1 .le. pd2 .and. jday .ge. pd1 .and. jday .le. pd2)  p2 = 2
      If (pd1 .gt. pd2 .and.(jday .ge. pd1 .or.  jday .le. pd2)) p2 = 2
      If (hd1 .le. hd2 .and. jday .ge. hd1 .and. jday .le. hd2)  p2 = 4
      If (hd1 .gt. hd2 .and.(jday .ge. hd1 .or.  jday .le. hd2)) p2 = 4
        
      phenofactor_c = p
      phenostatus_c = p2

      return
 
      end subroutine pheno_update_crop


!*************************************************************************
! Calculates Live Pool changes due to litter fall and turnover

      subroutine litter_cohort(cop, Clossacc)

      type(cohort),pointer :: cop
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator
      !--Local-----------------
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) 
      integer :: pft,i

      Closs(:,:,:) = 0.d0
      pft = cop%pft
        
      do i=1,N_CASA_LAYERS   
        if (i.eq.1) then        !only top CASA layer has leaf and wood litter -PK   
          Closs(CARBON,LEAF,i) = cop%C_leaf_lit*1.0d3
          Closs(CARBON,WOOD,i) = cop%C_wood_lit*1.0d3
        else    
          Closs(CARBON,LEAF,i) = 0.d0 
          Closs(CARBON,WOOD,i) = 0.d0
        end if
        Closs(CARBON,FROOT,i) = cop%C_root_lit*1.0d3
      enddo

      call accumulate_Clossacc(pft,Closs, Clossacc)

      end subroutine litter_cohort

!*************************************************************************

      subroutine accumulate_Clossacc(pft,Closs, Clossacc)
!@sum Accumulate litter into soil carbon pools.
      integer, intent(in) :: pft
      real*8,intent(in) :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort by depth.
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS)
!Litter accumulator.
      !---Local-----
      integer :: i

      !loop through CASA layers-->cumul litter per pool per layer -PK

      do i=1,N_CASA_LAYERS
        !* Accumulate *!
        Clossacc(CARBON,LEAF,i) = Clossacc(CARBON,LEAF,i)
     &       + Closs(CARBON,LEAF,i)
        Clossacc(CARBON,FROOT,i) = Clossacc(CARBON,FROOT,i)
     &       + Closs(CARBON,FROOT,i)
        Clossacc(CARBON,WOOD,i) = Clossacc(CARBON,WOOD,i)
     &       + Closs(CARBON,WOOD,i)

        !* NDEAD POOLS *!
        Clossacc(CARBON,SURFMET,i) = Clossacc(CARBON,SURFMET,i)
     &       + Closs(CARBON,LEAF,i) * solubfract(pft)
        Clossacc(CARBON,SOILMET,i) = Clossacc(CARBON,SOILMET,i)
     &       + Closs(CARBON,FROOT,i) * solubfract(pft)
        Clossacc(CARBON,SURFSTR,i) = Clossacc(CARBON,SURFSTR,i)
     &       + Closs(CARBON,LEAF,i) * (1-solubfract(pft))
        Clossacc(CARBON,SOILSTR,i) = Clossacc(CARBON,SOILSTR,i)
     &       + Closs(CARBON,FROOT,i) * (1-solubfract(pft))
        Clossacc(CARBON,CWD,i) = Clossacc(CARBON,CWD,i)
     &       + Closs(CARBON,WOOD,i)
      end do

      !* Return Clossacc *!
      end subroutine accumulate_Clossacc


!**********************************************************************
      subroutine litter_patch(pp, Clossacc)
!@sum litter_dead.  Update soil Tpools following litterfall.
      type(patch),pointer :: pp 
      real*8,intent(in) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !gC/m^2 Litter accumulator.
      !----Local------
      integer :: i

      !* NDEAD POOLS *!
       do i=1,N_CASA_LAYERS
        pp%Tpool(CARBON,SURFMET,i) = pp%Tpool(CARBON,SURFMET,i) 
     &     + Clossacc(CARBON,SURFMET,i)
        pp%Tpool(CARBON,SOILMET,i) = pp%Tpool(CARBON,SOILMET,i) 
     &     + Clossacc(CARBON,SOILMET,i)
        pp%Tpool(CARBON,SURFSTR,i) = pp%Tpool(CARBON,SURFSTR,i)
     &     + Clossacc(CARBON,SURFSTR,i)
        pp%Tpool(CARBON,SOILSTR,i) = pp%Tpool(CARBON,SOILSTR,i) 
     &     + Clossacc(CARBON,SOILSTR,i)
        pp%Tpool(CARBON,CWD,i) = pp%Tpool(CARBON,CWD,i) 
     &     + Clossacc(CARBON,CWD,i)
       end do   !loop through CASA layers-->total C per pool per layer -PK

       end subroutine litter_patch


!*********************************************************************
      real*8 function running_mean(dtsec,numd,var,var_mean) 
!@sum Function for expiring running mean for phenological climate statistics.
!@+   Daily average, taking into account time step of function call.
      real*8, intent(in) :: dtsec
      real*8, intent(in) :: numd !number of days for running mean
      real*8, intent(in) :: var
      real*8, intent(in) :: var_mean
      real*8 :: zweight

      zweight=exp(-2.d0/(numd*86400.d0/dtsec))
      running_mean=zweight*var_mean+(1.d0-zweight)*var  
      
      end function running_mean
!*************************************************************************

      end module phenology
