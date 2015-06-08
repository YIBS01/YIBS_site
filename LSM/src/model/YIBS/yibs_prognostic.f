#include "rundeck_opts.h"
      module yibs_prognostic
!@sum Routines for updating prognostic vegetation. 
!@+   These routines work on the ycell level or lower.
!@Author: Xu Yue

      use yibs_types

      implicit none
      save

      public yibs_growth
      public leaf_fall
      public calc_npp_resp
      public calc_daylength

      integer, parameter ::     nt_pft = 5                ! 5 general types of PFTs
!                                                          'BT', 'NT', 'C3G', 'C4G', 'shrub'
      integer, parameter ::     ns_lay = 4                ! 4 vertical soil layers
      integer, parameter ::     nd_dyn = 10               ! Interval for active growth calc. (days)
      real*8,  parameter ::   secs_360 = 31104000.d0      ! Number of seconds in 360 days         
      real*8,  parameter ::       forw = 0.0d0            ! Forward timestep weighting
      real*8,  parameter ::    dt_phen = 1.0d0/360.0d0    ! Time step for phenology (years)
      logical :: init=.true.
! PFT variables
      real*8, dimension(nt_pft), parameter :: lai_min = (/
     &   1.00d0,  1.00d0,  1.00d0,  1.00d0,  1.00d0/)     ! Minimum projected LAI

      real*8, dimension(nt_pft), parameter :: lai_max = (/
     &   9.00d0,  5.00d0,  3.00d0,  3.00d0,  3.00d0/)     ! Maximum projected LAI

      real*8, dimension(nt_pft), parameter :: laip0 = (/
     &   5.00d0,  4.00d0,  2.00d0,  2.00d0,  1.00d0/)     ! Initial LAI

      real*8, dimension(nt_pft), parameter :: tupp  = (/
     &   36.0d0,  26.0d0,  36.0d0,  45.0d0,  36.0d0/)     ! Upper temperature for 
                                                          ! photosynthesis (deg C)

      real*8, dimension(nt_pft), parameter :: tlow  = (/
     &    0.0d0, -10.0d0,   0.0d0,  13.0d0,   0.0d0/)     ! Lower temperature for 
                                                          ! photosynthesis (deg C)

      real*8, dimension(nt_pft), parameter :: a_wl = (/
     &   0.95d0,  0.85d0, 0.005d0, 0.005d0,  0.10d0/)     ! Allometric coefficient relating
!                                                           woody biomass to LAI (kg C/m2)
      real*8, dimension(nt_pft), parameter :: a_ws = (/
     &   10.0d0, 10.0d0,  1.0d0, 1.0d0, 10.0d0/)            ! Woody biomass as a multiple
!                                                           of live stem biomass
      real*8, dimension(nt_pft), parameter :: b_wl = (/
     &   1.667d0, 1.667d0, 1.667d0, 1.667d0, 1.667d0/)    ! Allometric exponent relating
!                                                           woody biomass to LAI
      real*8, dimension(nt_pft), parameter :: eta_sl = (/
     &   0.01d0,  0.01d0,  0.01d0, 0.01d0,  0.01d0/)      ! Live stemwood coefficient (kg C/m/LAI)

      real*8, dimension(nt_pft), parameter :: canht  = (/
     &   19.01d0, 20.8d0, 0.79d0, 1.26d0, 1.00d0/)       ! inital canopy height (m)

      real*8, dimension(nt_pft), parameter :: sigl   = (/
     &   0.0375d0, 0.1d0, 0.025d0, 0.05d0, 0.05d0/)       ! Specific density of leaf carbon 
!                                                           (kg C/ m2 leaf)
      real*8, dimension(nt_pft), parameter :: g_grow = (/
     &   15.0d0, 20.0d0, 20.0d0, 20.0d0, 20.0d0/)         ! Rate of leaf growth (/360days) 

      real*8, dimension(nt_pft), parameter :: gleaf0 = (/
     &   0.75d0, 0.25d0, 0.75d0, 0.75d0, 0.50d0/)         ! Min. turnover rate for leaves (/360days)

      real*8, dimension(nt_pft), parameter :: g_root = (/
     &   0.75d0, 0.25d0, 0.75d0, 0.75d0, 0.50d0/)         ! Turnover rate for root biomass (/360days)

      real*8, dimension(nt_pft), parameter :: g_wood = (/
     &   0.015d0, 0.01d0,  0.20d0,  0.20d0,  0.10d0/)      ! Turnover rate for woody biomass (/360days)

      real*8, dimension(nt_pft), parameter :: tleaf0 = (/
     &   5.0d0, -40.0d0, 5.0d0,  5.0d0, -40.0d0/)         ! Temperature below which leaves fall (C)

      real*8, dimension(nt_pft), parameter :: fsmc0  = (/
     &   0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0/)          ! Moisture avail. below which leaves fall

      real*8, dimension(nt_pft), parameter :: dgl_dt = (/
     &   9.0d0,  9.0d0,  9.0d0,  9.0d0,  9.0d0/)          ! Rate of change of leaf turnover rate 
!                                                           with temperature (C) 

      real*8, dimension(nt_pft), parameter :: dgl_dm = (/
     &   0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0/)          ! Rate of change of leaf turnover rate 
!                                                           with moisture avaialbity

      real*8, dimension(nt_pft), parameter :: rootd  = (/
     &   3.0d0,  1.0d0, 0.50d0, 0.50d0, 0.50d0/)          ! Rate of change of leaf turnover rate 

      real*8, dimension(nt_pft), parameter :: nl0  = (/
     &   0.046d0, 0.033d0, 0.073d0, 0.06d0, 0.06d0/)      ! Top leaf nitrogen concentration (kg N/kg C)

      real*8, dimension(nt_pft), parameter :: nr_nl = (/
     &   0.5d0,  0.75d0,  1.00d0,  1.00d0,  0.5d0/)     ! Ratio of root nitrogen concentration to
!                                                           leaf nitrogen concentration

      real*8, dimension(nt_pft), parameter :: ns_nl = (/
     &   0.1d0,  0.1d0,  1.00d0,  1.00d0,  0.1d0/)     ! Ratio of stem nitrogen concentration to
!                                                           leaf nitrogen concentration

      real*8, dimension(nt_pft), parameter :: r_grow = (/
     &   0.20d0,  0.20d0,  0.20d0,  0.20d0,  0.15d0/)     ! Growth respiration fraction

      real*8, dimension(nt_pft), parameter :: kpar = (/
     &   0.50d0,  0.50d0,  0.50d0,  0.50d0,  0.50d0/)     ! PAR Extinction coefficient

#ifndef PFT_MODEL_YIBS         
      integer, dimension(N_PFT), parameter :: ipft = (/   ! look-up table from YIBS to prog PFT
     &     5, 3, 5, 3, 1, 2, 1, 3
#ifdef PFT9_C4
     &    ,4
#endif
     &/)
#else
      integer, dimension(N_PFT), parameter :: ipft = (/   ! look-up table from YIBS to prog PFT
     &     1, 1, 2, 2, 1, 1, 1, 2, 5, 5, 3, 4, 3, 3, 4, 3/)
#endif
      
! Soil variables
      real*8, dimension(ns_lay), parameter :: dzsoil = (/
     &   0.1d0,  0.25d0, 0.65d0, 2.0d0/)                  ! Allometric coefficient relating

      real*8, dimension(ns_lay), parameter :: v_crit = (/
     &   0.24243d0, 0.24243d0, 0.24243d0, 0.24243d0/)     ! Volumetic soil mositure content
!                                                           above which evapotraspiration 
!                                                           is not sensitive to soil water
!                                                           (m3 H2O/ m3 soil)

      real*8, dimension(ns_lay), parameter :: v_wilt = (/
     &   0.13633d0, 0.13633d0, 0.13633d0, 0.13633d0/)     ! Volumetic soil mositure content
!                                                           below which stomata close 
!                                                           (m3 H2O/ m3 soil)

      real*8, dimension(ns_lay), parameter :: v_sat = (/
     &   0.45815d0, 0.45815d0, 0.45815d0, 0.45815d0/)     ! Volumetric soil mositure content at
!                                                           saturation (m3 H2O/ m3 soil)

      contains


      subroutine init_yibs_growth(pp)

      type(patch) :: pp
      type(cohort), pointer :: cop
      integer pft, n

      pp%l_area = .true.
      cop => pp%tallest
      do while(ASSOCIATED(cop))
         pft  = cop%pft
         n    = ipft(pft)
         cop%ht_p  = canht(n)
         cop%ht_o  = canht(n)
         cop%lai_p = laip0(n)
         cop%nstep = 0
         cop%g_leaf_ac = 0.0d0
         cop => cop%shorter
      end do 

      end subroutine init_yibs_growth


      subroutine yibs_growth(pp, config)

      type(patch) :: pp
      type(yibs_config) :: config
      !--Local-----
      type(cohort), pointer :: cop
      integer :: pft        ! YIBS pft
      integer :: n          ! prog pft
      real*8  :: ht         ! plant height (m)
      real*8  :: phen       ! phenological elongation factor [0,1] (unitless)
      real*8  :: dphen      ! changes in phenological elongation factor [-1 1] (unitless)
      real*8  :: lai_bal    ! Leaf area index in balanced growth state
      real*8  :: lai        ! Predicted leaf area index
      real*8  :: leaf       ! Leaf biomass (kg C/m2)
      real*8  :: root       ! Root biomass (kg C/m2)
      real*8  :: wood       ! Woody biomass (kg C/m2)
      real*8  :: gamma0     ! Inverse of prog timestep (360 days /nd_dyn)
      real*8  :: g_leaf_360     ! Mean leaf turnover rate for 360 days
      real*8  :: g_leaf_phen    ! Rate of leaf turnover including leaf phenology for 360 days
      real*8  :: g_leaf_dr      ! Mean leaf turnover rate for driving prognotic calculation   
      real*8  :: npp_dr         ! Mean NPP for driving prognotic calculation   
      real*8  :: npp_dr_out     ! Mean NPP for driving prognotic calculation   
      real*8  :: resp_w_dr      ! Mean wood respiration for driving prognotic calculation   
      real*8  :: resp_w_dr_out  ! Mean wood respiration for driving prognotic calculation   
      real*8  :: htpatch        ! Height for a patch
      real*8  :: laipatch       ! LAI for a patch
      real*8  :: area_tot       ! Total vegetation area of patch
      real*8  :: frac           ! Vegetation fraction
      real*8  :: frac_old       ! Old vegetation fraction
      real*8  :: ratio          ! Ratio of fractional coverage before-to-after the prog. calc.
      real*8  :: dcveg          ! Change in vegetation carbon during the timestep 
      real*8  :: c_veg          ! Total carbon content of the vegetation (kg C/m2)
      real*8  :: pc_s           ! Net carbon flux available for spreading 
      real*8  :: dh_1yr         ! Height change during 1 year
      logical :: ht_spinup      ! Ture: fast spinup for tree height
      
      ! Initialization
      htpatch  = 0.d0
      laipatch = 0.d0
      area_tot = 0.d0
      frac     = pp%area
      frac_old = frac

      if (config%do_init_activegrowth) then
         call init_yibs_growth(pp)
         config%do_init_activegrowth = .false.
      endif
      ht_spinup = .false.
#ifdef FAST_SPINUP
      if (mod(pp%cellptr%jday,365) .eq. 0) ht_spinup = .true.
#endif

      cop => pp%tallest

      do while(ASSOCIATED(cop))

         pft  = cop%pft
         n    = ipft(pft)
         ht   = cop%h

!-----------------------------------------------------------------------
! Calculate the inverse coupling timestep.
!-----------------------------------------------------------------------
         gamma0 = 360.d0/dble(nd_dyn)

!----------------------------------------------------------------------
! Diagnose the balanced-growth leaf area index and the associated leaf,
! wood, root and total vegetation carbon
!----------------------------------------------------------------------

         lai_bal = (a_ws(n)*eta_sl(n)*ht/a_wl(n))**(1.d0/(b_wl(n)-1.d0))
         leaf    = sigl(n)*lai_bal
         root    = leaf
         wood    = a_wl(n)*(lai_bal**b_wl(n))
         cop%C_leaf  = leaf + cop%dC_leaf
         cop%C_root  = root + cop%dC_root
         cop%C_wood  = wood + cop%dC_wood

!-----------------------------------------------------------------------
! Update LAI based on the leaf phenological state 
!-----------------------------------------------------------------------

         phen = cop%phenofactor
         cop%lai_p = phen*lai_bal
         dphen = phen-cop%phen_old
         g_leaf_360 = cop%g_leaf_ac/dt_phen
         if (phen .lt. cop%phen_old) then        ! Leaf fall
           g_leaf_phen = -dphen/dt_phen
         else
           g_leaf_phen = phen*g_leaf_360
         endif
         cop%g_leafp_ac = cop%g_leafp_ac + g_leaf_phen*dt_phen
         cop%g_leaf_ac = 0.0d0

!-----------------------------------------------------------------------
! Update accumulative time step and related carbon flux
!-----------------------------------------------------------------------
         cop%nstep = cop%nstep + 1

!-----------------------------------------------------------------------
! Update vegetation and terrestrial carbon storage.
!-----------------------------------------------------------------------
         If (cop%nstep .eq. nd_dyn) then 

            g_leaf_dr = cop%g_leafp_ac*gamma0
            npp_dr = cop%npp_ac*gamma0
            npp_dr_out = npp_dr
            resp_w_dr = cop%resp_w_ac*gamma0
            resp_w_dr_out =resp_w_dr
 
            call vegcarb(n, gamma0, g_leaf_dr, npp_dr, resp_w_dr,
     &                   leaf, root, wood, dcveg, pc_s,
     &                   cop%dC_leaf, cop%dC_root, cop%dC_wood,
     &                   cop%C_leaf_lit, cop%C_root_lit, cop%C_wood_lit)
            cop%dC_leaf = cop%dC_leaf/(nd_dyn*1.0d0)
            cop%dC_root = cop%dC_root/(nd_dyn*1.0d0)
            cop%dC_wood = cop%dC_wood/(nd_dyn*1.0d0)
            cop%C_leaf_lit = cop%C_leaf_lit/(gamma0*nd_dyn*1.0d0)
            cop%C_root_lit = cop%C_root_lit/(gamma0*nd_dyn*1.0d0)
            cop%C_wood_lit = cop%C_wood_lit/(gamma0*nd_dyn*1.0d0)

!-----------------------------------------------------------------------
! Diagnose the new value of Canopy Height, Leaf Area Index and Total
! Vegetation Carbon
!-----------------------------------------------------------------------
            ht = wood / (a_ws(n) * eta_sl(n))
     &         * (a_wl(n)/wood)**(1.0/b_wl(n))
            lai_bal = leaf / sigl(n)
            lai = phen * lai_bal
            c_veg = leaf + root + wood

            cop%ht_p = ht
            cop%lai_p= lai
            cop%g_leafp_ac = 0.0d0
             
            ratio    = frac_old/frac
            cop%npp_ac=ratio*(npp_dr-npp_dr_out)/gamma0
            cop%resp_w_ac = ratio*(resp_w_dr-resp_w_dr_out)/gamma0
            cop%nstep = 0

c            print*, 'Successfully run vegcarb ...'

         Endif

         if (ht_spinup) then
         dh_1yr  = cop%ht_p - cop%ht_o
         if (abs(dh_1yr) .gt. 0.1d0) cop%ht_p = cop%ht_p + dh_1yr/2.d0
         cop%ht_p = Max(cop%ht_p, 0.5d0)
         cop%ht_o = cop%ht_p
c         print*, 'Height fast spinup ...'
         endif
 
#ifdef ACTIVE_GROWTH
         cop%lai  = cop%lai_p
         cop%h    = cop%ht_p
#endif
         area_tot = area_tot + frac
         htpatch  = htpatch  + cop%ht_p
         laipatch = laipatch + cop%lai_p

         cop => cop%shorter
      end do
      pp%area_p = area_tot
      pp%ht_p   = htpatch
      pp%lai_p  = laipatch
#ifdef ACTIVE_GROWTH
      pp%lai    = laipatch
      pp%h      = htpatch
#endif
   
      return
      end subroutine yibs_growth


      subroutine leaf_fall(pft, tstar, soilm, fsmc, g_leaf)

      integer :: pft        ! IN YIBS pft
      real*8  :: tstar      ! IN Surface temperature
      real*8  :: soilm(6)   ! IN Soil moisture
      real*8  :: fsmc       ! OUT Soil moisture availability factor
      real*8  :: g_leaf     ! OUT Rate of leaf turnover for 360 days
      integer :: n          ! Prog pft
      real*8  :: fm, ft     ! Soil moisture and leaf temperature amplifiers of leaf turnover
      real*8  :: sthu(ns_lay)     ! Unfrozen soil mositure content of each
!                                   layer as a fraction of saturation
      real*8  :: f_root(ns_lay)   ! Fraction of roots in each soil layer

      n    = ipft(pft)
!-----------------------------------------------------------------------
! Calculates the fraction of the plant roots within each soil layer
!-----------------------------------------------------------------------
      call root_frac (ns_lay, dzsoil, rootd(n),f_root)

!-----------------------------------------------------------------------
! Calculates the soil moisture availability factor 
!-----------------------------------------------------------------------
      sthu(1) = soilm(1)                  ! 0.1 m ==> 0.1 m
      sthu(2) = (soilm(2)+soilm(3))/2.d0  ! 0.27, 0.57 m ==> 0.35 m
      sthu(3) = soilm(4)                  ! 1.08 m ==> 1 m
      sthu(4) = (soilm(5)+soilm(6))/2.d0  ! 2, 3.5 m ==> 3 m 
      call smc_ext(ns_lay, f_root, sthu, fsmc)

      ft = 1.0d0
      fm = 1.0d0
      if (tstar .lt. tleaf0(n)) then
         ft = 1.0d0 + dgl_dt(n)*(tleaf0(n)-tstar)
      else if (fsmc .lt. fsmc0(n)) then 
         fm = 1.0d0 + dgl_dm(n)*(fsmc0(n)-fsmc)
      endif
      g_leaf = gleaf0(n)*ft*fm
 
      return
      end subroutine leaf_fall

      
      subroutine root_frac(ns, dz, rootd, f_root)     

      integer :: ns             ! IN Number of soil layers
      real*8  :: dz(ns)         ! IN Soil layer thickness (m)
      real*8  :: rootd          ! IN Rootdepth (m)
      real*8  :: f_root(ns)     ! OUT Fraction of roots in each soil layer
      real*8  :: z1, z2         ! Depth of the top and bottom of the soil layer (m)
      real*8  :: ztot           ! Total depth of soil (m)
      real*8  :: ftot           ! Normalization factor
      integer :: nn
      
      z2 = 0.0d0
      ztot = 0.0d0
      f_root(:) = 0.0d0
      do nn = 1, ns
         z1 = z2
         z2 = z2 + dz(nn)
         ztot = ztot + dz(nn)
         f_root(nn) = exp(-z1/rootd)-exp(-z2/rootd)
      enddo
      ftot = 1.0d0 - exp(-ztot/rootd)
      do nn = 1, ns
         f_root(nn) = f_root(nn)/ftot
      enddo

      return
      end subroutine root_frac

 
      subroutine smc_ext(ns, f_root, sthu, fsmc)

      integer :: ns             ! IN Number of soil layers
      real*8  :: f_root(ns)     ! IN Fraction of roots in each soil layer
      real*8  :: sthu(ns)       ! IN Unfrozen soil mositure content of each
!                                    layer as a fraction of saturation
      real*8  :: fsmc           ! OUT Soil mositure availability
      integer :: nn
      real*8  :: fsmc_l
      
      fsmc = 0.0d0
      do nn = 1, ns
         if (abs(v_crit(nn)-v_wilt(nn)) .gt. 0.0d0) then
            fsmc_l = (sthu(nn)*v_sat(nn) - v_wilt(nn))
     &             / (v_crit(nn) - v_wilt(nn))
         else
            fsmc_l = 0.0d0
         endif
         fsmc_l = Min( Max(fsmc_l, 0.0d0), 1.0d0)
         fsmc   = fsmc + f_root(nn)*fsmc_l
      enddo

      return
      end subroutine smc_ext

      
      subroutine calc_npp_resp(pft, ht, lai, tl, fsmc, gpp, rd, 
     &                         npp, resp_w, resp_r, resp_l, resp_p)

      integer :: pft         ! IN YIBS pft
      real*8  :: ht          ! IN Canopy height (m)
      real*8  :: lai         ! IN Canopy LAI (m2/m2)
      real*8  :: tl          ! In Leaf temperature (degree C)
      real*8  :: fsmc        ! In Soil mositure availability
      real*8  :: gpp         ! IN Gross primary productivity (kg C /m2/s)
      real*8  :: rd          ! IN Canopy dark respiration (kg C /m2/s)
      real*8  :: npp         ! OUT Net primary productivity (kg C /m2/s)
      real*8  :: resp_w      ! OUT Wood respiration rate (kg C /m2/s)
      real*8  :: resp_l      ! OUT Leaf respiration rate (kg C /m2/s)
      real*8  :: resp_r      ! OUT Root respiration rate (kg C /m2/s)
      real*8  :: resp_p      ! OUT Plant respiration rate (kg C /m2/s)
      real*8  :: resp_p_m    ! Plant maintenance respiration (kgC/m2/s)
      real*8  :: resp_p_g    ! Plant growth respiration (kgC/m2/s)
      integer :: n           ! Prog. PFT
      real*8  :: lai_bal     ! Leaf area index in balanced growth state
      real*8  :: leaf        ! Leaf biomass (kg C/m2)
      real*8  :: root        ! Root biomass (kg C/m2)
      real*8  :: wood        ! Wood biomass (kg C/m2)
      real*8  :: n_leaf      ! Nitrogen contents of the leaf (kg N/m2)
      real*8  :: n_root      ! Nitrogen contents of the root (kg N/m2)
      real*8  :: n_stem      ! Nitrogen contents of the stem (kg N/m2)
      real*8  :: denom       
      
      n = ipft(pft)
      lai_bal = (a_ws(n)*eta_sl(n)*ht/a_wl(n))**(1.d0/(b_wl(n)-1.d0))
      leaf    = sigl(n)*lai_bal
      root    = leaf
      wood    = a_wl(n)*(lai_bal**b_wl(n))

!-----------------------------------------------------------------------
! Calculate the total nitrogen content of the leaf, root and stem
!-----------------------------------------------------------------------
      n_leaf = nl0(n) * sigl(n) * lai
      n_root = nr_nl(n) * nl0(n) * root
      n_stem = ns_nl(n) * nl0(n) * eta_sl(n) * ht * lai

!-----------------------------------------------------------------------
! Calculate the plant maintenance respiration rate, and 
! the wood maintenance respiration rate in kg C/m2/sec
!-----------------------------------------------------------------------
      If (n .eq. 4) rd = rd*2.0          ! C4 plant has larger respiration
      denom  = (1.d0 + EXP (0.3d0 * (tl - tupp(n))))
     &       * (1.d0 + EXP (0.3d0 * (tlow(n) - tl)))
      rd   = rd/denom

      resp_p_m = rd * (n_leaf*fsmc + n_stem + n_root) / n_leaf
      resp_l   = rd * fsmc
      resp_r   = rd * n_root / n_leaf
      resp_w   = rd * n_stem / n_leaf

!-----------------------------------------------------------------------
! Calculate the total plant respiration and the Net Primary Productivity
!-----------------------------------------------------------------------
      resp_p_g = r_grow(n) * (gpp - resp_p_m)
      If (resp_p_g .lt. 0.) resp_p_g = 0.d0
      resp_p   = resp_p_m + resp_p_g
      npp      = gpp - resp_p

      return
      end subroutine calc_npp_resp

  
      subroutine vegcarb(n, gamma0, g_leaf, npp, resp_w,
     &                   leaf, root, wood, dcveg, pc_s, 
     &                   dleaf, droot, dwood,
     &                   dleaf_lit, droot_lit, dwood_lit)

      integer :: n          ! IN prog pft
      real*8  :: gamma0     ! IN Inverse of prog timestep (360 days /nd_dyn)
      real*8  :: g_leaf     ! IN Turnover rate for leaf and fine root biomass (/ 360 days)   
      real*8  :: npp        ! INOUT Net primary productivity (kg C/m2/360days)
      real*8  :: resp_w     ! INOUT Wood maintenance respiration (kg C/m2/360days)  
      real*8  :: leaf       ! INOUT Leaf biomass (kg C/m2)
      real*8  :: root       ! INOUT Root biomass (kg C/m2)
      real*8  :: wood       ! INOUT Woody biomass (kg C/m2)
      real*8  :: dcveg      ! OUT Change in vegetation carbon during the timestep (kg C/m2/timestep)
      real*8  :: pc_s       ! OUT Net carbon flux available for spreading (kg C/m2/360days)
      real*8  :: dleaf      ! OUT Change in leaf carbon during the timestep (kg C/m2/timestep)
      real*8  :: droot      ! OUT Change in root carbon during the timestep (kg C/m2/timestep)
      real*8  :: dwood      ! OUT Change in wood carbon during the timestep (kg C/m2/timestep)
      real*8  :: dleaf_lit  ! OUT Litterfall from leaf (kg C/m2/360days)
      real*8  :: droot_lit  ! OUT Litterfall from leaf (kg C/m2/360days)
      real*8  :: dwood_lit  ! OUT Litterfall from leaf (kg C/m2/360days)
      real*8  :: lai        ! Leaf area index
      real*8  :: lit_c_l    ! Local rate of Carbon Litter production (kg C/m2/360days)
      real*8  :: pc         ! Net carbon flux available to vegetation (kg C/m2/360days)
      real*8  :: pc_g       ! Net carbon flux available for growth (kg C/m2/360days)
      real*8  :: fpar       ! PAR interception factor
      real*8  :: lamg       ! Fraction of NPP available for spreading
      real*8  :: dlamg_dlai, dlit_dlai, dfpar_dlai, dnpp_dlai
      real*8  :: dpc_dlai, dpcg_dlai, drespw_dlai
      real*8  :: dl_dw      ! Rate of change of leaf carbon with wood carbon
      real*8  :: dr_dw      ! Rate of change of root carbon with wood carbon
      real*8  :: dlai_dw    ! Rate of change of LAI with wood carbon
      real*8  :: numer, denom, dlai
      real*8  :: wood_min, wood_max

      
      lai    = leaf/sigl(n)
!----------------------------------------------------------------------
! Calculate the local production rate for carbon litter
!----------------------------------------------------------------------
      lit_c_l = g_leaf*leaf + g_root(n)*root + g_wood(n)*wood
      dleaf_lit = g_leaf*leaf/lit_c_l
      droot_lit = g_root(n)*root/lit_c_l
      dwood_lit = g_wood(n)*wood/lit_c_l
!----------------------------------------------------------------------
! Diagnose the net local carbon flux into the vegetation
!----------------------------------------------------------------------
      pc = npp - lit_c_l
     
!----------------------------------------------------------------------
! Variables required for the implicit and equilibrium calculations
!----------------------------------------------------------------------
      dlit_dlai = (g_leaf*leaf + g_root(n)*root)/lai       
     &          + b_wl(n)*g_wood(n)*wood/lai
      
      fpar = (1.d0 - exp(-kpar(n)*lai))/kpar(n)
      dfpar_dlai = exp(-kpar(n)*lai)
      dnpp_dlai  = npp*dfpar_dlai/fpar 
     &           + (1-r_grow(n))*resp_w*(dfpar_dlai/fpar-b_wl(n)/lai)
      lamg = 1.d0 - (lai - lai_min(n))/(lai_max(n) - lai_min(n))
      dlamg_dlai = -1.0d0/(lai_max(n) - lai_min(n))
      pc_g = lamg*npp - lit_c_l
      dpcg_dlai = lamg*dnpp_dlai + dlamg_dlai*npp - dlit_dlai
      dpc_dlai  = dnpp_dlai - dlit_dlai

!----------------------------------------------------------------------
! Update vegetation carbon contents
!----------------------------------------------------------------------
      dcveg = leaf + root + wood
      
!----------------------------------------------------------------------
! Calculate the increment to the wood carbon
!----------------------------------------------------------------------
      dl_dw = leaf/(b_wl(n)*wood)
      dr_dw = dl_dw
      dlai_dw = dl_dw/sigl(n)

      numer = pc_g
      denom = (1+dl_dw+dr_dw)*gamma0-forw*dlai_dw*dpcg_dlai
      denom = max(denom,1.0d-6)
      dwood = numer/denom

c      print*, 'dwood:', dwood, dwood/wood*100.
c      print*, 'pc_g:' , lamg, npp, lit_c_l

!----------------------------------------------------------------------
! Ensure that the local leaf area index does not drop below its
! minimum value or exceed its maximum value.
!----------------------------------------------------------------------
      wood_min = a_wl(n)*lai_min(n)**b_wl(n)
      wood_max = a_wl(n)*lai_max(n)**b_wl(n)
      dwood = max((wood_min-wood),dwood)
      dwood = min((wood_max-wood),dwood)

!----------------------------------------------------------------------
! Diagnose the increments to leaf and root carbon
!----------------------------------------------------------------------
      dleaf = sigl(n)*((wood+dwood)/a_wl(n))**(1.0/b_wl(n)) - leaf
      droot = dleaf

!----------------------------------------------------------------------
! Update carbon contents
!----------------------------------------------------------------------
      leaf = leaf + dleaf
      root = root + droot
      wood = wood + dwood

      dcveg = leaf + root + wood - dcveg

!----------------------------------------------------------------------
! Diagnose the carbon available for spreading and apply implicit
! corrections to the driving fluxes.
!----------------------------------------------------------------------
      dlai = leaf/sigl(n) - lai
      pc_s = pc + forw*dpc_dlai*dlai - dcveg*gamma0

      fpar = (1.0d0 - exp(-kpar(n)*lai))/kpar(n)
      dfpar_dlai = exp(-kpar(n)*lai)
      drespw_dlai = resp_w*b_wl(n)/lai

      npp = npp + forw*dnpp_dlai*dlai
      resp_w = resp_w + forw*drespw_dlai*dlai
#ifndef ACTIVE_GROWTH
      dcveg = 0.0d0
#endif
      dleaf_lit = dleaf_lit*(npp - dcveg)
      droot_lit = droot_lit*(npp - dcveg)
      dwood_lit = dwood_lit*(npp - dcveg)

      return
      end subroutine vegcarb  


      real*8 function calc_daylength(jday, lat)
      integer, intent(in) :: jday             ! day of the year
      real*8, intent(in)  :: lat              ! latitude of cell
      real*8 :: dec, pie, plat

      pie = atan(1.0d0)*4.0
      plat= lat*pie/180.d0
      dec = asin(sin(23.45/180.*pie)*sin(360.*(jday-81.)/365./180.*pie))
      calc_daylength = 120.d0/15.d0*180.d0/pie   ! number of sun minutes
     &       *acos(-sin(plat)*sin(dec)/cos(plat)/cos(dec))
      
      end function calc_daylength


      end module yibs_prognostic
