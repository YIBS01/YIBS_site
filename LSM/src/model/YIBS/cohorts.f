      module cohorts
!@sum Routines to organize cohorts within an ycell.

      use yibs_const
      use yibs_types

      implicit none


      contains
      !*********************************************************************
      subroutine insert_cohort(pp,pft,h, LAI, fracroot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     phenofactor_c, phenofactor_d, phenofactor, phenostatus, 
     &     betad_10d, CB_d,
     &     turnover_amp, llspan) !KIM -7 vars for phenology
!@sum insert_cohort Insert new cohort into a canopy patch.

      type(patch),pointer :: pp
      integer :: pft
      real*8, optional, intent(in) :: h, LAI
      real*8, optional, intent(in) :: fracroot(N_DEPTH),
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     phenofactor_c, phenofactor_d, phenofactor, phenostatus, 
     &     betad_10d, CB_d,
     &     turnover_amp, llspan
      !------------------
      type(cohort),pointer :: cop, csp, newc

      call cohort_construct(newc, pp, pft)

      if ( present(h) ) then
        call assign_cohort(newc,pft,h, LAI, fracroot,
     &       Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &       phenofactor_c, phenofactor_d, phenofactor,phenostatus,
     &       betad_10d, CB_d,
     &       turnover_amp, llspan)
      endif

        if (ASSOCIATED(pp%shortest)) then !A. There are other cohorts.
          if (pp%shortest%h.ge.newc%h) then !newc is shortest
            pp%shortest%shorter => newc  !changed = to => -PK 9/28/07 
            newc%taller => pp%shortest
            pp%shortest => newc
          else if (pp%tallest%h.lt.newc%h) then !newc is tallest
            pp%tallest%taller => newc  
            newc%shorter => pp%tallest
            pp%tallest => newc
          else !newc is neither tallest nor shortest
            cop => pp%shortest
            do while (cop%h.lt.newc%h) !find next taller cohort
              cop => cop%taller
            end do
            newc%taller => cop
            newc%shorter => cop%shorter
            newc%shorter%taller => newc
            cop%shorter => newc
          end if
          !Now parse through csp's
          csp => newc%taller
          if (ASSOCIATED(csp)) then
            do while (csp%pft.ne.newc%pft)
              if (ASSOCIATED(csp%taller).and.(csp%pft.ne.newc%pft)) then
                csp => csp%taller
              else
                exit !exit loop
              end if
            end do
            if (csp%pft.eq.newc%pft) then !1.newc is not tallest csp
              newc%csptaller => csp
              newc%cspshorter => csp%cspshorter
              csp%cspshorter => newc
            else !2. no taller con-specifics
              nullify(newc%csptaller)
            end if
          else  !3. no taller con-specifics
            nullify(newc%csptaller)
          end if
          if (.NOT.ASSOCIATED(newc%cspshorter)) then !Case 1 did not hold
            csp => newc%shorter
            if (ASSOCIATED(csp)) then
              do while (csp%pft.ne.newc%pft)
                if (ASSOCIATED(csp%shorter).and.
     &               (csp%pft.ne.newc%pft)) then
                  csp => csp%shorter
                else
                  exit
                end if
              end do
              if (csp%pft.eq.newc%pft) then !4. newc is not shortest csp
                newc%cspshorter => csp
                newc%csptaller => csp%csptaller
                csp%csptaller => newc
              else !5. no shorter con-specifics
                nullify(newc%cspshorter)
              end if
            else !6. no shorter con-specifics
              nullify(newc%cspshorter)
            end if
          end if
        else !B. newc is the only cohort
          pp%tallest => newc
          pp%shortest => newc
        end if

      end subroutine insert_cohort
      !*********************************************************************
      
      subroutine assign_cohort(cop,pft,h, LAI, fracroot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     phenofactor_c, phenofactor_d, phenofactor, phenostatus, 
     &     betad_10d, CB_d,
     &     turnover_amp, llspan)
!     &     stressH2O, stressH2Ol)
!@sum assign_cohort  Assign values to a cohort.

      use yibs_pfts
      use yibs_prognostic, only: a_ws, a_wl, b_wl, eta_sl, sigl, ipft
      type(cohort) :: cop
      integer :: pft, nn
      real*8,optional :: h, LAI, fracroot(:),
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     phenofactor_c, phenofactor_d, phenofactor, phenostatus,
     &     betad_10d, CB_d,
     &     turnover_amp, llspan
      real*8 :: lai_bal
!     &     stressH2O, stressH2Ol(:)
      cop%pft = pft
      cop%LAI = LAI
      cop%h = h
Cxyue
      cop%ht_p  = h
      cop%ht_o  = h
      cop%lai_p = lai
      cop%nstep = 0.d0
      cop%g_leaf_ac = 0.0d0
      nn = ipft(pft)
      lai_bal=(a_ws(nn)*eta_sl(nn)*h/a_wl(nn))**(1.d0/(b_wl(nn)-1.d0))
      cop%C_leaf = sigl(nn)*lai_bal
      cop%C_root = leaf
      cop%C_wood = a_wl(nn)*(lai_bal**b_wl(nn))
      cop%dC_leaf = 0.d0
      cop%dC_root = 0.d0
      cop%dC_wood = 0.d0
      cop%C_leaf_lit = 0.d0
      cop%C_root_lit = 0.d0
      cop%C_wood_lit = 0.d0
Cxyue
      cop%fracroot(:) = fracroot(:)
      cop%Ci = Ci
      cop%GCANOPY =  GCANOPY 
      cop%GPP =  GPP 
      cop%NPP =  NPP 
      cop%R_auto = 0.0
      cop%R_root = 0.0 !PK 5/15/07
      cop%phenofactor_c = phenofactor_c
      cop%phenofactor_d = phenofactor_d
      cop%phenofactor = phenofactor
      cop%phenostatus = phenostatus
      cop%phenostatus_c = phenostatus
      cop%phenostatus_d = phenostatus
      cop%betad_10d = betad_10d
      cop%CB_d = CB_d
      cop%turnover_amp = turnover_amp
      cop%llspan = llspan
      if (cop%llspan.eq.-999.d0 .and.
     &   (pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
     &   (pfpar(cop%pft)%leaftype.eq.BROADLEAF))
     &   cop%llspan=pfpar(cop%pft)%lrage*12.d0
      cop%C_total = 0.d0
      cop%C_growth = 0.d0
      cop%C_growth_flux = 0.d0

      end subroutine assign_cohort
      !*********************************************************************

      subroutine zero_cohort(cop)
!@sum Zero all real variables in cohort record.      
      use yibs_pfts
      type(cohort),pointer :: cop

      cop%h = 0.0
      cop%LAI = 0.0
      cop%fracroot(:) = 0.0

      cop%Ci =  0.0127D0        !Initial value not zero.
      cop%GCANOPY = 0.0
      cop%GPP = 0.0
      cop%IPP = 0.0
      cop%MTP = 0.0
      cop%NPP = 0.0
      cop%R_auto = 0.0
      cop%R_root = 0.0 !PK 5/15/07

      !* PHENOLOGY/GROWTH *!
      !KIM - starting in the middle of winter for cold-dec.
      cop%phenofactor_c = 0.d0
      cop%phenofactor_d = 1.d0
      cop%phenofactor = 0.d0
      cop%phenostatus = 1.d0
      cop%phenostatus_c = 1.d0
      cop%phenostatus_d = 1.d0
      cop%betad_10d = 1.d0
      cop%CB_d = 0.d0
      cop%turnover_amp = 1.d0
      cop%llspan = -999.d0
      cop%Sacclim = 25.d0 !NK - force mild average temperatures default.
      cop%stressH2O = 1.d0 !Default no stress.
      cop%stressH2Ol(:) = 1.d0 !Default no stress.
      cop%C_total = 0.d0
      cop%C_growth = 0.d0
      cop%C_growth_flux = 0.d0

      end subroutine zero_cohort


      !*********************************************************************
       
      subroutine reorganize_cohorts(pp)
!@sum Place holder.
      type(patch),pointer :: pp

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not need in GISS replication test.
      
      end subroutine reorganize_cohorts
      !*********************************************************************


      subroutine cohort_construct(cop, parent_patch, pnum)
!@sum cohort_construct  Create a cohort with default values. if optional values
!@+ are provided - set them
!@auth I.Aleinov.
      ! This function may eventually be combined with assign_cohort
      ! for better performance
      type(cohort),pointer :: cop
      integer, optional :: pnum
      type(patch), optional, target :: parent_patch

      ! allocate memory
      allocate( cop )
      allocate( cop%fracroot(N_DEPTH) )
      allocate( cop%stressH2Ol(N_DEPTH) )

      ! just in case...
      nullify( cop%height_dz )
      nullify( cop%fp_dz )
      nullify( cop%height )
      nullify( cop%fp )

      ! set pointers if any
      nullify(cop%pptr )
      nullify(cop%taller )
      nullify(cop%shorter )
      nullify(cop%csptaller )
      nullify(cop%cspshorter )
      if ( present(parent_patch) ) then
        cop%pptr => parent_patch
      endif

      ! set variables
      cop%pft = -1              ! = -1 if pft not set
      if ( present(pnum) ) cop%pft = pnum
      cop%n = 0.0

      call zero_cohort(cop)

      end subroutine cohort_construct
      !*********************************************************************


      subroutine cohort_destruct(cop)
!@sum cohort_destruct Deallocate memory used by cohort
!@auth I.Aleinov.
      type(cohort),pointer :: cop

      ! we may want ot collapse hole between "taller" and "shorter"
      ! here if this functionality is needed

      ! deallocate all memory
      deallocate( cop%stressH2Ol )
      deallocate( cop%fracroot )
      deallocate( cop )
      nullify( cop )

      end subroutine cohort_destruct


      !*********************************************************************
      subroutine calc_CASArootfrac(copfracroot,fracrootCASA)  !PK 11/06
!@sum calc_CASArootfrac  Maps fracroot(N_DEPTH) to fracrootCASA(N_CASA_LAYERS)
!@+   ifdef customization required dependent on thicknesses of CASA layers
!@+   and GCM layers 
      !type(cohort),intent(in) :: cop
      real*8,pointer :: copfracroot(:)
      real*8,intent(out) :: fracrootCASA(N_CASA_LAYERS)

      if (N_CASA_LAYERS == 1) then
         fracrootCASA = 1.d0  !if there is no explicit depth structure
      else
#ifdef NCASA2
      !***scheme for N_CASA_LAYERS=2 (layers: 0-30, 30-100 cm)*** 
         fracrootCASA(1) = copfracroot(1) + copfracroot(2)  !CASA layer 1 --> GISS GCM layers 1,2
         fracrootCASA(2) = copfracroot(3) + copfracroot(4)  !CASA layer 2 --> GISS layers 3,4
     &                + copfracroot(5)                    !need to add 5th GISS layer (mainly for trees) -PK 6/26/07
#endif
      end if
                               
      end subroutine calc_CASArootfrac
      !*********************************************************************

      end module cohorts
