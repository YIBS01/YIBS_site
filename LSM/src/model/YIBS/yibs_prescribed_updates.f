      module yibs_prescribed_updates
!@sum Routines for updating prescribed (Matthews) vegetation. 
!@+   These routines work on the ycell level or lower.

!#define DEBUG TRUE  !NYK

      use yibs_types

      implicit none
      private

      public ycell_vegupdate

      contains

!************************************************************************

      subroutine ycell_vegupdate(ecp, hemi, jday
     &     ,do_giss_phenology, do_giss_lai, do_giss_albedo!,mixed_veg
     &     ,laidata, hdata, albedodata, cropsdata, init )
!@sum Main routine for prescribed vegetation structure updates
!@+   at the ycell level (and down). 
!@+   DAILY TIME STEP ASSUMED. Coordinate-dependent
!@+   All var parameters except ycell are optional. 
!@+   Var arrays have pointer attribute to provide a way to tell the 
!@+   program that an argument is actually optional and not missing
!@+   (see how it is used in yibs_prescribe_vegupdate)
      use yibs_pfts, only: COVEROFFSET
      use patches, only : summarize_patch
      use ycells,only : summarize_ycell
      use yibs_prescr_veg, only : prescr_veg_albedo
      implicit none
      type(ycelltype),pointer :: ecp
      integer,intent(in) :: jday
      integer,intent(in) :: hemi
      logical, intent(in) :: do_giss_phenology
      logical, intent(in) :: do_giss_lai
      logical, intent(in) :: do_giss_albedo
!      logical, intent(in) :: mixed_veg
      real*8,  pointer :: laidata(:)  !Array of length N_PFT
      real*8,  pointer :: hdata(:)  !Array of length N_PFT
      real*8,  pointer :: albedodata(:,:)
      real*8,  pointer :: cropsdata
      logical, intent(in) :: init
      !----Local------
      type(patch),pointer :: pp
      
      !* 1. Update crops to get right patch/cover distribution.
      !*     NOTE:  CARBON CONSERVATION NEEDS TO BE CALCULATED FOR CHANGING VEG/CROP COVER ##
      !* 2. Update height to get any height growth (with GISS veg, height is static)
      !* 3. Update LAI, and accumulate litter from new LAI and growth/senescence.
      !*       Cohort litter is accumulated to the patch level.
      !*    3a. If external LAI, then litter is calculated based on that external LAI change.
      !*    3b. If GISS prescribed LAI, then new LAI is calculated, and then litter.
      !* 4. Update albedo based on new vegetation structure.

      !* VEGETATION STRUCTURE AND LITTER *!
      ! veg structure update with external data if provided

      ecp%jday = jday

      if (.not.do_giss_lai) then

        if ( associated(hdata) )
     &       call ycell_update_height(ecp, hdata)!, mixed_veg)
        !print *, "update hdata: ", associated(hdata) !##debug

        if ( associated(laidata) )
     &       call ycell_update_lai_poolslitter(ecp,laidata,
     &       init,jday, hemi)!,mixed_veg)
        !print *, "update laidata: ",associated(laidata) !##debug

      endif
      ! or veg structure from prescribed GISS LAI phenology 
      if ( do_giss_phenology ) then !do_giss_phenology is redundant with do_giss_lai.
        if ( hemi<-2 .or. jday <-2 )
     &       call stop_model("ycell_vegupdate: needs hemi,jday",255)
        pp => ecp%oldest
        do while (ASSOCIATED(pp))
        !* LAI, SENESCEFRAC *!
          if (do_giss_lai) 
     &         call prescr_phenology(jday,hemi, pp, do_giss_lai)
         !call summarize_patch(pp) !* Redundant because summarize_ycell is called.
          pp => pp%younger
        end do
      endif


      !* ALBEDO *!
      if ( associated(albedodata) ) then
        call ycell_update_albedo(ecp, albedodata)
        !print *, "update albedodata from array"
      else
        pp => ecp%oldest
        do while (ASSOCIATED(pp))
          ! update if have vegetation or not prognostic albedo
          if ( ASSOCIATED(pp%tallest).and.do_giss_albedo )
     &         call prescr_veg_albedo(hemi, pp%tallest%pft+COVEROFFSET, 
     &         jday, pp%albedo)
          pp => pp%younger
        end do
        !print *, "update albedodata hemi"
      endif

      call summarize_ycell(ecp)

      end subroutine ycell_vegupdate

!************************************************************************
      subroutine ycell_update_lai_poolslitter( ecp,
     i    laidata,init,jday, hemi)!,mixed_VEG)
      use yibs_pfts
      use phenology, only : litter_patch
      use yibs_prescr_veg, only : prescr_calc_lai
      implicit none
      type(ycelltype) :: ecp  !?pointer? ##
      integer,intent(in) :: jday
      integer,intent(in) :: hemi
      real*8,intent(in) :: laidata(N_PFT) !@var LAI for all PFT's 
      logical,intent(in) :: init!,mixed_VEG
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch
      type(cohort), pointer :: cop !@var current cohort
      real*8 :: laipatch, lai_old,lai_new
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.

      laipatch = 0.d0
      lai_new = 0.d0
      lai_old = 0.d0
!      Clossacc(:,:,:) = 0.d0  !?Moved to inside patch loop - NK
      pp => ecp%oldest      
      do while ( associated(pp) )
        Clossacc(:,:,:) = 0.d0
        laipatch = 0.d0
        cop => pp%tallest
        do while ( associated(cop) )
#ifndef ACTIVE_GROWTH
          cop%LAI = laidata(cop%pft)
#endif
          lai_new = cop%LAI
          !* Update biomass pools, senescefrac, accumulate litter
          call prescr_veglitterupdate_cohort(cop,lai_new,Clossacc)

          laipatch = laipatch + cop%lai
          cop => cop%shorter
        enddo
        call litter_patch(pp,Clossacc) 
        pp%LAI = laipatch
        pp => pp%younger
      enddo

      end subroutine ycell_update_lai_poolslitter


      subroutine ycell_update_height( ecp,
     i    hdata)!,mixed_VEG)
!@sum sets prescribed canopy height for ycell subgrid fractions.
      use yibs_prognostic, only: a_ws, a_wl, b_wl, eta_sl, sigl, ipft
      use yibs_pfts
      type(ycelltype) :: ecp
      real*8,intent(in) :: hdata(N_PFT) !@var LAI for all PFT's 
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch
      type(cohort), pointer :: cop !@var current cohort
      real*8 :: cpool(N_BPOOLS)
      real*8 :: h, lai_bal
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      integer :: nn

      cpool(:) = 0.d0
      pp => ecp%oldest      

      do while ( associated(pp) )
        Clossacc(:,:,:) = 0.d0
        cop => pp%tallest
        do while ( associated(cop) )
#ifndef ACTIVE_GROWTH
          cop%h = hdata(cop%pft)
#endif
          h  = cop%h
          nn = ipft(cop%pft)
          lai_bal = (a_ws(nn)*eta_sl(nn)*h/a_wl(nn))
     &            **(1.d0/(b_wl(nn)-1.d0))
          cop%C_leaf = sigl(nn)*lai_bal
          cop%C_root = leaf
          cop%C_wood = a_wl(nn)*(lai_bal**b_wl(nn))
 
          cop => cop%shorter
          print *,'After cop=>cop%shorter'
        enddo
        pp => pp%younger
      enddo
      end subroutine ycell_update_height


      subroutine ycell_update_albedo( ecp,
     i    albedodata)
!@sum sets prescribed albedo in vegetated patches of the cell (skips bare soil)
!@+   This subroutine assumes one cohort per patch !!!
      type(ycelltype) :: ecp
      real*8,intent(in) :: albedodata(N_BANDS,N_PFT) !@var albedo for all PFTs 
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch

      pp => ecp%oldest      
      do while ( associated(pp) )
        ! update albedo of vegetated patches only
        if ( associated(pp%tallest) ) then ! have vegetation
          if ( pp%tallest%pft > N_PFT .or.  pp%tallest%pft < 1 )
     &         call stop_model("ycell_update_albedo: bad pft",255)
          pp%albedo(1:N_BANDS) = albedodata(1:N_BANDS, pp%tallest%pft)
        endif
        pp => pp%younger
      enddo

      end subroutine ycell_update_albedo

!******************************************************************

      subroutine prescr_phenology(jday,hemi,pp,do_giss_lai)
!@sum Prescribed (Matthews 1983) phenology (LAI seasonality)
!@+   and associated litterfall.
!@+   DAILY TIME STEP.
!@+   Calculate new LAI, biomass poos, and senescefrac 
!@+      for given jday, for prescribed vegetation structure.
      use yibs_pfts
      use yibs_prescr_veg, only : prescr_calc_lai
     &     ,prescr_veg_albedo
      use phenology, only : litter_patch
      implicit none
      integer,intent(in) :: jday !Day of year.
      integer,intent(in) :: hemi !@var hemi -1: S.hemisphere, 1: N.hemisphere
      type(patch),pointer :: pp
      logical, intent(in) :: do_giss_lai
      !-------local-----
      type(cohort),pointer :: cop
      real*8 :: laipatch
      real*8 :: cpool(N_BPOOLS)
      real*8 :: lai_new
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.

      if (ASSOCIATED(pp)) then
        !* LAI *AND* BIOMASS - carbon pools *!
        laipatch = 0.d0         !Initialize for summing
        cpool(:) = 0.d0
        Clossacc(:,:,:) = 0.d0
        
        cop => pp%tallest
        do while (ASSOCIATED(cop))
          if ( do_giss_lai ) then  !This if statement is now redundant
            lai_new = prescr_calc_lai(cop%pft+COVEROFFSET, jday, hemi)
            call prescr_veglitterupdate_cohort(cop,
     &           lai_new,Clossacc)
            laipatch = laipatch + cop%LAI
         endif
          cop => cop%shorter
        end do
        call litter_patch(pp,Clossacc) !Update Tpools following litter accumulation. Daily time step
        pp%LAI = laipatch

      endif
      end subroutine prescr_phenology

!******************************************************************************
      subroutine prescr_veglitterupdate_cohort(
     &     cop,lai_new,Clossacc)
!@sum Given new LAI, update cohort biomass pools, litter, senescefrac.
      use phenology, only : litter_cohort
      use yibs_prognostic, only: a_ws, a_wl, b_wl, eta_sl, sigl, ipft
      implicit none
      type(cohort),pointer :: cop
      real*8,intent(in) :: lai_new
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
       !-------local-----
      real*8 :: h, lai_bal
      integer :: nn

      cop%LAI = lai_new  
      h  = cop%h
      nn = ipft(cop%pft)
      lai_bal=(a_ws(nn)*eta_sl(nn)*h/a_wl(nn))**(1.d0/(b_wl(nn)-1.d0))
      cop%C_leaf = sigl(nn)*lai_bal
      cop%C_root = leaf
      cop%C_wood = a_wl(nn)*(lai_bal**b_wl(nn))
 
      call litter_cohort(cop, Clossacc)
 
      end subroutine prescr_veglitterupdate_cohort
!******************************************************************************

      end module yibs_prescribed_updates
