      module yibs
      
!@sum Contains main routines to perform single time step YIBS model simulation
!@sum on a single grid cell, ycell, whose parameters are set previously by
!@sum yibs_driver routines and the subroutine yibs_model called by a 
!@sum main program.

!@auth N.Y. Kiang

      use yibs_const
      use yibs_types

      implicit none
      private
      save

      public yibs_integrate 
      public update_veg_structure

      contains
      !*********************************************************************

      subroutine yibs_integrate(dtsec, ecp, update_day, config)
!@sum Main routine to control YIBS biophysics/biogeochemistry. 
!@+   (Patch ecological dynamics TBA)
      use cohorts
      use patches
      use biophysics, only : photosynth_cond
      use soilbgc, only : soil_bgc
      use phenology, only : clim_stats, pheno_update
      use ycells, only : summarize_ycell

      implicit none
      real*8 :: dtsec  !dt in seconds
      !type(timestruct),pointer :: tt !Time in year.fraction, Greenwich Mean Time
      type(ycelltype) :: ecp
      logical :: update_day
      type(yibs_config) :: config 

      !-----local--------
      integer :: patchnum
      type(patch),pointer :: pp


      call clim_stats(dtsec,ecp,config,update_day)

      !* Dynamic phenology
      if (update_day) then 
        call update_veg_structure(ecp, config)
      endif

      !* Biophysics
      patchnum = 0
      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 
        patchnum = patchnum + 1
        !call patch_print(771,pp," ff ")
        call photosynth_cond(dtsec, pp)

        if (config%do_soilresp) then
          call soil_bgc(dtsec, pp)
          pp%carb_tot = patch_carbon(pp, pp_Csoil=pp%carb_soil)
        endif
        pp%CO2flux = -pp%NPP + pp%Soil_resp
        
        pp%age = pp%age + dtsec

          !*********** DIAGNOSTICS FOR PLOTTING ********************!
#ifdef  YIBS_1D_DIAG         
        call summarize_patch(pp)
#endif
          !*********************************************************!

        pp => pp%younger 
      end do

      call summarize_ycell(ecp)

      end subroutine yibs_integrate


      subroutine update_veg_structure(ecp, config)
!@sum Update prognostic vegetation structure (seasonal) at the end of day
      use phenology, only : pheno_update
      use ycells, only : summarize_ycell
     &     ,ycell_carbon
      use yibs_prognostic, only: yibs_growth
      implicit none
      type(ycelltype) :: ecp
      type(yibs_config) :: config 
      !-----local--------
      type(patch),pointer :: pp
      real*8 c_before, c_after


      c_before = ycell_carbon(ecp)
      !* Loop through patches
      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 
        
c        if (config%do_phenology_activegrowth) then
          call pheno_update(pp)

! yibs dynamic growth by xyue
          call yibs_growth(pp, config)
c        endif

        pp => pp%younger 
      end do

      call summarize_ycell(ecp)

      ecp%daylength(1) = ecp%daylength(2)
      ecp%daylength(2) = 0.d0
      ecp%fall_old     = ecp%fall

      c_after = ycell_carbon(ecp)

#ifdef DEBUG
      if ( abs(c_after-c_before) > 1.d-10 ) then
        write(904,*) "dC_cell ", c_after-c_before, c_before
      endif
#endif

      end subroutine update_veg_structure


      end module yibs
