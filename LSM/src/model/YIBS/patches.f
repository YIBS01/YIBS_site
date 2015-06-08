      module patches
!@sum Routines to organize patches in an ycell.

      use yibs_types

      implicit none


      contains
      !*********************************************************************
      subroutine insert_patch(gp, area, soil_type)
!@sum insert_patch Insert patch at youngest end of patch list.
!@+   Blank patch with no cohorts.
      implicit none
      type(ycelltype) :: gp
      real*8, intent(in) :: area
      integer, intent(in) :: soil_type
      !----
      type(patch),  pointer :: pp !Would be good to return pointer to patch.

      !* If !ASSOCIATED(gp) ERROR *!

      ! create patch
      call patch_construct(pp, gp, area, soil_type)

      !* If there exist patches, then insert.  If no patches, then allocate.
      if (ASSOCIATED(gp%youngest)) then
        !allocate(gp%youngest%younger)
        gp%youngest%younger => pp
        pp%older => gp%youngest
        gp%youngest => pp
      else !first patch on ycell
        gp%oldest => pp
        gp%youngest => pp
      end if

      end subroutine insert_patch

      !*********************************************************************
      subroutine assign_patch(pp,
     &     Ci_ini, CNC_ini, pft, Tpool_ini)
!@sum assign_patch   Initialize some canopy and soil variables for a patch.
      use yibs_const  
      
      !Eventually may want to include all patch variables as optional vars.
      type(patch),pointer :: pp
      integer, intent(in) :: pft
      real*8 :: Ci_ini, CNC_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS) ::
     &         Tpool_ini  !in g/m2 -PK
      !----Local-------
      integer :: i, n

      pp%Ci = Ci_ini
      pp%GCANOPY = CNC_ini

      !Assign soil (dead) pools.
      call assign_patch_soilcarbon(pp, pft, Tpool_ini)

      end subroutine assign_patch

      !*********************************************************************
      subroutine assign_patch_soilcarbon(pp, pft, Tpool_ini)
!@sum assign_patch_soilcarbon  Initialize.
      !Eventually may want to include all patch variables as optional vars.
      !in assign_patch and get rid of this subroutine
     
      use yibs_const  
      
      type(patch),pointer :: pp
      integer, intent(in) :: pft
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS) ::
     &         Tpool_ini  !gC/m2, soil pools only
      !----Local-------
      integer :: i, n

      !Assign soil (dead) pools.
      do n=1,N_CASA_LAYERS 
       do i=NLIVE+1,NPOOLS
        pp%Tpool(CARBON,i,n)   = Tpool_ini(pft,CARBON,i-NLIVE,n)  
        pp%Tpool_o(CARBON,i,n) = Tpool_ini(pft,CARBON,i-NLIVE,n)  
       end do
      end do

      end subroutine assign_patch_soilcarbon
      !*********************************************************************
      subroutine delete_patch(gp, pp)
!@sum Delete patch pointed to by pp
      !* NOTE:  THIS DOES NOT AUTOMATICALLY UPDATE YCELL SUMMMARY VALUES
      implicit none
      type(ycelltype) :: gp
      type(patch), pointer :: pp

      if (.NOT.ASSOCIATED(pp%older)) then !is oldest
        gp%oldest => gp%oldest%younger
        nullify(gp%oldest%older)
      else if(.NOT.ASSOCIATED(pp%younger)) then !is youngest
        gp%youngest => gp%youngest%older
        nullify(gp%youngest%younger)
      else !pp points somewhere in the middle of patch list
        pp%older%younger => pp%younger
        pp%younger%older => pp%older
      end if
      call patch_destruct(pp)

      end subroutine delete_patch

 !**************************************************************************
      real*8 function shc_patch(pp) Result(shc)
!@sum Return GISS GCM specific heat capacity for patch, but with subgrid
!@+   calculation of shc.  See comments in ycell_update_shc_mosaicveg.
!@+   Awkward legacy code.
      use yibs_const
      use yibs_pfts, only: COVEROFFSET, alamax, alamin
      !use yibs_prescr_veg, only : GISS_shc
      use allometryfn, only : GISS_shc
      implicit none
      type(patch),pointer :: pp
      !-----
      type(cohort),pointer :: cop
      real*8 :: lai
      integer :: anum

      shc = 0.d0
      lai = 0.d0
      cop => pp%shortest
      do while (associated(cop))
         anum = cop%pft + COVEROFFSET
         lai = lai + .5d0*(alamax(anum) + alamin(anum)) 
         !shc = shc + GISS_shc(cop%pft) !Alternative way.
         cop => cop%shorter
      enddo

      !shc = (.010d0+.002d0*lai+.001d0*lai**2)*shw*rhow
      shc = GISS_shc(lai)  !Note: original GISS shc is at cell level, not patch.

      end function shc_patch
!**************************************************************************

      subroutine summarize_patch(pp)
!@sum summarize_patch Calculate patch-level summary values of cohort pools.
!@+   Intensive properties (e.g. geometry, LMA) are averages weighted by
!@+   total number of individuals (may want to do by biomass of cohort)
!@+   Extensive properties (e.g. biomass, Ntot) are totals per m2 ground
      use cohorts, only: calc_CASArootfrac  !PK 7/07
      !use canopyrad
      implicit none
      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: cop
      integer :: ia  !array index
      integer :: pft
      real*8 :: fracrootCASA(N_CASA_LAYERS)  !to map fracroot to fracrootCASA

      !* Zero out cohort summary variables *!
      call zero_patch_cohortsum(pp)
      pp%Tpool(:,1:NLIVE,:) = 0.d0 !###

      if ( .not. associated(pp%tallest) ) return ! no cohorts in this patch

      do ia=1,N_COVERTYPES
        pp%LAIpft(ia) = 0.d0
        pp%HTpft(ia)  = 0.d0
        pp%GPPpft(ia) = 0.d0
        pp%IPPpft(ia) = 0.d0
        pp%MTPpft(ia) = 0.d0
        pp%Phenfpft(ia) = 0.d0
        pp%Cdead(ia)  = 0.d0
        pp%Clive(ia)  = 0.d0
      end do

      cop => pp%tallest
      do while(ASSOCIATED(cop))
        pft = cop%pft

        !assign root fractions for CASA layers -PK 11/06 
        call calc_CASArootfrac(cop%fracroot,fracrootCASA)

        pp%LAI = pp%LAI + cop%LAI  !Total
        pp%LAIpft(pft) = pp%LAIpft(pft) + cop%lai  !Total
        pp%HTpft(pft)  = pp%HTpft(pft) + cop%h   !Total
        pp%Phenfpft(pft) = pp%Phenfpft(pft) + cop%phenofactor  !Total
        do ia = NLIVE+1,NPOOLS
           pp%Cdead(pft) = pp%Cdead(pft) + pp%Tpool(CARBON,ia,1)*1.d-3
        enddo
        pp%Clive(pft) = pp%Clive(pft) + cop%C_leaf
        pp%Clive(pft) = pp%Clive(pft) + cop%C_wood
        pp%Clive(pft) = pp%Clive(pft) + cop%C_root

        pp%Ci = pp%Ci + cop%Ci * cop%LAI !wtd average
        pp%GCANOPY = pp%GCANOPY + cop%GCANOPY !total
        pp%GPP = pp%GPP + cop%GPP     !total
        pp%GPP0= pp%GPP0+ cop%GPP0    !total
        pp%GPPpft(pft) = pp%GPPpft(pft) + cop%GPP  !Total
        pp%IPPpft(pft) = pp%IPPpft(pft) + cop%IPP  !Total
        pp%MTPpft(pft) = pp%MTPpft(pft) + cop%MTP  !Total
        pp%IPP = pp%IPP + cop%IPP     !total
        pp%MTP = pp%MTP + cop%MTP     !total
        pp%NPP = pp%NPP + cop%NPP     !total
        pp%resp_r = pp%resp_r + cop%resp_r
        pp%resp_l = pp%resp_l + cop%resp_l
        pp%resp_w = pp%resp_w + cop%resp_w
        pp%resp_p = pp%resp_p + cop%resp_p
        pp%R_auto = pp%R_auto + cop%R_auto        !total
        pp%R_root = pp%R_root + cop%R_root        !total 

        pp%betad = pp%betad + cop%stressH2O*cop%LAI
        do ia=1,N_DEPTH
          pp%betadl(ia) = pp%betadl(ia) + cop%stressH2Ol(ia)*cop%LAI
        end do
        !*- - - - - end cohort summary variables - - - - - - - - - - - - -

        !* CASA pools * These are redundant but not 1-1 with the cohort pools.
       do ia=1,N_CASA_LAYERS  !loop over all CASA layers -PK 7/07
        if (ia.eq.1) then  !only top CASA layer has leaves & wood
         pp%Tpool(CARBON,LEAF,ia) = pp%Tpool(CARBON,LEAF,ia)
     &        + cop%C_leaf * 1.0d3   
         pp%Tpool(CARBON,WOOD,ia) = pp%Tpool(CARBON,WOOD,ia) 
     &        + cop%C_wood * 1.0d3
        else    
         pp%Tpool(CARBON,LEAF,ia) = 0.d0  
         pp%Tpool(CARBON,WOOD,ia) = 0.d0
        end if 
         pp%Tpool(CARBON,FROOT,ia) = pp%Tpool(CARBON,FROOT,ia)
     &        + cop%C_root * 1.0d3   
       end do  !N_CASA_LAYERS

       pp%C_total = pp%C_total + cop%C_total
       pp%C_growth = pp%C_growth + cop%C_growth

        cop => cop%shorter
      end do    !loop through cohorts

      !* ------- DO AVERAGES ----------------------------------------------*!

      !* GEOMETRY - WEIGHTED AVERAGES
c        CALL  sum_roots_cohorts2patch(pp) !froot and C_froot

      !* FLUXES 
      if (pp%LAI.gt.0.d0) then
        pp%Ci = pp%Ci / pp%LAI  !wtd average
        
        pp%betad = pp%betad/pp%LAI
        pp%betadl = pp%betadl/pp%LAI
      else
        pp%Ci = 0.d0
        pp%betad = 0.d0
        pp%betadl(:) = 0.0d0
      endif
      
      end subroutine summarize_patch

      !*********************************************************************

      subroutine zero_patch_cohortsum(pp)
!@sum Zero only the patch's cohort summary variables
      implicit none
      type(patch),pointer :: pp

      pp%LAI = 0.d0
      pp%h = 0.d0
      pp%fracroot = 0.d0
      pp%Ci = 0.d0
      pp%GCANOPY = 0.d0
      pp%GPP = 0.d0
      pp%GPP0= 0.d0
      pp%IPP = 0.d0
      pp%MTP = 0.d0
      pp%NPP = 0.d0
      pp%resp_r = 0.0d0
      pp%resp_l = 0.0d0
      pp%resp_w = 0.0d0
      pp%resp_p = 0.0d0
      pp%R_auto = 0.d0
      pp%R_root = 0.d0  !PK 5/15/07
      pp%betad = 0.d0
      pp%betadl(:) = 0.d0
      ! diags and hacks
      pp%C_total = 0.d0
      pp%C_growth = 0.d0
      end subroutine zero_patch_cohortsum
      !*********************************************************************

      subroutine zero_patch(pp)
!@sum Initialize patch, zeroing variables.
      implicit none
      type(patch),pointer :: pp
      !----Local----

            pp%age=0.d0
            pp%area = 0.d0
            pp%Reproduction(:) = 0.d0

            !* Structural variables *!
            call zero_patch_cohortsum(pp)

            !* Flux variables for GCM/EWB - patch total
            pp%z0 = 0.0d0    !## Get GISS z0 ######!
            pp%albedo = 0.0d0!## Get GISS albveg ##!
            pp%TRANS_SW = 1.0d0 !## Calculate for zero veg.##!
            pp%GCANOPY = 0.d0 !
            pp%CO2flux = 0.d0 
            pp%Ci = 0.0127D0  !Initial value not zero.
            pp%betad = 0.d0
            pp%betadl = 0.d0
            pp%Soil_resp = 0.d0

            !* Variables calculated by GCM/EWB - downscaled from grid cell
            pp%Soilmoist(:) = 0.0d0 !## Get GISS soil moisture layers ##!
!            pp%N_deposit = 0.d0

            !* Variables for biophysics and biogeochemistry

            !* DIAGNOSTIC SUMMARIES
            !* Biomass pools - patch total
            pp%LAIpft(:) = 0.d0
            pp%HTpft(:)  = 0.d0
            pp%GPPpft(:) = 0.d0
            pp%IPPpft(:) = 0.d0
            pp%MTPpft(:) = 0.d0
            pp%Phenfpft(:) = 0.d0
            pp%Cdead(:)  = 0.d0
            pp%Clive(:)  = 0.d0
            pp%Tpool(:,:,:) = 0.d0

            !* Soil type
            pp%soil_type = -1    ! set to undefined soil type (maybe use -1?)

            !* Activity diagnostics - can be summed by month, year, etc.

      end subroutine zero_patch

      !*********************************************************************
      !*********************************************************************
      subroutine sum_roots_cohorts2patch(pp)
!@sum Calculate patch-level depth- and mass-weighted average
!@sum of fine roots depth fractions, fracroot.
      implicit none
      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: cop
      integer :: n
      real*8 :: fracroot(N_DEPTH)
      real*8 :: frootC_total

      !Initialize
      do n=1,N_DEPTH  
        fracroot(n) = 0.0
      end do 
      frootC_total = 0.0

      cop => pp%tallest
      do while(ASSOCIATED(cop)) 
        !* NYK - haven't set up fracroot yet.
        do n=1,N_DEPTH
          fracroot(n) = fracroot(n) + cop%fracroot(n)
        end do
        cop => cop%shorter
      end do
      if (frootC_total > 0.0) then
        pp%fracroot = fracroot/frootC_total
      end if
!      pp%C_froot = frootC_total !Already calculated outside of this routine.

      end subroutine sum_roots_cohorts2patch
      !*********************************************************************

      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      
      subroutine reorganize_patches(ycell)
!@sum Place holder.
      implicit none
      type(ycelltype) :: ycell


      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not needed in GISS replication test.
      
      end subroutine reorganize_patches


      !*********************************************************************

      subroutine patch_construct(pp, parent_ycell, area, soil_type)
!@sum patch_construct  Allocate memory for a new patch and point to ycell
!@+   but nullify older and younger and cohort pointers.
      use cohorts, only : cohort_construct
      implicit none
      type(patch), pointer :: pp
      type(ycelltype), intent(in), target :: parent_ycell
      real*8, intent(in) :: area
      integer, intent(in) :: soil_type
      !---

      if ( area > 1.d0 ) then 
        print *, "patch area > 1, area=",area
        call stop_model("patch_construct: area >1",255)
      endif
      ! allocate memory
      allocate( pp )
      allocate( pp%LAIpft(N_COVERTYPES) )  
      allocate( pp%HTpft(N_COVERTYPES) )   
      allocate( pp%GPPpft(N_COVERTYPES) )   
      allocate( pp%IPPpft(N_COVERTYPES) )  
      allocate( pp%MTPpft(N_COVERTYPES) ) 
      allocate( pp%Phenfpft(N_COVERTYPES) )  
      allocate( pp%Cdead(N_COVERTYPES) )  
      allocate( pp%Clive(N_COVERTYPES) )  
      allocate( pp%betadl(N_DEPTH) )
      allocate( pp%fracroot(N_DEPTH) )
      allocate( pp%Reproduction(N_PFT) ) 

      ! set pointers
      pp%cellptr => parent_ycell
      nullify( pp%older )
      nullify( pp%younger )
      nullify( pp%tallest )
      nullify( pp%shortest )

      ! set variables
      call zero_patch(pp)

      pp%area = area
      pp%soil_type = soil_type
      
      end subroutine patch_construct

      !*********************************************************************

      subroutine patch_destruct(pp)
!@sum patch_destruct  Deallocate memory for a patch.
      use cohorts, only : cohort_destruct
      implicit none
      type(patch), pointer :: pp
      !---
      type(cohort),pointer :: cop, cop_tmp

      ! destroy other cohorts
      cop => pp%tallest
      do while( associated(cop) )
        cop_tmp => cop
        cop => cop%shorter
        call cohort_destruct(cop_tmp)
      enddo
      !nullify(pp%tallest)
      !nullify(pp%shortest)

      ! deallocate memory
!      deallocate( pp%Soilmoist )
      deallocate( pp )
      nullify( pp )

      end subroutine patch_destruct


      !*********************************************************************
      subroutine patch_copy(ppin,ppout)
!@sum patch_copy  Copy values in ppin into ppout, only patch-level variables.
!@+   Does not copy cohorts or ycell pointer.
!@+   Patch ppout should already be allocated.
      type(patch),intent(in) :: ppin
      type(patch),intent(out) :: ppout

!      if (associated(ppin)) then
!         if (associated(ppout)) then
         ppout%age = ppin %age
         ppout%area = ppin%area
         ppout%soil_type = ppin%soil_type
            ppout%Tpool(:,:,:) = ppin%Tpool(:,:,:)
         !* Other: nm, albedo, fluxes if desired.
!         else
!            call stop_model("patch_copy: ppout not allocated",255)
!         endif
         !If blank ppin, then do nothing.  May want to stop_model.
 !     endif
      end subroutine patch_copy


      subroutine patch_extract_pfts(pp, vfraction)
!@sum patch_extract_pfts  Extract cohort density fraction from a patch.
!@+   This is a place holder routine, as the calculation is meaningless.
      use yibs_pfts
      use yibs_const
      type(patch) :: pp
      real*8 :: vfraction(:)
      !---
      type(cohort),pointer :: cop
      real*8 :: density

      vfraction(:) = 0.d0
      density = 0.d0
      cop => pp%tallest
      do while( associated(cop) )
        vfraction( cop%pft + COVEROFFSET ) = 
     &       vfraction( cop%pft + COVEROFFSET ) + cop%n
        density = density + cop%n
        cop => cop%shorter
      enddo

      if ( density>0.d0 ) then
        vfraction(:) = vfraction(:)/density
      else
        if ( pp%soil_type == 1 ) then
          vfraction(COVER_SAND) = 1.d0
        else
          vfraction(COVER_DIRT) = 1.d0
        endif
      endif

      end subroutine patch_extract_pfts
     

!**************************************************************************

      subroutine print_Tpool(Tpool)
!@sum Print carbon pools in Tpool.
      real*8 :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS)
      print*,'Tpool(CARBON,LEAF,:)',Tpool(CARBON,LEAF,:)
      print*,'Tpool(CARBON,FROOT,:)',Tpool(CARBON,FROOT,:)
      print*,'Tpool(CARBON,WOOD,:)',Tpool(CARBON,WOOD,:)
      print*,'Tpool(CARBON,SURFMET,:)',Tpool(CARBON,SURFMET,:)
      print*,'Tpool(CARBON,SURFSTR,:)',Tpool(CARBON,SURFSTR,:)
      print*,'Tpool(CARBON,SOILMET,:)',Tpool(CARBON,SOILMET,:)
      print*,'Tpool(CARBON,SOILSTR,:)',Tpool(CARBON,SOILSTR,:)
      print*,'Tpool(CARBON,CWD,:)',Tpool(CARBON,CWD,:)
      print*,'Tpool(CARBON,SOILMIC,:)',Tpool(CARBON,SURFMIC,:)
      print*,'Tpool(CARBON,SURFMIC,:)',Tpool(CARBON,SOILMIC,:)
      print*,'Tpool(CARBON,SLOW,:)',Tpool(CARBON,SLOW,:)
      print*,'Tpool(CARBON,PASSIVE,:)',Tpool(CARBON,PASSIVE,:)
      end subroutine print_Tpool

!**************************************************************************


      subroutine patch_set_pft(pp, pft)
!@sum patch_set_pft  For mosaicked veg, insert a cohort of pft into a patch.
!!!   this subroutine works only with 0 or 1 coherts per patch !
      use cohorts, only : insert_cohort
      type(patch), pointer :: pp
      integer :: pft
      !---

      if ( .not. associated(pp%tallest) ) then
        ! no cohorts, will insert one
        call insert_cohort(pp, pft, 0.d0)
      else if ( associated(pp%tallest, pp%shortest) ) then
        ! one cohort - just change pft
        pp%tallest%pft = pft
      else
        ! more then one cohort
        call stop_model("patch_set_pft: called for >1 cohorts",255)
      endif
       
      
      end subroutine patch_set_pft

      function patch_has_pft( pp, pft )
!@sum Logical if a pft exists in a patch.
      logical patch_has_pft
      type(patch), intent(in) :: pp
      integer, intent(in) :: pft
      !---
      type(cohort),pointer :: cop
      
      patch_has_pft = .false.
      cop => pp%tallest
      do while( associated(cop) )
        if ( cop%pft == pft ) then
          patch_has_pft = .true.
          return
        endif
        cop => cop%shorter
      enddo

      end function patch_has_pft


      subroutine patch_delete_cohort( pp, cop_del )
!@sum patch_delete_cohort  Delete a cohort from a patch
      use cohorts, only : cohort_destruct
      type(patch), pointer :: pp
      type(cohort), pointer :: cop_del
      !---
      type(cohort),pointer :: cop, cop_next
      
      cop => pp%tallest
      do while( associated(cop) )

        cop_next => cop%shorter
        if ( associated( cop, cop_del ) ) then
          if ( associated( cop%taller ) ) then
            cop%taller%shorter => cop%shorter
          else
            pp%tallest => cop%shorter
          endif
          if ( associated( cop%shorter ) ) then
            cop%shorter%taller => cop%taller
          else
            pp%shortest => cop
          endif
          call cohort_destruct( cop )
        endif

        cop => cop_next
      enddo


      end subroutine patch_delete_cohort

!----------------------------------------------------------------------
      real*8 function patch_carbon(pp,
     &     pp_Cfol, pp_Cstem, pp_Croot, pp_Csoil)
!@sum patch_carbon (kg-C m-2) - Return biomass + soil carbon per land area
!@+   Also optionally return component carbon pools (kg-C m-2).      
      type(patch), pointer :: pp
      real*8,optional :: pp_Cfol
      real*8,optional :: pp_Cstem
      real*8,optional :: pp_Croot
      real*8,optional :: pp_Csoil
      !--Local-----------------
      type(cohort), pointer :: cop
      real*8 :: cohort_kgC, cohort_kgCm2

      integer :: i,k
      real*8 :: soil_gCm2, dc_1yr
      logical :: cpool_spinup      ! Ture: fast spinup for carbon pool

      cohort_kgC = 0.d0
      cohort_kgCm2 = 0.d0
      soil_gCm2 = 0.d0
      cpool_spinup = .false.
#ifdef FAST_SPINUP
      if (mod(pp%cellptr%jday,365) .eq. 0) cpool_spinup = .true.
#endif

      if (present(pp_Cfol)) pp_Cfol = 0.d0
      if (present(pp_Cstem)) pp_Cstem = 0.d0
      if (present(pp_Croot)) pp_Croot = 0.d0

      !* Biomass.
      cop => pp%tallest
      do while (associated(cop))
         cohort_kgCm2 = cohort_kgCm2 + cop%C_leaf
         cohort_kgCm2 = cohort_kgCm2 + cop%C_root
         cohort_kgCm2 = cohort_kgCm2 + cop%C_wood

         if (present(pp_Cfol)) pp_Cfol   = pp_Cfol + cop%C_leaf
         if (present(pp_Cstem)) pp_Cstem = pp_Cstem + cop%C_wood
         if (present(pp_Croot)) pp_Croot = pp_Croot + cop%C_root

         cop => cop%shorter
      end do

      if (cpool_spinup) then
         do i=1,N_CASA_LAYERS
            do k=(NLIVE+1),NPOOLS
               dc_1yr  = pp%Tpool(CARBON,k,i) - pp%Tpool_o(CARBON,k,i)
               if (dc_1yr .lt. 0.0d0) then
                  dc_1yr = Max(-pp%Tpool(CARBON,k,i)+0.1d0, dc_1yr/2.d0)
               endif
               pp%Tpool(CARBON,k,i)   = pp%Tpool(CARBON,k,i) + dc_1yr
               pp%Tpool_o(CARBON,k,i) = pp%Tpool(CARBON,k,i)
            enddo
         enddo
      endif

      !* Soil.
      do i=1,N_CASA_LAYERS
         do k=(NLIVE+1),NPOOLS
            soil_gCm2 = soil_gCm2 + pp%Tpool(CARBON,k,i)
         enddo
      enddo
      if (present(pp_Csoil)) pp_Csoil = soil_gCm2*1.d-3

      !* Total patch carbon per m2-ground.
      patch_carbon = cohort_kgCm2 + soil_gCm2*1.d-3

      end function patch_carbon
!----------------------------------------------------------------------

      end module patches
