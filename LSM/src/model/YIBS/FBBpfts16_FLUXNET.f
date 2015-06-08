      module FarquharBBpspar
!@sum Module YIBS plant functional type (PFT) parameters for YIBS 16 PFTs tailored to
!@+   Fluxnet sites in YIBS_1D runs.
!@+   for Farqhuar-von Caemmerer (1982) photosynthesis and 
!@+   Ball-Berry (1985) stomatal conductance.

      use yibs_const, only : N_PFT

      implicit none

      !======DECLARED TYPES====== !
      type pspartype
      integer :: pst            !Photosynth type.  1-C3, 2=C4
      real*8 :: PARabsorb       !Leaf PAR absorptance (fraction)

      !Photosynthesis/Conductance - Farquhar/Ball-Berry parameters
      real*8 :: Vcmax           !Maximum photosynthetic capacity (umol m-2 s-1)
      real*8 :: m               !Slope of Ball-Berry equation
      real*8 :: b               !Intercept of Ball-Berry equation (mol m-2 s-1)
      real*8 :: Nleaf           !g-N/m2[leaf]
      end type pspartype


      type psdrvtype
      real*8 :: ca              !Surface CO2 mole fraction (umol mol-1)
      real*8 :: ci              !Leaf internal CO2 mole fraction (umol mol-1)
      real*8 :: Tc              !Canopy (foliage) temperature (Celsius)
      real*8 :: Pa              !Atmospheric pressure (Pa)
      real*8 :: rh              !Relative humidity (fraction)
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8 :: o3              !Surface O3 mole conc. (nmol m-3)
      real*8 :: mol_to_ma       !unit convertor for atmospheric conductance
      real*8 :: mol_to_ml       !unit convertor for leaf conductance
      real*8 :: Fo3             !ozone flux to stomata (nmol O3 m-2 s-1)
      real*8 :: dFo3            !excess ozone flux to stomata (nmol O3 m-2 s-1)
#endif
      end type psdrvtype


      !=======CONSTANTS========!
!      integer,parameter :: N_PFT = 16 !In yibs_const.f

!*********************************************************************
!* YIBS PFTs
!* 1.  evergreen broadleaf early successional
!* 2.  evergreen broadleaf late successional
!* 3.  evergreen needleleaf early successional
!* 4.  evergreen needleleaf late successional
!* 5.  cold deciduous broadleaf early successional
!* 6.  cold deciduous broadleaf late successional
!* 7.  drought deciduous broadleaf
!* 8.  decidous needleleaf
!* 9.  cold adapted shrub
!* 10.  arid adapted shrub
!* 11.  C3 grass - perennial
!* 12.  C4 grass - perennial
!* 13.  C3 grass - annual
!* 14.  arctic C3 grass
!* 15.  crops - C4 herbaceous
!* 16.  crops - broadleaf woody
!********************************************************************* 
      type(pspartype),parameter :: pftpar(N_PFT) = !PFT parameters for YIBS veg types
     &!     pft PARabsorb Vcmax Kc Ko KcQ10 KoQ10 Gammastar  m b !Rdc RdH
     &     (/
     &     pspartype(1          !1. EVERGREEN BROADLEAF EARLY SUCCESSIONAL
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BET temperate & tropical, Table 3.1 (Oleson, et al 2004)
     &     ,40.d0               !Vmax25, CLM BET tropical, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2.7d0),             !Nleaf (gN/m2-leaf).

     &     pspartype(1          !2. EVERGREEN BROADLEAF LATE SUCCESSIONAL
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BET & BDT temperate & tropical, Table 3.1 (Oleson, et al 2004)
     &     ,40.d0               !Tapajo KM67 !25.d0 test
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2.7d0),             !Nleaf (gN/m2-leaf).

     &     pspartype(1          !3. EVERGREEN NEEDLELEAF EARLY SUCCESSIONAL
     &     ,.93d0               !from leaf VIS 1-albedo,CLM NET & NDT temperate & boreal, Table 3.1 (Oleson, et al 2004)
     &     ,43.d0               !Vmax25, CLM NET temperate, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2.9d0),             !Nleaf (gN/m2-leaf)

     &     pspartype(1          !4. EVERGREEN NEEDLELEAF LATE SUCCESSIONAL
     &     ,.93d0               !from leaf VIS 1-albedo,CLM NET & NDT temperate & boreal, Table 3.1 (Oleson, et al 2004)
     &     ,43.d0               !Vmax25
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,1.8d0),             !Nleaf (gN/m2-leaf).

     &     pspartype(1          !5. COLD DECIDUOUS BROADLEAF EARLY SUCCESSIONAL
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDT temperate, Table 3.1 (Oleson, et al 2004)
     &     ,45.d0               !Vmax25
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,1.07d0),             !Nleaf (gN/m2-leaf), Oak, derived from C:N 27.2 and SLA 34.5 m2/kg-C in Tatarinov & Cienciala (2006) for BIOME-BGC.

     &     pspartype(1          !6. COLD DECIDUOUS BROADLEAF LATE SUCCESSIONAL
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDT temperate, Table 3.1 (Oleson, et al 2004)
     &     ,45.d0               !Vmax25
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,1.07d0),            !Nleaf (gN/m2-leaf).  Q. ilex, Mediavilla & Escudero(2003)

     &     pspartype(1          !7. DROUGHT DECIDUOUS BROADLEAF 
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDT temperate & tropical, Table 3.1 (Oleson, et al 2004)
     &     ,45.d0               !Vcmax25, 
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,3.1d0),             !Nleaf (gN/m2-leaf), Q. ilex, Mediavilla & Escudero(2003)

     &     pspartype(1          !8. DECIDUOUS NEEDLELEAF
     &     ,.93d0               !from leaf VIS 1-albedo,CLM NDT boreal, Table 3.1 (Oleson, et al 2004)
     &     ,43.d0               !Vmax25
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,1.2d0),             !Nleaf (gN/m2-leaf). 

     &     pspartype(1          !9. COLD ADAPTED SHRUB (TUNDRA)
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDS boreal, Table 3.1 (Oleson, et al 2004)
     &     ,33.d0               !Vmax25
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,1.6d0),             !Nleaf (gN/m2-leaf).

     &     pspartype(1          !10. ARID ADAPTED SHRUB
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDS temperate, Table 3.1 (Oleson, et al 2004)
     &     ,38.d0               !Vmax25
     &     ,9.d0                !m, CLM
     &     ,.002d0              !b, CLM
     &     ,2.38d0),            !Nleaf (gN/m2-leaf).

     &     pspartype(1          !11. GRASSC3 - perennial
     &    ,.86d0                !PARabsorb, Collatz et al. (1991)
     &     ,43.d0               !Vmax25, CLM 
     &     ,9.d0               !m
     &     ,.002d0              !b
     &     ,2.46d0),            !Nleaf (gN/m2-leaf). Ponca winter wheat.

     &     pspartype(2          !12. GRASSC4 - perennial
     &     ,.9d0                !leaf VIS albedo,CLM C4 grass, Table 3.1 (Oleson, et al 2004)
     &     ,24.d0               !Vmax25
     &     ,5d0                 !m
     &     ,.002d0              !b
     &     ,1.8d0),             !Nleaf (gN/m2-leaf).

     &     pspartype(1          !13. GRASSC3 - annual
     &    ,.86d0                !PARabsorb, Collatz et al. (1991)
     &     ,43.d0               !Vmax25
     &     ,9.d0               !m
     &     ,.002d0              !b
     &     ,2.46d0),            !Nleaf (gN/m2-leaf) 

     &     pspartype(1          !14. GRASSC3 - arctic
     &     ,.89d0               !from leaf VIS 1-albedo,CLM C3 grass, Table 3.1 (Oleson, et al 2004)
     &     ,43d0                !Vmax25
     &     ,9.d0                !m
     &     ,.002d0              !b
     &     ,2.46d0),         

     &     pspartype(2          !15. CROPS - C4
     &     ,.89d0               !from leaf VIS 1-albedo,CLM Crop1 & Crop2, Table 3.1 (Oleson, et al 2004)
     &     ,40d0                !Vmax25
     &     ,5.d0                !m
     &     ,.002d0              !b
     &     ,2.5d0),             !Nleaf (gN/m2-leaf)

     &     pspartype(1          !16. CROPS - BROADLEAF WOODY
     &     ,.89d0               !from leaf VIS 1-albedo,CLM BDT, Table 3.1 (Oleson, et al 2004)
     &     ,40.d0               !Vmax25
     &     ,11.d0               !m
     &     ,.008d0              !b 
     &     ,3.1d0)              !Nleaf (gN/m2-leaf),
     &/)



!****************************************************************************
      end module FarquharBBpspar
