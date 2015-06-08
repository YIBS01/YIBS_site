!#define DEBUG 1

      module respauto_physio
!@sum Routines to calculate autotrophic respiration and physiological status
!@+   like water stress, and diagnostics like ci.
!@auth N.Y.Kiang

      
      use yibs_const
      use yibs_pfts
      use FarquharBBpspar !eventually need to consolidate with yibs_pfts

      implicit none

      public Rdark
      public water_stress3

      contains
!################# AUTOTROPHIC RESPIRATION ######################################

      real*8 function Rdark(Vcmax)
!@sum Rdark  Leaf dark respiration, Rd (umol m-2_leaf s-1)
!@+   From S. von Caemmerer (2000) Biochemical Models of Leaf Photosynthesis,
!@+   CSIRO book.
!@auth N.Y.Kiang

      real*8, intent(in) :: Vcmax

      Rdark = 0.011d0 * Vcmax !von Caemmerer book.
      
      end function Rdark

!---------------------------------------------------------------------!
      function water_stress3(pft, nlayers, thetarel, 
     &     fracroot, fice, betadl) Result(betad)
!@sum Rodriguez-Iturbe et al. (2001) water stress function.
!@+   Version if input is relative soil water content (saturated fraction).
!@auth N.Y.Kiang
      implicit none
      integer,intent(in) :: pft  !Plant functional type number.
      integer,intent(in) :: nlayers !Number of soil layers
!      real*8,intent(in) ::  thetas(:) !Soil vol. water (vol.water/vol.soil)
      real*8,intent(in) ::  thetarel(:) !Relative soil vol. water (vol.water/vol. saturated)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(out) :: betadl(:) !Water stress in layers
      real*8 :: betad !Stress value, 0-1, 1=no stress, weighted by layers

      !Local vars
      integer :: k
      real*8 :: s  !Normalized soil moisture, s=thetas/thetasat
      real*8 :: betak  !Stress value for layer k

      !2. Rodriguez-Iturbe, Laio, & Porporato (2001 set) water stress fn 
      betad = 0.d0
      do k = 1,nlayers
        s = thetarel(k)
        if (s.ge.pfpar(pft)%sstar) then
          betak = 1.d0  !No stress
        else if ((s.lt.pfpar(pft)%sstar).and.
     &         (s.gt.(pfpar(pft)%swilt))) then
          betak = (s-pfpar(pft)%swilt)/
     &         (pfpar(pft)%sstar-pfpar(pft)%swilt) !Just linear
        else
          betak = 0.d0
        end if
        betadl(k) = (1.d0-fice(k))*fracroot(k)*betak
        betad = betad +  (1.d0-fice(k))*fracroot(k)*betak
      end do
      if (betad < EPS2) betad=0.d0

      end function water_stress3

!***************************************************************************
      end module respauto_physio
