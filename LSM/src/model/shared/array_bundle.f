!@sum subroutines for packing/unpacking a selection of arrays
!@+   to/from a single bundle
!@auth I. Aleinov, D. Gueyffier

      module array_bundle
      implicit none
      private

      public lookup_str
      public ab_init, ab_add, ab_bundle, ab_unbundle, ab_copy, 
     &     get_bounds

      integer, parameter :: N_LOOKUP_MAX=256

      type lookup_record
        integer :: lm,km
        real*8, pointer :: src(:), dest(:)
        real*8, pointer :: src_w(:,:), dest_w(:,:)
      end type lookup_record

      type lookup_str
        private
        integer :: si_0, si_1, sj_0, sj_1  ! source bounds
        integer :: di_0, di_1, dj_0, dj_1  ! destination bounds
        integer  :: n_lookup=0
        type (lookup_record) :: lr(N_LOOKUP_MAX)
      end type lookup_str

      contains

      subroutine get_bounds(lstr, si_0, si_1, sj_0, sj_1,
     &    di_0, di_1, dj_0, dj_1)
!@sum accessor, returns domain bouds for both source and dest. grids
      implicit none
      type (lookup_str) :: lstr
      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      end subroutine get_bounds

      subroutine ab_init( lstr, si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 )
      implicit none
      type (lookup_str) :: lstr
      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      
      lstr%si_0 = si_0
      lstr%si_1 = si_1
      lstr%sj_0 = sj_0
      lstr%sj_1 = sj_1
      lstr%di_0 = di_0
      lstr%di_1 = di_1
      lstr%dj_0 = dj_0
      lstr%dj_1 = dj_1

      lstr%n_lookup = 0
      end subroutine ab_init


      subroutine ab_add( lstr, src, dest, shp, flag, src_w, dest_w )
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(*), target :: src, dest
      integer :: shp(:)
      character*(*) :: flag
      real*8, dimension(:,:), target, optional :: src_w, dest_w

      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      integer :: sim, sjm, dim, djm, lm, km

      lstr%n_lookup = lstr%n_lookup + 1
      if ( lstr%n_lookup > N_LOOKUP_MAX )
     &     call stop_model("ab_add: increase N_LOOKUP_MAX", 255)

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      sim = si_1 - si_0 + 1
      sjm = sj_1 - sj_0 + 1
      dim = di_1 - di_0 + 1
      djm = dj_1 - dj_0 + 1

      select case( flag )
      case('ij')
        lm = 1
        km = 1
      case('lij')
        lm = shp(1)
        km = 1
      case('ijk')
        lm = 1
        km = shp(3)
      case('lijk')
        lm = shp(1)
        km = shp(4)
      case default
        call stop_model("ab_add: unexpected flag", 255)
      end select

      lstr%lr(lstr%n_lookup)%lm = lm
      lstr%lr(lstr%n_lookup)%km = km

      lstr%lr(lstr%n_lookup)%src => src(1:lm*sim*sjm*km)
      lstr%lr(lstr%n_lookup)%dest => dest(1:lm*dim*djm*km)

      if ( present(src_w) .and. present(dest_w) ) then
        lstr%lr(lstr%n_lookup)%src_w => src_w(:,:)
        lstr%lr(lstr%n_lookup)%dest_w => dest_w(:,:)
      else
        if ( present(src_w) ) call stop_model(
     &       "ab_add: use both weights or none", 255)
        nullify( lstr%lr(lstr%n_lookup)%src_w )
        nullify( lstr%lr(lstr%n_lookup)%dest_w )
      endif

      end subroutine ab_add


      subroutine ab_bundle( lstr, buf_s, buf_d )
      !USE DOMAIN_DECOMP_ATM, only : agrid=>grid  !remove : for debugging only
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), pointer :: buf_s, buf_d

      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      integer im,jm,km,lm
      integer i,j,k,l,m,n,ind
      

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      im = si_1 - si_0 + 1
      jm = sj_1 - sj_0 + 1

      n = 0
      do k=1,lstr%n_lookup
        n = n+lstr%lr(k)%km*lstr%lr(k)%lm
      enddo

      allocate( buf_s(n, si_0:si_1, sj_0:sj_1) )
      allocate( buf_d(n, di_0:di_1, dj_0:dj_1) )

      n = 0
      do m = 1,lstr%n_lookup
        km = lstr%lr(m)%km
        lm = lstr%lr(m)%lm
        do k=0,km-1
          do l=0,lm-1
            n = n+1

            do j=0,jm-1
              do i=0,im-1
                ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                buf_s(n,i+si_0,j+sj_0) = lstr%lr(m)%src(ind)
              enddo
            enddo

            if ( associated(lstr%lr(m)%src_w) ) then
              do j=0,jm-1
                do i=0,im-1
                  ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                  buf_s(n,i+si_0,j+sj_0) = buf_s(n,i+si_0,j+sj_0)
     &                 * lstr%lr(m)%src_w(i+1,j+1)
                enddo
              enddo
            endif

          enddo
        enddo
      enddo

      end subroutine ab_bundle


      subroutine ab_unbundle( lstr, buf_s, buf_d )
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), pointer :: buf_s, buf_d

      integer :: di_0, di_1, dj_0, dj_1
      integer im,jm,km,lm
      integer i,j,k,l,m,n,ind

      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      im = di_1 - di_0 + 1
      jm = dj_1 - dj_0 + 1

      n = 0
      do m = 1,lstr%n_lookup
        km = lstr%lr(m)%km
        lm = lstr%lr(m)%lm
        do k=0,km-1
          do l=0,lm-1
            n = n+1

            do j=0,jm-1
              do i=0,im-1
                ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                lstr%lr(m)%dest(ind) = buf_d(n,i+di_0,j+dj_0)
              enddo
            enddo

          enddo
        enddo
      enddo

      n = 0
      do m = 1,lstr%n_lookup
        km = lstr%lr(m)%km
        lm = lstr%lr(m)%lm
        do k=0,km-1
          do l=0,lm-1
            n = n+1
            
            if ( associated(lstr%lr(m)%dest_w ) ) then
              do j=0,jm-1
                do i=0,im-1
                  ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                  if ( lstr%lr(m)%dest_w(i+1,j+1) .ne. 0.d0 ) then
                    lstr%lr(m)%dest(ind) = lstr%lr(m)%dest(ind)
     &                   / lstr%lr(m)%dest_w(i+1,j+1)
                  else
                    lstr%lr(m)%dest(ind) = 0.d0
                  endif
                enddo
              enddo
            endif

          enddo
        enddo
      enddo

      deallocate( buf_s )
      deallocate( buf_d )     

      end subroutine ab_unbundle


      subroutine ab_copy( lstr )
      implicit none
      type (lookup_str) :: lstr

      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1
      integer im,jm,km,lm,nm
      integer m

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      im = si_1 - si_0 + 1
      jm = sj_1 - sj_0 + 1

      ! make sure we copy arrays of the same size
      if ( im .ne. di_1 - di_0 + 1 .or. jm .ne. dj_1 - dj_0 + 1 )
     &     call stop_model("ab_copy: dims of arrays differ", 255)

      do m = 1,lstr%n_lookup
        km = lstr%lr(m)%km
        lm = lstr%lr(m)%lm
        nm = lm*im*jm*km
        lstr%lr(m)%dest(1:nm) = lstr%lr(m)%src(1:nm)
      enddo

      end subroutine ab_copy

      end module array_bundle

