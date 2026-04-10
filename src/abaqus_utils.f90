!! Imitations of Abaqus utility routines for use by UMAT subroutines.
!! Keep the original Abaqus names — UMATs call these directly.
module indr_abaqus_utils
   use stdlib_kinds, only: dp
   use indr_linalg, only: spectral_decomposition_of_symmetric
   implicit none
   private
   public :: ROTSIG, SINV, SPRINC, SPRIND, XIT

contains

   ! Rotate a stress (LSTR=1) or strain (LSTR=0) vector by rotation matrix R.
   subroutine ROTSIG(S, R, SPRIME, LSTR, NDI, NSHR)
      implicit none
      integer,  intent(in)  :: LSTR, NDI, NSHR
      real(dp), intent(in)  :: R(3,3)
      real(dp), intent(in)  :: S(NDI+NSHR)
      real(dp), intent(out) :: SPRIME(NDI+NSHR)
      real(dp) :: a(6), b(3,3)
      integer  :: ntens
      ntens = NDI + NSHR
      a(:) = 0.0_dp
      a(1:ntens) = S(:)
      if (LSTR == 1) b = reshape([a(1),a(4),a(5),a(4),a(2),a(6),a(5),a(6),a(3)], [3,3])
      if (LSTR == 0) b = reshape([a(1),a(4)/2,a(5)/2,a(4)/2,a(2),a(6)/2,a(5)/2,a(6)/2,a(3)], [3,3])
      b = matmul(matmul(R, b), transpose(R))
      if (LSTR == 1) a = [b(1,1),b(2,2),b(3,3),b(1,2),b(1,3),b(2,3)]
      if (LSTR == 0) a = [b(1,1),b(2,2),b(3,3),2*b(1,2),2*b(1,3),2*b(2,3)]
      SPRIME = a(1:ntens)
   end subroutine ROTSIG

   ! Return two stress invariants: mean stress p (SINV1) and deviatoric norm q (SINV2).
   subroutine SINV(STRESS, SINV1, SINV2, NDI, NSHR)
      implicit none
      integer,  intent(in)  :: NDI, NSHR
      real(dp), intent(in)  :: STRESS(NDI+NSHR)
      real(dp), intent(out) :: SINV1, SINV2
      real(dp) :: devia(NDI+NSHR)
      real(dp), parameter :: sq2 = 1.4142135623730950488_dp
      if (NDI /= 3) stop 'stopped because ndi/=3 in sinv'
      SINV1 = (STRESS(1) + STRESS(2) + STRESS(3)) / 3.0_dp
      devia(1:3) = STRESS(1:3) - SINV1
      devia(4:3+NSHR) = STRESS(4:3+NSHR) * sq2
      SINV2 = sqrt(1.5_dp * dot_product(devia, devia))
   end subroutine SINV

   ! Return principal values of a stress (LSTR=1) or strain (LSTR=2) vector.
   subroutine SPRINC(S, PS, LSTR, NDI, NSHR)
      implicit none
      integer,  intent(in)  :: LSTR, NDI, NSHR
      real(dp), intent(in)  :: S(NDI+NSHR)
      real(dp), intent(out) :: PS(NDI+NSHR)
      real(dp) :: A(3,3), AN(3,3), r(6)
      if (NDI /= 3) stop 'stopped because ndi/=3 in sprinc'
      r(1:3) = S(1:3)
      if (LSTR == 1 .and. NSHR > 0) r(4:3+NSHR) = S(4:3+NSHR)
      if (LSTR == 2 .and. NSHR > 0) r(4:3+NSHR) = S(4:3+NSHR) / 2.0_dp
      A = reshape([r(1),r(4),r(5),r(4),r(2),r(6),r(5),r(6),r(3)], [3,3])
      call spectral_decomposition_of_symmetric(A, PS, AN, 3)
   end subroutine SPRINC

   ! Return principal values and directions of a stress (LSTR=1) or strain (LSTR=2) vector.
   subroutine SPRIND(S, PS, AN, LSTR, NDI, NSHR)
      implicit none
      integer,  intent(in)  :: LSTR, NDI, NSHR
      real(dp), intent(in)  :: S(NDI+NSHR)
      real(dp), intent(out) :: PS(3), AN(3,3)
      real(dp) :: A(3,3), r(6)
      if (NDI /= 3) stop 'stopped because ndi/=3 in sprind'
      r(1:3) = S(1:3)
      if (LSTR == 1 .and. NSHR > 0) r(4:3+NSHR) = S(4:3+NSHR)
      if (LSTR == 2 .and. NSHR > 0) r(4:3+NSHR) = S(4:3+NSHR) / 2.0_dp
      A = reshape([r(1),r(4),r(5),r(4),r(2),r(6),r(5),r(6),r(3)], [3,3])
      call spectral_decomposition_of_symmetric(A, PS, AN, 3)
   end subroutine SPRIND

   ! Abaqus stop utility — called by UMATs to abort cleanly.
   subroutine XIT
      stop 'stopped because umat called XIT'
   end subroutine XIT

end module indr_abaqus_utils
