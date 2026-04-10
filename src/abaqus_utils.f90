!! Imitations of Abaqus utility routines for use by UMAT subroutines.
!! Original Abaqus names (ROTSIG, SINV, etc.) are preserved so that UMATs
!! can call them without modification. Modern alternatives (calc_rot_sig,
!! calc_rot_eps) are also exported for new code.
module indr_abaqus_utils
   use stdlib_kinds, only: dp
   use indr_linalg, only: spectral_decomposition_of_symmetric
   implicit none
   private
   public :: ROTSIG, SINV, SPRINC, SPRIND, XIT, &
             calc_rot_sig, calc_rot_eps

contains

   ! Apply rotation R to a symmetric 3×3 tensor T: T' = R·T·Rᵀ
   pure function rotate_tensor(T, R) result(T_rot)
      implicit none
      real(dp), intent(in) :: T(3,3), R(3,3)
      real(dp) :: T_rot(3,3)
      T_rot = matmul(matmul(R, T), transpose(R))
   end function rotate_tensor

   ! Rotate a stress Voigt vector S by rotation matrix R.
   ! Stress mapping: no factor-of-2 on shear components.
   pure function calc_rot_sig(S, R) result(S_rot)
      implicit none
      real(dp), intent(in) :: S(:), R(3,3)
      real(dp) :: S_rot(size(S))
      real(dp) :: a(6), T(3,3)
      integer  :: ntens
      ntens = size(S)
      a = 0.0_dp
      a(1:ntens) = S
      T = reshape([a(1),a(4),a(5), a(4),a(2),a(6), a(5),a(6),a(3)], [3,3])
      T = rotate_tensor(T, R)
      a = [T(1,1),T(2,2),T(3,3), T(1,2),T(1,3),T(2,3)]
      S_rot = a(1:ntens)
   end function calc_rot_sig

   ! Rotate a strain Voigt vector E by rotation matrix R.
   ! Strain mapping: engineering shear (γ=2ε) is halved into tensor form,
   ! then doubled back after rotation.
   pure function calc_rot_eps(E, R) result(E_rot)
      implicit none
      real(dp), intent(in) :: E(:), R(3,3)
      real(dp) :: E_rot(size(E))
      real(dp) :: a(6), T(3,3)
      integer  :: ntens
      ntens = size(E)
      a = 0.0_dp
      a(1:ntens) = E
      T = reshape([a(1),       a(4)/2.0_dp, a(5)/2.0_dp, &
                   a(4)/2.0_dp, a(2),       a(6)/2.0_dp, &
                   a(5)/2.0_dp, a(6)/2.0_dp, a(3)], [3,3])
      T = rotate_tensor(T, R)
      a = [T(1,1), T(2,2), T(3,3), &
           2.0_dp*T(1,2), 2.0_dp*T(1,3), 2.0_dp*T(2,3)]
      E_rot = a(1:ntens)
   end function calc_rot_eps

   ! Abaqus ROTSIG interface — thin wrapper around calc_rot_sig / calc_rot_eps.
   ! LSTR=1: stress rotation.  LSTR=0: strain rotation.
   subroutine ROTSIG(S, R, SPRIME, LSTR, NDI, NSHR)
      implicit none
      integer,  intent(in)  :: LSTR, NDI, NSHR
      real(dp), intent(in)  :: R(3,3)
      real(dp), intent(in)  :: S(NDI+NSHR)
      real(dp), intent(out) :: SPRIME(NDI+NSHR)
      if (LSTR == 1) SPRIME = calc_rot_sig(S, R)
      if (LSTR == 0) SPRIME = calc_rot_eps(S, R)
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
      A = reshape([r(1),r(4),r(5), r(4),r(2),r(6), r(5),r(6),r(3)], [3,3])
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
      A = reshape([r(1),r(4),r(5), r(4),r(2),r(6), r(5),r(6),r(3)], [3,3])
      call spectral_decomposition_of_symmetric(A, PS, AN, 3)
   end subroutine SPRIND

   ! Abaqus stop utility — called by UMATs to abort cleanly.
   subroutine XIT
      stop 'stopped because umat called XIT'
   end subroutine XIT

end module indr_abaqus_utils
