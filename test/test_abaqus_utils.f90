program test_abaqus_utils
   use stdlib_kinds, only: dp
   use indr_abaqus_utils, only: SINV, ROTSIG, SPRINC
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-10_dp

   call test_sinv_hydrostatic()
   call test_sinv_uniaxial()
   call test_rotsig_identity_stress()
   call test_rotsig_identity_strain()
   call test_sprinc_diagonal()

   if (nfail == 0) then
      print *, 'PASS  test_abaqus_utils'
   else
      print *, 'FAIL  test_abaqus_utils — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   subroutine test_sinv_hydrostatic()
      ! Hydrostatic stress [p,p,p,0,0,0] -> sinv1=p, sinv2=0
      real(dp) :: stress(6), sinv1, sinv2
      real(dp), parameter :: p = 150.0_dp
      stress = [p, p, p, 0.0_dp, 0.0_dp, 0.0_dp]
      call SINV(stress, sinv1, sinv2, 3, 3)
      if (abs(sinv1 - p) > tol) then
         print *, 'FAIL  test_sinv_hydrostatic: sinv1 expected', p, 'got', sinv1
         nfail = nfail + 1
      end if
      if (abs(sinv2) > tol) then
         print *, 'FAIL  test_sinv_hydrostatic: sinv2 expected 0, got', sinv2
         nfail = nfail + 1
      end if
   end subroutine test_sinv_hydrostatic

   subroutine test_sinv_uniaxial()
      ! Uniaxial [s,0,0,0,0,0] -> sinv1=s/3, sinv2=|s|
      real(dp) :: stress(6), sinv1, sinv2
      real(dp), parameter :: s = 300.0_dp
      stress = [s, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call SINV(stress, sinv1, sinv2, 3, 3)
      if (abs(sinv1 - s/3.0_dp) > tol) then
         print *, 'FAIL  test_sinv_uniaxial: sinv1 expected', s/3.0_dp, 'got', sinv1
         nfail = nfail + 1
      end if
      if (abs(sinv2 - s) > tol) then
         print *, 'FAIL  test_sinv_uniaxial: sinv2 expected', s, 'got', sinv2
         nfail = nfail + 1
      end if
   end subroutine test_sinv_uniaxial

   subroutine test_rotsig_identity_stress()
      ! Identity rotation leaves a stress vector unchanged (LSTR=1)
      real(dp) :: S(6), R(3,3), Sprime(6)
      integer :: i
      S = [100.0_dp, 200.0_dp, 300.0_dp, 10.0_dp, 20.0_dp, 30.0_dp]
      R = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      call ROTSIG(S, R, Sprime, 1, 3, 3)
      do i = 1, 6
         if (abs(Sprime(i) - S(i)) > tol) then
            print *, 'FAIL  test_rotsig_identity_stress: component', i, &
                     'expected', S(i), 'got', Sprime(i)
            nfail = nfail + 1
         end if
      end do
   end subroutine test_rotsig_identity_stress

   subroutine test_rotsig_identity_strain()
      ! Identity rotation leaves a strain vector unchanged (LSTR=0)
      real(dp) :: S(6), R(3,3), Sprime(6)
      integer :: i
      S = [0.01_dp, 0.02_dp, 0.03_dp, 0.004_dp, 0.006_dp, 0.008_dp]
      R = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      call ROTSIG(S, R, Sprime, 0, 3, 3)
      do i = 1, 6
         if (abs(Sprime(i) - S(i)) > tol) then
            print *, 'FAIL  test_rotsig_identity_strain: component', i, &
                     'expected', S(i), 'got', Sprime(i)
            nfail = nfail + 1
         end if
      end do
   end subroutine test_rotsig_identity_strain

   subroutine test_sprinc_diagonal()
      ! Diagonal [100,200,300,0,0,0] -> principal values are 100, 200, 300 (any order)
      real(dp) :: S(6), PS(6)
      real(dp) :: expected(3)
      logical :: found(3)
      integer :: i, j
      S = [100.0_dp, 200.0_dp, 300.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      call SPRINC(S, PS, 1, 3, 3)
      expected = [100.0_dp, 200.0_dp, 300.0_dp]
      found = .false.
      do i = 1, 3
         do j = 1, 3
            if (.not. found(j) .and. abs(PS(i) - expected(j)) < tol) then
               found(j) = .true.
               exit
            end if
         end do
      end do
      if (.not. all(found)) then
         print *, 'FAIL  test_sprinc_diagonal: principal values wrong, got', PS(1:3)
         nfail = nfail + 1
      end if
   end subroutine test_sprinc_diagonal

end program test_abaqus_utils
