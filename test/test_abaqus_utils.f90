program test_abaqus_utils
   use stdlib_kinds, only: dp
   use indr_abaqus_utils, only: SINV, ROTSIG, SPRINC, calc_rot_sig, calc_rot_eps
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-10_dp

   call test_sinv_hydrostatic()
   call test_sinv_uniaxial()
   call test_calc_rot_sig_identity()
   call test_calc_rot_eps_identity()
   call test_calc_rot_sig_90deg()
   call test_calc_rot_eps_90deg()
   call test_rotsig_delegates_stress()
   call test_rotsig_delegates_strain()
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

   subroutine test_calc_rot_sig_identity()
      ! Identity rotation leaves a stress vector unchanged
      real(dp) :: S(6), R(3,3), S_rot(6)
      integer :: i
      S = [100.0_dp, 200.0_dp, 300.0_dp, 10.0_dp, 20.0_dp, 30.0_dp]
      R = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      S_rot = calc_rot_sig(S, R)
      do i = 1, 6
         if (abs(S_rot(i) - S(i)) > tol) then
            print *, 'FAIL  test_calc_rot_sig_identity: component', i, &
                     'expected', S(i), 'got', S_rot(i)
            nfail = nfail + 1
         end if
      end do
   end subroutine test_calc_rot_sig_identity

   subroutine test_calc_rot_eps_identity()
      ! Identity rotation leaves a strain vector unchanged
      real(dp) :: E(6), R(3,3), E_rot(6)
      integer :: i
      E = [0.01_dp, 0.02_dp, 0.03_dp, 0.004_dp, 0.006_dp, 0.008_dp]
      R = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      E_rot = calc_rot_eps(E, R)
      do i = 1, 6
         if (abs(E_rot(i) - E(i)) > tol) then
            print *, 'FAIL  test_calc_rot_eps_identity: component', i, &
                     'expected', E(i), 'got', E_rot(i)
            nfail = nfail + 1
         end if
      end do
   end subroutine test_calc_rot_eps_identity

   subroutine test_calc_rot_sig_90deg()
      ! 90-degree rotation about z-axis swaps x and y stress components.
      ! R = [[0,-1,0],[1,0,0],[0,0,1]]
      ! Stress [s11,s22,s33,s12,s13,s23] -> [s22,s11,s33,-s12,-s23,s13]
      real(dp) :: S(6), R(3,3), S_rot(6), expected(6)
      S = [100.0_dp, 200.0_dp, 300.0_dp, 40.0_dp, 50.0_dp, 60.0_dp]
      R = reshape([0.0_dp,1.0_dp,0.0_dp, -1.0_dp,0.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      S_rot    = calc_rot_sig(S, R)
      expected = [S(2), S(1), S(3), -S(4), -S(6), S(5)]
      if (any(abs(S_rot - expected) > tol)) then
         print *, 'FAIL  test_calc_rot_sig_90deg: got', S_rot, 'expected', expected
         nfail = nfail + 1
      end if
   end subroutine test_calc_rot_sig_90deg

   subroutine test_calc_rot_eps_90deg()
      ! Same 90-degree rotation applied to strain vector — same index swap pattern.
      real(dp) :: E(6), R(3,3), E_rot(6), expected(6)
      E = [0.1_dp, 0.2_dp, 0.3_dp, 0.04_dp, 0.05_dp, 0.06_dp]
      R = reshape([0.0_dp,1.0_dp,0.0_dp, -1.0_dp,0.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      E_rot    = calc_rot_eps(E, R)
      expected = [E(2), E(1), E(3), -E(4), -E(6), E(5)]
      if (any(abs(E_rot - expected) > tol)) then
         print *, 'FAIL  test_calc_rot_eps_90deg: got', E_rot, 'expected', expected
         nfail = nfail + 1
      end if
   end subroutine test_calc_rot_eps_90deg

   subroutine test_rotsig_delegates_stress()
      ! ROTSIG(LSTR=1) must agree with calc_rot_sig
      real(dp) :: S(6), R(3,3), via_rotsig(6), via_func(6)
      S = [100.0_dp, 200.0_dp, 300.0_dp, 10.0_dp, 20.0_dp, 30.0_dp]
      R = reshape([0.0_dp,1.0_dp,0.0_dp, -1.0_dp,0.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      call ROTSIG(S, R, via_rotsig, 1, 3, 3)
      via_func = calc_rot_sig(S, R)
      if (any(abs(via_rotsig - via_func) > tol)) then
         print *, 'FAIL  test_rotsig_delegates_stress: ROTSIG/calc_rot_sig disagree'
         nfail = nfail + 1
      end if
   end subroutine test_rotsig_delegates_stress

   subroutine test_rotsig_delegates_strain()
      ! ROTSIG(LSTR=0) must agree with calc_rot_eps
      real(dp) :: E(6), R(3,3), via_rotsig(6), via_func(6)
      E = [0.01_dp, 0.02_dp, 0.03_dp, 0.004_dp, 0.006_dp, 0.008_dp]
      R = reshape([0.0_dp,1.0_dp,0.0_dp, -1.0_dp,0.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      call ROTSIG(E, R, via_rotsig, 0, 3, 3)
      via_func = calc_rot_eps(E, R)
      if (any(abs(via_rotsig - via_func) > tol)) then
         print *, 'FAIL  test_rotsig_delegates_strain: ROTSIG/calc_rot_eps disagree'
         nfail = nfail + 1
      end if
   end subroutine test_rotsig_delegates_strain

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
