program test_usolver
   ! Tests for USOLVER: solver for mixed stress/strain control
   !
   ! KK * u = rhs
   !   is(i) = 0  ->  u(i) is prescribed (strain-controlled)
   !   is(i) = 1  ->  rhs(i) is prescribed (stress-controlled)
   use stdlib_kinds, only: dp
   use mod_inc_driver_funcs, only: USOLVER
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-10_dp

   call test_pure_strain()
   call test_pure_stress()
   call test_mixed_one_stress()

   if (nfail == 0) then
      print *, 'PASS  test_usolver'
   else
      print *, 'FAIL  test_usolver — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   subroutine test_pure_strain()
      ! All strain-controlled: USOLVER computes rhs = KK * u
      ! KK = [[3,1],[1,3]],  u = [1,2]  ->  rhs = [5,7]
      real(dp) :: KK(2,2), u(2), rhs(2)
      integer  :: is(2)
      KK  = reshape([3.0_dp,1.0_dp,1.0_dp,3.0_dp], [2,2])
      u   = [1.0_dp, 2.0_dp]
      rhs = [0.0_dp, 0.0_dp]
      is  = [0, 0]
      call USOLVER(KK, u, rhs, is, 2)
      if (abs(rhs(1) - 5.0_dp) > tol .or. abs(rhs(2) - 7.0_dp) > tol) then
         print *, 'FAIL  test_pure_strain: rhs expected [5,7] got', rhs
         nfail = nfail + 1
      end if
   end subroutine test_pure_strain

   subroutine test_pure_stress()
      ! All stress-controlled: USOLVER computes u = KK^{-1} * rhs
      ! KK = [[3,1],[1,3]], det=8, KK^{-1} = (1/8)*[[3,-1],[-1,3]]
      ! rhs = [5,7]  ->  u = (1/8)*[15-7, -5+21] = [1,2]
      real(dp) :: KK(2,2), u(2), rhs(2)
      integer  :: is(2)
      KK  = reshape([3.0_dp,1.0_dp,1.0_dp,3.0_dp], [2,2])
      u   = [0.0_dp, 0.0_dp]
      rhs = [5.0_dp, 7.0_dp]
      is  = [1, 1]
      call USOLVER(KK, u, rhs, is, 2)
      if (abs(u(1) - 1.0_dp) > tol .or. abs(u(2) - 2.0_dp) > tol) then
         print *, 'FAIL  test_pure_stress: u expected [1,2] got', u
         nfail = nfail + 1
      end if
   end subroutine test_pure_stress

   subroutine test_mixed_one_stress()
      ! Mixed: component 1 stress-controlled, component 2 strain-controlled
      ! KK = [[3,1],[1,3]], u(2)=2 prescribed, rhs(1)=5 prescribed
      ! rhs1(1) = 5 - 2*KK(1,2) = 5 - 2 = 3
      ! u(1) = rhs1(1) / KK(1,1) = 3/3 = 1
      ! rhs(2) = KK(2,1)*u(1) + KK(2,2)*u(2) = 1 + 6 = 7
      real(dp) :: KK(2,2), u(2), rhs(2)
      integer  :: is(2)
      KK   = reshape([3.0_dp,1.0_dp,1.0_dp,3.0_dp], [2,2])
      u    = [0.0_dp, 2.0_dp]
      rhs  = [5.0_dp, 0.0_dp]
      is   = [1, 0]
      call USOLVER(KK, u, rhs, is, 2)
      if (abs(u(1) - 1.0_dp) > tol) then
         print *, 'FAIL  test_mixed_one_stress: u(1) expected 1 got', u(1)
         nfail = nfail + 1
      end if
      if (abs(rhs(2) - 7.0_dp) > tol) then
         print *, 'FAIL  test_mixed_one_stress: rhs(2) expected 7 got', rhs(2)
         nfail = nfail + 1
      end if
   end subroutine test_mixed_one_stress

end program test_usolver
