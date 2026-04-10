program test_maps
   use stdlib_kinds, only: dp
   use indr_maps, only: map2T, map2stress, map2D, map2stran
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-12_dp

   call test_stress_roundtrip()
   call test_strain_roundtrip()
   call test_map2stress_known()
   call test_map2stran_known()

   if (nfail == 0) then
      print *, 'PASS  test_maps'
   else
      print *, 'FAIL  test_maps — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   subroutine test_stress_roundtrip()
      ! map2stress(map2T(v)) should recover v (stress/tensor pair, no factor-of-2)
      real(dp) :: v(6), T(3,3), v2(6)
      integer :: i
      v = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      T  = map2T(v, 6)
      v2 = map2stress(T, 6)
      do i = 1, 6
         if (abs(v2(i) - v(i)) > tol) then
            print *, 'FAIL  test_stress_roundtrip: component', i, &
                     'expected', v(i), 'got', v2(i)
            nfail = nfail + 1
         end if
      end do
   end subroutine test_stress_roundtrip

   subroutine test_strain_roundtrip()
      ! map2stran(map2D(v)) should recover v (strain/tensor pair, factor-of-2 on shear)
      real(dp) :: v(6), D(3,3), v2(6)
      integer :: i
      v = [1.0_dp, 2.0_dp, 3.0_dp, 8.0_dp, 10.0_dp, 12.0_dp]
      D  = map2D(v, 6)
      v2 = map2stran(D, 6)
      do i = 1, 6
         if (abs(v2(i) - v(i)) > tol) then
            print *, 'FAIL  test_strain_roundtrip: component', i, &
                     'expected', v(i), 'got', v2(i)
            nfail = nfail + 1
         end if
      end do
   end subroutine test_strain_roundtrip

   subroutine test_map2stress_known()
      ! map2stress on a diagonal tensor gives the diagonal elements
      real(dp) :: T(3,3), v(6)
      T = reshape([100.0_dp,0.0_dp,0.0_dp, 0.0_dp,200.0_dp,0.0_dp, 0.0_dp,0.0_dp,300.0_dp], [3,3])
      v = map2stress(T, 6)
      if (abs(v(1) - 100.0_dp) > tol .or. abs(v(2) - 200.0_dp) > tol .or. &
          abs(v(3) - 300.0_dp) > tol .or. abs(v(4)) > tol .or. &
          abs(v(5)) > tol .or. abs(v(6)) > tol) then
         print *, 'FAIL  test_map2stress_known: got', v
         nfail = nfail + 1
      end if
   end subroutine test_map2stress_known

   subroutine test_map2stran_known()
      ! map2stran on a diagonal tensor gives the diagonal elements (no factor-of-2 on zeros)
      real(dp) :: T(3,3), v(6)
      T = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,2.0_dp,0.0_dp, 0.0_dp,0.0_dp,3.0_dp], [3,3])
      v = map2stran(T, 6)
      if (abs(v(1) - 1.0_dp) > tol .or. abs(v(2) - 2.0_dp) > tol .or. &
          abs(v(3) - 3.0_dp) > tol .or. abs(v(4)) > tol .or. &
          abs(v(5)) > tol .or. abs(v(6)) > tol) then
         print *, 'FAIL  test_map2stran_known: got', v
         nfail = nfail + 1
      end if
   end subroutine test_map2stran_known

end program test_maps
