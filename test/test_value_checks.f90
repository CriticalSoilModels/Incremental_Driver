program test_value_checks
   use stdlib_kinds, only: dp
   use indr_value_checks, only: set_zero_with_tol
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-14_dp

   call test_values_above_default_tol()
   call test_values_below_default_tol()
   call test_custom_tolerance()
   call test_mixed_array()

   if (nfail == 0) then
      print *, 'PASS  test_value_checks'
   else
      print *, 'FAIL  test_value_checks — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   subroutine test_values_above_default_tol()
      ! Values well above default tolerance (1e-99) should be unchanged
      real(dp) :: arr(3), result(3)
      integer :: i
      arr = [1.0_dp, -2.5_dp, 3.14e-50_dp]
      result = set_zero_with_tol(arr)
      do i = 1, 3
         if (abs(result(i) - arr(i)) > tol) then
            print *, 'FAIL  test_values_above_default_tol: index', i, &
                     'expected', arr(i), 'got', result(i)
            nfail = nfail + 1
         end if
      end do
   end subroutine test_values_above_default_tol

   subroutine test_values_below_default_tol()
      ! Values below default tolerance (1e-99) should be zeroed
      real(dp) :: arr(2), result(2)
      ! 1e-300 is below the default 1e-99 threshold
      arr = [1.0e-300_dp, -1.0e-200_dp]
      result = set_zero_with_tol(arr)
      if (abs(result(1)) > tol) then
         print *, 'FAIL  test_values_below_default_tol: index 1 expected 0 got', result(1)
         nfail = nfail + 1
      end if
      if (abs(result(2)) > tol) then
         print *, 'FAIL  test_values_below_default_tol: index 2 expected 0 got', result(2)
         nfail = nfail + 1
      end if
   end subroutine test_values_below_default_tol

   subroutine test_custom_tolerance()
      ! Custom tolerance of 1e-5: values below it should be zeroed
      real(dp) :: arr(3), result(3)
      real(dp), parameter :: custom_tol = 1.0e-5_dp
      arr = [1.0_dp, 1.0e-6_dp, -1.0e-10_dp]
      result = set_zero_with_tol(arr, ext_tol=custom_tol)
      ! arr(1)=1.0 is above threshold -> unchanged
      if (abs(result(1) - 1.0_dp) > tol) then
         print *, 'FAIL  test_custom_tolerance: index 1 expected 1.0 got', result(1)
         nfail = nfail + 1
      end if
      ! arr(2)=1e-6 is below 1e-5 threshold -> zero
      if (abs(result(2)) > tol) then
         print *, 'FAIL  test_custom_tolerance: index 2 expected 0 got', result(2)
         nfail = nfail + 1
      end if
      ! arr(3)=1e-10 is below 1e-5 threshold -> zero
      if (abs(result(3)) > tol) then
         print *, 'FAIL  test_custom_tolerance: index 3 expected 0 got', result(3)
         nfail = nfail + 1
      end if
   end subroutine test_custom_tolerance

   subroutine test_mixed_array()
      ! Mixed: some values zeroed, some unchanged
      real(dp) :: arr(4), result(4)
      arr = [100.0_dp, 1.0e-300_dp, -50.0_dp, 1.0e-200_dp]
      result = set_zero_with_tol(arr)
      if (abs(result(1) - 100.0_dp) > tol) then
         print *, 'FAIL  test_mixed_array: index 1 expected 100 got', result(1)
         nfail = nfail + 1
      end if
      if (abs(result(2)) > tol) then
         print *, 'FAIL  test_mixed_array: index 2 expected 0 got', result(2)
         nfail = nfail + 1
      end if
      if (abs(result(3) - (-50.0_dp)) > tol) then
         print *, 'FAIL  test_mixed_array: index 3 expected -50 got', result(3)
         nfail = nfail + 1
      end if
      if (abs(result(4)) > tol) then
         print *, 'FAIL  test_mixed_array: index 4 expected 0 got', result(4)
         nfail = nfail + 1
      end if
   end subroutine test_mixed_array

end program test_value_checks
