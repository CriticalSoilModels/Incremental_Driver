program test_material_state
   !! Tests for material_state_t: verify fields are readable/writable,
   !! that allocatable statev round-trips correctly, and that struct
   !! assignment deep-copies the allocatable component.
   use stdlib_kinds, only: dp
   use indr_step_params, only: material_state_t
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-14_dp

   call test_field_assignment()
   call test_statev_alloc()
   call test_struct_copy()

   if (nfail == 0) then
      print *, 'PASS  test_material_state'
   else
      print *, 'FAIL  test_material_state — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   subroutine test_field_assignment()
      !! Verify scalar and fixed-size array fields can be set and read back.
      type(material_state_t) :: s
      real(dp), parameter :: identity33(3,3) = reshape( &
         [1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])

      s%sig      = [100.0_dp, 50.0_dp, 50.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      s%eps      = [0.01_dp, -0.005_dp, -0.005_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      s%time     = [0.5_dp, 10.0_dp]
      s%dt       = 0.1_dp
      s%temp     = 20.0_dp
      s%F_start  = identity33
      s%F_end    = identity33

      if (abs(s%sig(1) - 100.0_dp) > tol) then
         print *, 'FAIL  test_field_assignment: sig(1)'; nfail = nfail + 1
      end if
      if (abs(s%eps(1) - 0.01_dp) > tol) then
         print *, 'FAIL  test_field_assignment: eps(1)'; nfail = nfail + 1
      end if
      if (abs(s%time(1) - 0.5_dp) > tol) then
         print *, 'FAIL  test_field_assignment: time(1)'; nfail = nfail + 1
      end if
      if (abs(s%time(2) - 10.0_dp) > tol) then
         print *, 'FAIL  test_field_assignment: time(2)'; nfail = nfail + 1
      end if
      if (abs(s%dt - 0.1_dp) > tol) then
         print *, 'FAIL  test_field_assignment: dt'; nfail = nfail + 1
      end if
      if (abs(s%temp - 20.0_dp) > tol) then
         print *, 'FAIL  test_field_assignment: temp'; nfail = nfail + 1
      end if
      if (abs(s%F_start(1,1) - 1.0_dp) > tol .or. abs(s%F_start(1,2)) > tol) then
         print *, 'FAIL  test_field_assignment: F_start'; nfail = nfail + 1
      end if
   end subroutine test_field_assignment

   subroutine test_statev_alloc()
      !! Verify allocatable statev can be allocated and read back.
      type(material_state_t) :: s
      integer :: i

      allocate(s%statev(5))
      do i = 1, 5
         s%statev(i) = real(i, dp) * 0.1_dp
      end do

      if (size(s%statev) /= 5) then
         print *, 'FAIL  test_statev_alloc: size'; nfail = nfail + 1
      end if
      if (abs(s%statev(3) - 0.3_dp) > tol) then
         print *, 'FAIL  test_statev_alloc: statev(3)'; nfail = nfail + 1
      end if
   end subroutine test_statev_alloc

   subroutine test_struct_copy()
      !! Verify intrinsic assignment deep-copies the allocatable statev component.
      type(material_state_t) :: s1, s2

      s1%sig   = 0.0_dp
      s1%eps   = 0.0_dp
      s1%time  = [1.0_dp, 5.0_dp]
      s1%dt    = 0.5_dp
      s1%temp  = 15.0_dp
      s1%F_start = 0.0_dp;  s1%F_start(1,1) = 1.0_dp
      s1%F_start(2,2) = 1.0_dp;  s1%F_start(3,3) = 1.0_dp
      s1%F_end = s1%F_start
      allocate(s1%statev(3), source=[1.0_dp, 2.0_dp, 3.0_dp])

      s2 = s1

      if (.not. allocated(s2%statev)) then
         print *, 'FAIL  test_struct_copy: statev not allocated after copy'; nfail = nfail + 1
         return
      end if
      if (size(s2%statev) /= 3) then
         print *, 'FAIL  test_struct_copy: statev size'; nfail = nfail + 1
      end if
      if (abs(s2%statev(2) - 2.0_dp) > tol) then
         print *, 'FAIL  test_struct_copy: statev value'; nfail = nfail + 1
      end if

      ! Mutate s2%statev and verify s1 is unaffected (deep copy, not alias)
      s2%statev(1) = 999.0_dp
      if (abs(s1%statev(1) - 1.0_dp) > tol) then
         print *, 'FAIL  test_struct_copy: statev is aliased (shallow copy bug!)'; nfail = nfail + 1
      end if
   end subroutine test_struct_copy

end program test_material_state
