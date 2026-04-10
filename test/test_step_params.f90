program test_step_params
   ! Tests for step_config_t: verify fields are readable/writable and struct
   ! assignment copies all data correctly.
   use stdlib_kinds, only: dp
   use indr_step_params, only: step_config_t
   use indr_constants, only: voigt_len
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-14_dp

   call test_field_assignment()
   call test_struct_copy()

   if (nfail == 0) then
      print *, 'PASS  test_step_params'
   else
      print *, 'FAIL  test_step_params — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   subroutine test_field_assignment()
      ! Verify all fields can be set and read back.
      type(step_config_t) :: cfg
      integer :: i

      cfg%n_inc           = 50
      cfg%max_iter        = 25
      cfg%ifstress        = [1, 0, 1, 0, 1, 0]
      cfg%delta_load      = [(real(i,dp)*10.0_dp, i=1,9)]
      cfg%delta_load_circ = [(real(i,dp)*1.0_dp,  i=1,6)]
      cfg%phase0          = [(real(i,dp)*0.1_dp,  i=1,6)]
      cfg%delta_time      = 2.5_dp
      cfg%delta_temp      = 10.0_dp
      cfg%load_type       = '*LinearLoad'
      cfg%coord_sys       = '*Roscoe'
      cfg%exit_cond       = '*StressCondition'
      cfg%has_exit_cond   = .true.
      cfg%import_file     = 'data.dat'
      cfg%n_import        = 4
      cfg%columns_in_file = [1,2,3,4,5,6,7]
      cfg%import_factor   = [(real(i,dp)*0.5_dp, i=1,7)]
      cfg%dfgrd0          = reshape([(real(i,dp)*0.1_dp, i=1,9)], [3,3])
      cfg%dfgrd1          = reshape([(real(i,dp)*0.2_dp, i=1,9)], [3,3])
      cfg%cMt             = 0.0_dp
      cfg%cMe             = 0.0_dp
      do i = 1, 6
         cfg%cMt(i,i) = real(i, dp)
         cfg%cMe(i,i) = real(i, dp) * 2.0_dp
      end do
      cfg%mbinc           = [(real(i,dp)*3.0_dp, i=1,6)]

      if (cfg%n_inc /= 50) then
         print *, 'FAIL  test_field_assignment: n_inc'; nfail = nfail + 1
      end if
      if (cfg%max_iter /= 25) then
         print *, 'FAIL  test_field_assignment: max_iter'; nfail = nfail + 1
      end if
      if (cfg%ifstress(1) /= 1 .or. cfg%ifstress(2) /= 0) then
         print *, 'FAIL  test_field_assignment: ifstress'; nfail = nfail + 1
      end if
      if (abs(cfg%delta_time - 2.5_dp) > tol) then
         print *, 'FAIL  test_field_assignment: delta_time'; nfail = nfail + 1
      end if
      if (trim(cfg%load_type) /= '*LinearLoad') then
         print *, 'FAIL  test_field_assignment: load_type'; nfail = nfail + 1
      end if
      if (.not. cfg%has_exit_cond) then
         print *, 'FAIL  test_field_assignment: has_exit_cond'; nfail = nfail + 1
      end if
      if (cfg%n_import /= 4) then
         print *, 'FAIL  test_field_assignment: n_import'; nfail = nfail + 1
      end if
   end subroutine test_field_assignment

   subroutine test_struct_copy()
      ! Verify that cfg2 = cfg1 copies all data (no hidden aliasing).
      type(step_config_t) :: cfg1, cfg2

      cfg1%n_inc     = 100
      cfg1%delta_time = 5.0_dp
      cfg1%load_type = '*CirculatingLoad'
      cfg1%ifstress  = [1,1,0,0,0,0]

      cfg2 = cfg1

      if (cfg2%n_inc /= 100) then
         print *, 'FAIL  test_struct_copy: n_inc'; nfail = nfail + 1
      end if
      if (abs(cfg2%delta_time - 5.0_dp) > tol) then
         print *, 'FAIL  test_struct_copy: delta_time'; nfail = nfail + 1
      end if
      if (trim(cfg2%load_type) /= '*CirculatingLoad') then
         print *, 'FAIL  test_struct_copy: load_type'; nfail = nfail + 1
      end if
      if (cfg2%ifstress(3) /= 0) then
         print *, 'FAIL  test_struct_copy: ifstress'; nfail = nfail + 1
      end if

      ! Mutate cfg2 and verify cfg1 is unaffected
      cfg2%n_inc = 999
      if (cfg1%n_inc /= 100) then
         print *, 'FAIL  test_struct_copy: struct copy is aliased!'; nfail = nfail + 1
      end if
   end subroutine test_struct_copy

end program test_step_params
