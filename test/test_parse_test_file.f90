program test_parse_test_file
   !! Tests for parse_test_file in indr_test_parser.
   !! Each test writes a minimal test.inp to a scratch file, parses it,
   !! and checks the returned step records against expected values.
   use stdlib_kinds,       only: dp
   use indr_test_parser,   only: step_record_t, parse_test_file
   use indr_types,         only: StressAlignment
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-12_dp

   call test_single_linear_load()
   call test_heading_is_read()
   call test_repetition_expansion()

   if (nfail == 0) then
      print *, 'PASS  test_parse_test_file'
   else
      print *, 'FAIL  test_parse_test_file — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   ! ---------------------------------------------------------------------------
   ! test_single_linear_load
   ! Minimal file with one *LinearLoad step; checks count, fields, write_freq.
   ! ---------------------------------------------------------------------------
   subroutine test_single_linear_load()
      integer                          :: fid, n_steps
      type(step_record_t), allocatable :: steps(:)
      type(StressAlignment)            :: align
      character(len=260)               :: heading

      open(newunit=fid, status='scratch')
      write(fid, '(a)') 'output.txt # test heading'
      write(fid, '(a)') '*LinearLoad'
      write(fid, '(a)') '10 5 1.0 : 5'
      write(fid, '(a)') '*Cartesian'
      write(fid, '(a)') '0 0.001'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '*End'
      rewind(fid)

      call parse_test_file(fid, heading, steps, n_steps, align)
      close(fid)

      call check_int(n_steps, 1, 'single_linear: n_steps', nfail)

      if (n_steps < 1) return

      call check_int(steps(1)%config%n_inc, 10, 'single_linear: n_inc', nfail)
      call check_int(steps(1)%config%max_iter, 5, 'single_linear: max_iter', nfail)
      call check_real(steps(1)%config%delta_time, 1.0_dp, tol, 'single_linear: delta_time', nfail)
      call check_int(steps(1)%write_freq, 5, 'single_linear: write_freq', nfail)

      if (trim(steps(1)%config%load_type) /= '*LinearLoad') then
         print *, 'FAIL  single_linear: load_type expected "*LinearLoad" got "', &
                  trim(steps(1)%config%load_type), '"'
         nfail = nfail + 1
      end if

      if (trim(steps(1)%config%coord_sys) /= '*Cartesian') then
         print *, 'FAIL  single_linear: coord_sys expected "*Cartesian" got "', &
                  trim(steps(1)%config%coord_sys), '"'
         nfail = nfail + 1
      end if

      call check_int(steps(1)%config%ifstress(1), 0, 'single_linear: ifstress(1)', nfail)
      call check_real(steps(1)%config%delta_load(1), 0.001_dp, tol, &
                      'single_linear: delta_load(1)', nfail)
      call check_real(steps(1)%config%delta_load(2), 0.0_dp, tol, &
                      'single_linear: delta_load(2)', nfail)

      ! align should be inactive when no *ImportFile is present
      if (align%active) then
         print *, 'FAIL  single_linear: align%active should be .false.'
         nfail = nfail + 1
      end if
   end subroutine test_single_linear_load

   ! ---------------------------------------------------------------------------
   ! test_heading_is_read
   ! Checks that the first line is captured as heading and multiple steps parse.
   ! ---------------------------------------------------------------------------
   subroutine test_heading_is_read()
      integer                          :: fid, n_steps
      type(step_record_t), allocatable :: steps(:)
      type(StressAlignment)            :: align
      character(len=260)               :: heading

      open(newunit=fid, status='scratch')
      write(fid, '(a)') 'myresults.txt # two-step test'
      write(fid, '(a)') '*LinearLoad'
      write(fid, '(a)') '5 10 0.5 : 1'
      write(fid, '(a)') '*Cartesian'
      write(fid, '(a)') '1 100.0'
      write(fid, '(a)') '1 100.0'
      write(fid, '(a)') '1 100.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '*LinearLoad'
      write(fid, '(a)') '20 50 2.0 : 4'
      write(fid, '(a)') '*Roscoe'
      write(fid, '(a)') '0 0.002'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '*End'
      rewind(fid)

      call parse_test_file(fid, heading, steps, n_steps, align)
      close(fid)

      call check_int(n_steps, 2, 'two_step: n_steps', nfail)

      if (index(trim(heading), 'myresults.txt') == 0) then
         print *, 'FAIL  two_step: heading does not contain "myresults.txt", got: "', &
                  trim(heading), '"'
         nfail = nfail + 1
      end if

      if (n_steps < 2) return

      call check_int(steps(1)%config%n_inc, 5, 'two_step: step1 n_inc', nfail)
      call check_int(steps(1)%write_freq, 1,   'two_step: step1 write_freq', nfail)
      call check_int(steps(2)%config%n_inc, 20, 'two_step: step2 n_inc', nfail)
      call check_int(steps(2)%write_freq, 4,    'two_step: step2 write_freq', nfail)
      call check_real(steps(2)%config%delta_time, 2.0_dp, tol, 'two_step: step2 delta_time', nfail)
   end subroutine test_heading_is_read

   ! ---------------------------------------------------------------------------
   ! test_repetition_expansion
   ! *Repetition with 2 steps x 3 repetitions = 6 total steps.
   ! ---------------------------------------------------------------------------
   subroutine test_repetition_expansion()
      integer                          :: fid, n_steps
      type(step_record_t), allocatable :: steps(:)
      type(StressAlignment)            :: align
      character(len=260)               :: heading

      open(newunit=fid, status='scratch')
      write(fid, '(a)') 'output.txt'
      write(fid, '(a)') '*Repetition'
      write(fid, '(a)') '2 3'           ! 2 steps, 3 repetitions
      write(fid, '(a)') '*LinearLoad'
      write(fid, '(a)') '10 5 1.0 : 1'
      write(fid, '(a)') '*Cartesian'
      write(fid, '(a)') '0 0.001'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '*LinearLoad'
      write(fid, '(a)') '10 5 1.0 : 1'
      write(fid, '(a)') '*Cartesian'
      write(fid, '(a)') '0 -0.001'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '0 0.0'
      write(fid, '(a)') '*End'
      rewind(fid)

      call parse_test_file(fid, heading, steps, n_steps, align)
      close(fid)

      call check_int(n_steps, 6, 'repetition: n_steps (2 x 3)', nfail)

      if (n_steps < 6) return

      ! Odd steps (1,3,5) should have positive delta_load(1)
      call check_real(steps(1)%config%delta_load(1),  0.001_dp, tol, 'repetition: step1 delta_load(1)', nfail)
      call check_real(steps(3)%config%delta_load(1),  0.001_dp, tol, 'repetition: step3 delta_load(1)', nfail)
      call check_real(steps(5)%config%delta_load(1),  0.001_dp, tol, 'repetition: step5 delta_load(1)', nfail)

      ! Even steps (2,4,6) should have negative delta_load(1)
      call check_real(steps(2)%config%delta_load(1), -0.001_dp, tol, 'repetition: step2 delta_load(1)', nfail)
      call check_real(steps(4)%config%delta_load(1), -0.001_dp, tol, 'repetition: step4 delta_load(1)', nfail)
      call check_real(steps(6)%config%delta_load(1), -0.001_dp, tol, 'repetition: step6 delta_load(1)', nfail)

      ! All steps should have n_inc = 10
      call check_int(steps(1)%config%n_inc, 10, 'repetition: step1 n_inc', nfail)
      call check_int(steps(6)%config%n_inc, 10, 'repetition: step6 n_inc', nfail)
   end subroutine test_repetition_expansion

   ! ---------------------------------------------------------------------------
   ! Assertion helpers
   ! ---------------------------------------------------------------------------

   subroutine check_int(got, expected, label, nfail)
      integer,      intent(in)    :: got, expected
      character(*), intent(in)    :: label
      integer,      intent(inout) :: nfail
      if (got /= expected) then
         print *, 'FAIL  ', trim(label), ': got', got, 'expected', expected
         nfail = nfail + 1
      end if
   end subroutine check_int

   subroutine check_real(got, expected, tol, label, nfail)
      real(dp),     intent(in)    :: got, expected, tol
      character(*), intent(in)    :: label
      integer,      intent(inout) :: nfail
      if (abs(got - expected) > tol) then
         print *, 'FAIL  ', trim(label), ': got', got, 'expected', expected
         nfail = nfail + 1
      end if
   end subroutine check_real

end program test_parse_test_file
