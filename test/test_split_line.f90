! Regression tests for splitaLine.
! Each test is independent; n_fail accumulates failures.
! Exit code 1 on any failure so fpm test reports FAILED.
program test_split_line
   use indr_parser, only: splitaLine
   implicit none

   character(len=40) :: aline, left, right
   logical :: ok
   integer :: n_fail

   n_fail = 0

   ! --- Test 1: separator present in the middle ---
   aline = 'alpha:beta'
   call splitaLine(aline, ':', left, right, ok)
   call assert_true(ok,                      'present: ok is .true.',  n_fail)
   call assert_str(trim(left),  'alpha',     'present: left',          n_fail)
   call assert_str(trim(right), 'beta',      'present: right',         n_fail)

   ! --- Test 2: separator absent ---
   aline = 'alpha'
   ok = .true.   ! pre-set opposite to prove it gets cleared
   call splitaLine(aline, ':', left, right, ok)
   call assert_true(.not. ok,           'absent: ok is .false.',        n_fail)
   call assert_str(trim(left), 'alpha', 'absent: left is the full line', n_fail)

   ! --- Test 3: separator is the first character ---
   aline = ':beta'
   call splitaLine(aline, ':', left, right, ok)
   call assert_true(ok,                 'first: ok is .true.', n_fail)
   call assert_str(trim(left),  '',     'first: left is empty', n_fail)
   call assert_str(trim(right), 'beta', 'first: right',         n_fail)

   ! --- Test 4: separator is the last character ---
   aline = 'alpha:'
   call splitaLine(aline, ':', left, right, ok)
   call assert_true(ok,                 'last: ok is .true.', n_fail)
   call assert_str(trim(left),  'alpha', 'last: left',          n_fail)
   call assert_str(trim(right), '',      'last: right is empty', n_fail)

   ! --- Report ---
   if (n_fail == 0) then
      print *, 'PASS  test_split_line'
   else
      print *, 'FAIL  test_split_line:', n_fail, 'assertion(s) failed'
      stop 1
   end if

contains

   subroutine assert_true(cond, label, n_fail)
      logical,          intent(in)    :: cond
      character(*),     intent(in)    :: label
      integer,          intent(inout) :: n_fail
      if (.not. cond) then
         print *, 'FAIL  ', trim(label)
         n_fail = n_fail + 1
      end if
   end subroutine

   subroutine assert_str(got, expected, label, n_fail)
      character(*), intent(in)    :: got, expected, label
      integer,      intent(inout) :: n_fail
      if (got /= expected) then
         print *, 'FAIL  ', trim(label), ': got "', got, '" expected "', expected, '"'
         n_fail = n_fail + 1
      end if
   end subroutine

end program test_split_line
