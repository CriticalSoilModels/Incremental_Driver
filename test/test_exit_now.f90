! Regression tests for EXITNOW.
! EXITNOW evaluates a string condition (e.g. 's1 > 100.0') against the
! current stress, strain, and state variable arrays.
program test_exit_now
   use indr_parser, only: EXITNOW
   implicit none

   integer, parameter :: nstatv = 10
   real(8) :: stress(6), stran(6), statev(nstatv)
   character(len=40) :: cond
   logical :: result
   integer :: n_fail

   n_fail = 0

   ! Baseline: zero state
   stress = 0.0d0
   stran  = 0.0d0
   statev = 0.0d0

   ! --- Test 1: stress component condition — true branch ---
   stress(1) = 150.0d0
   cond   = 's1 > 100.0'
   result = EXITNOW(cond, stress, stran, statev, nstatv)
   call assert_true(result, 's1=150 > 100 should be .true.', n_fail)

   ! --- Test 2: stress component condition — false branch ---
   stress(1) = 50.0d0
   cond   = 's1 > 100.0'
   result = EXITNOW(cond, stress, stran, statev, nstatv)
   call assert_true(.not. result, 's1=50 > 100 should be .false.', n_fail)

   ! --- Test 3: strain component, less-than condition ---
   stress = 0.0d0
   stran(1) = 0.001d0
   cond   = 'e1 < 0.005'
   result = EXITNOW(cond, stress, stran, statev, nstatv)
   call assert_true(result, 'e1=0.001 < 0.005 should be .true.', n_fail)

   ! --- Test 4: state variable component ---
   stran = 0.0d0
   statev(2) = 5.0d0
   cond   = 'v2 > 10.0'
   result = EXITNOW(cond, stress, stran, statev, nstatv)
   call assert_true(.not. result, 'v2=5 > 10 should be .false.', n_fail)

   ! --- Test 5: condition with no comparison operator → always .false. ---
   ! EXITNOW reaches goto 555 when neither '<' nor '>' is found.
   statev = 0.0d0
   cond   = 'no_operator_here'
   result = EXITNOW(cond, stress, stran, statev, nstatv)
   call assert_true(.not. result, 'bad cond string should return .false.', n_fail)

   ! --- Report ---
   if (n_fail == 0) then
      print *, 'PASS  test_exit_now'
   else
      print *, 'FAIL  test_exit_now:', n_fail, 'assertion(s) failed'
      stop 1
   end if

contains

   subroutine assert_true(cond_val, label, n_fail)
      logical,      intent(in)    :: cond_val
      character(*), intent(in)    :: label
      integer,      intent(inout) :: n_fail
      if (.not. cond_val) then
         print *, 'FAIL  ', trim(label)
         n_fail = n_fail + 1
      end if
   end subroutine

end program test_exit_now
