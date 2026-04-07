! Regression tests for get_increment — *LinearLoad branch.
!
! get_increment converts a step description (load type, total deltaLoad, ninc)
! into a per-increment ddstress/dstran pair. For *LinearLoad the formula is
! trivially linear, so expected values are exact.
program test_get_increment
   use mod_inc_driver_funcs, only: get_increment
   implicit none

   character(len=40) :: keywords(10)
   integer  :: ifstress(6), ninc
   real(8)  :: time(2), deltaTime
   real(8)  :: deltaLoadCirc(6), phase0(6), deltaLoad(9), deltaTemp
   real(8)  :: dtime, ddstress(6), dstran(6), dTemp
   real(8)  :: Qb33(3,3), dfgrd0(3,3), dfgrd1(3,3), drot(3,3)

   real(8), parameter :: tol = 1.0d-10
   integer :: n_fail

   n_fail = 0

   ! Fixed inputs that don't vary between tests
   keywords      = ''
   time          = 0.0d0
   deltaLoadCirc = 0.0d0
   phase0        = 0.0d0
   deltaTemp     = 0.0d0

   ! Identity matrices for the intent(inout) deformation-gradient args
   dfgrd0 = 0.0d0;  dfgrd0(1,1)=1.0d0; dfgrd0(2,2)=1.0d0; dfgrd0(3,3)=1.0d0
   dfgrd1 = dfgrd0
   drot   = dfgrd0

   ! ===========================================================================
   ! Test 1: *LinearLoad — pure stress control
   !   deltaLoad(1)=60, ninc=10 → ddstress(1) = 6.0, dstran = 0
   ! ===========================================================================
   keywords(2) = '*LinearLoad'
   ninc        = 10
   deltaTime   = 1.0d0
   ifstress    = [1, 1, 1, 1, 1, 1]
   deltaLoad   = [60.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]

   call get_increment(keywords, time, deltaTime, ifstress, ninc, &
                      deltaLoadCirc, phase0, deltaLoad, deltaTemp, &
                      dtime, ddstress, dstran, dTemp, Qb33, dfgrd0, dfgrd1, drot)

   call assert_close(dtime,       0.1d0,  tol, 'stress_ctrl: dtime',       n_fail)
   call assert_close(ddstress(1), 6.0d0,  tol, 'stress_ctrl: ddstress(1)', n_fail)
   call assert_close(ddstress(2), 0.0d0,  tol, 'stress_ctrl: ddstress(2)', n_fail)
   call assert_close(dstran(1),   0.0d0,  tol, 'stress_ctrl: dstran(1)',   n_fail)
   call assert_close(dstran(2),   0.0d0,  tol, 'stress_ctrl: dstran(2)',   n_fail)

   ! ===========================================================================
   ! Test 2: *LinearLoad — pure strain control
   !   deltaLoad(1)=0.01, ninc=10 → dstran(1) = 0.001, ddstress = 0
   ! ===========================================================================
   keywords(2) = '*LinearLoad'
   ninc        = 10
   deltaTime   = 1.0d0
   ifstress    = [0, 0, 0, 0, 0, 0]
   deltaLoad   = [0.01d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
   dfgrd0 = 0.0d0;  dfgrd0(1,1)=1.0d0; dfgrd0(2,2)=1.0d0; dfgrd0(3,3)=1.0d0
   dfgrd1 = dfgrd0; drot = dfgrd0

   call get_increment(keywords, time, deltaTime, ifstress, ninc, &
                      deltaLoadCirc, phase0, deltaLoad, deltaTemp, &
                      dtime, ddstress, dstran, dTemp, Qb33, dfgrd0, dfgrd1, drot)

   call assert_close(dtime,        0.1d0,   tol, 'strain_ctrl: dtime',        n_fail)
   call assert_close(dstran(1),    0.001d0, tol, 'strain_ctrl: dstran(1)',    n_fail)
   call assert_close(dstran(2),    0.0d0,   tol, 'strain_ctrl: dstran(2)',    n_fail)
   call assert_close(ddstress(1),  0.0d0,   tol, 'strain_ctrl: ddstress(1)', n_fail)

   ! ===========================================================================
   ! Test 3: *LinearLoad — mixed control
   !   Component 1: stress-ctrl, load=100, ninc=20 → ddstress(1) = 5.0
   !   Component 2: strain-ctrl, load=0.005, ninc=20 → dstran(2) = 0.00025
   !   dtime = 2.0/20 = 0.1
   ! ===========================================================================
   keywords(2) = '*LinearLoad'
   ninc        = 20
   deltaTime   = 2.0d0
   ifstress    = [1, 0, 0, 0, 0, 0]
   deltaLoad   = [100.0d0, 0.005d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
   dfgrd0 = 0.0d0;  dfgrd0(1,1)=1.0d0; dfgrd0(2,2)=1.0d0; dfgrd0(3,3)=1.0d0
   dfgrd1 = dfgrd0; drot = dfgrd0

   call get_increment(keywords, time, deltaTime, ifstress, ninc, &
                      deltaLoadCirc, phase0, deltaLoad, deltaTemp, &
                      dtime, ddstress, dstran, dTemp, Qb33, dfgrd0, dfgrd1, drot)

   call assert_close(dtime,       0.1d0,     tol, 'mixed: dtime',       n_fail)
   call assert_close(ddstress(1), 5.0d0,     tol, 'mixed: ddstress(1)', n_fail)
   call assert_close(dstran(1),   0.0d0,     tol, 'mixed: dstran(1)',   n_fail)
   call assert_close(dstran(2),   0.00025d0, tol, 'mixed: dstran(2)',   n_fail)
   call assert_close(ddstress(2), 0.0d0,     tol, 'mixed: ddstress(2)', n_fail)

   ! --- Report ---
   if (n_fail == 0) then
      print *, 'PASS  test_get_increment'
   else
      print *, 'FAIL  test_get_increment:', n_fail, 'assertion(s) failed'
      stop 1
   end if

contains

   subroutine assert_close(got, expected, tol, label, n_fail)
      real(8),      intent(in)    :: got, expected, tol
      character(*), intent(in)    :: label
      integer,      intent(inout) :: n_fail
      ! Relative tolerance for large values, absolute for near-zero
      if (abs(got - expected) > tol * max(abs(expected), 1.0d0)) then
         print *, 'FAIL  ', trim(label), ': got', got, 'expected', expected
         n_fail = n_fail + 1
      end if
   end subroutine

end program test_get_increment
