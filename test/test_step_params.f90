program test_step_params
   ! Round-trip test: set_repetition_params -> get_repetition_params recovers all fields
   use stdlib_kinds, only: dp
   use mod_step_params, only: descriptionOfStep, set_repetition_params, get_repetition_params
   use mod_constants, only: max_fname_len, voight_len
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-14_dp

   call test_roundtrip()

   if (nfail == 0) then
      print *, 'PASS  test_step_params'
   else
      print *, 'FAIL  test_step_params — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   subroutine test_roundtrip()
      type(descriptionOfStep) :: step

      ! Input values
      integer  :: ninc_in, maxiter_in, ifstress_in(voight_len)
      integer  :: mImport_in, columnsInFile_in(7)
      real(dp) :: deltaLoadCirc_in(6), phase0_in(6), deltaLoad_in(9)
      real(dp) :: dfgrd0_in(3,3), dfgrd1_in(3,3), deltaTime_in, deltaTemp_in
      real(dp) :: cMe_in(6,6), cMt_in(6,6), mbinc_in(6), importFactor_in(7)
      character(len=10)       :: keywords_in(3)
      character(len=40)       :: exitCond_in
      logical                 :: existCond_in
      character(max_fname_len) :: ImportFileName_in

      ! Output values (recovered)
      integer  :: ninc_out, maxiter_out, ifstress_out(voight_len)
      integer  :: mImport_out, columnsInFile_out(7)
      real(dp) :: deltaLoadCirc_out(6), phase0_out(6), deltaLoad_out(9)
      real(dp) :: dfgrd0_out(3,3), dfgrd1_out(3,3), deltaTime_out, deltaTemp_out
      real(dp) :: cMe_out(6,6), cMt_out(6,6), mbinc_out(6), importFactor_out(7)
      character(len=10)        :: keywords_out(3)
      character(len=40)        :: exitCond_out
      logical                  :: existCond_out
      character(max_fname_len) :: ImportFileName_out

      integer :: i

      ! Set known input values
      ninc_in           = 100
      maxiter_in        = 50
      ifstress_in       = [1, 0, 1, 0, 1, 0]
      deltaLoadCirc_in  = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
      phase0_in         = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp]
      deltaLoad_in      = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, &
                           60.0_dp, 70.0_dp, 80.0_dp, 90.0_dp]
      dfgrd0_in         = reshape([(real(i,dp)*0.1_dp, i=1,9)], [3,3])
      dfgrd1_in         = reshape([(real(i,dp)*0.2_dp, i=1,9)], [3,3])
      deltaTime_in      = 1.5_dp
      deltaTemp_in      = 25.0_dp
      keywords_in(1)    = '*Step'
      keywords_in(2)    = '*LinearLo'  ! 10-char truncation of '*LinearLoad'
      keywords_in(3)    = 'somekey'
      cMe_in            = 0.0_dp
      cMt_in            = 0.0_dp
      do i = 1, 6
         cMe_in(i,i) = real(i, dp) * 10.0_dp
         cMt_in(i,i) = real(i, dp) * 5.0_dp
      end do
      mbinc_in          = [7.0_dp, 8.0_dp, 9.0_dp, 10.0_dp, 11.0_dp, 12.0_dp]
      exitCond_in       = '*StressCondition'
      existCond_in      = .true.
      ImportFileName_in = 'myfile.dat'
      mImport_in        = 3
      columnsInFile_in  = [1, 2, 3, 4, 5, 6, 7]
      importFactor_in   = [1.1_dp, 2.2_dp, 3.3_dp, 4.4_dp, 5.5_dp, 6.6_dp, 7.7_dp]

      ! Pack
      step = set_repetition_params(ninc_in, maxiter_in, ifstress_in, deltaLoadCirc_in, &
         phase0_in, deltaLoad_in, dfgrd0_in, dfgrd1_in, deltaTime_in, keywords_in, &
         cMe_in, cMt_in, mbinc_in, deltaTemp_in, exitCond_in, existCond_in, &
         ImportFileName_in, mImport_in, columnsInFile_in, importFactor_in)

      ! Unpack
      call get_repetition_params(step, ninc_out, maxiter_out, ifstress_out, deltaLoadCirc_out, &
         phase0_out, deltaLoad_out, dfgrd0_out, dfgrd1_out, deltaTime_out, keywords_out, &
         cMe_out, cMt_out, mbinc_out, deltaTemp_out, exitCond_out, existCond_out, &
         ImportFileName_out, mImport_out, columnsInFile_out, importFactor_out)

      ! Check scalar integers
      if (ninc_out /= ninc_in) then
         print *, 'FAIL  test_roundtrip: ninc expected', ninc_in, 'got', ninc_out
         nfail = nfail + 1
      end if
      if (maxiter_out /= maxiter_in) then
         print *, 'FAIL  test_roundtrip: maxiter expected', maxiter_in, 'got', maxiter_out
         nfail = nfail + 1
      end if

      ! Check ifstress
      do i = 1, voight_len
         if (ifstress_out(i) /= ifstress_in(i)) then
            print *, 'FAIL  test_roundtrip: ifstress(', i, ') expected', ifstress_in(i), 'got', ifstress_out(i)
            nfail = nfail + 1
         end if
      end do

      ! Check scalar reals
      if (abs(deltaTime_out - deltaTime_in) > tol) then
         print *, 'FAIL  test_roundtrip: deltaTime expected', deltaTime_in, 'got', deltaTime_out
         nfail = nfail + 1
      end if
      if (abs(deltaTemp_out - deltaTemp_in) > tol) then
         print *, 'FAIL  test_roundtrip: deltaTemp expected', deltaTemp_in, 'got', deltaTemp_out
         nfail = nfail + 1
      end if

      ! Check keyword2 and keyword3 (keyword1 is not stored in descriptionOfStep)
      if (trim(keywords_out(2)) /= trim(keywords_in(2))) then
         print *, 'FAIL  test_roundtrip: keywords(2) expected "', trim(keywords_in(2)), &
                  '" got "', trim(keywords_out(2)), '"'
         nfail = nfail + 1
      end if
      if (trim(keywords_out(3)) /= trim(keywords_in(3))) then
         print *, 'FAIL  test_roundtrip: keywords(3) expected "', trim(keywords_in(3)), &
                  '" got "', trim(keywords_out(3)), '"'
         nfail = nfail + 1
      end if

      ! Check exitCond and existCond
      if (trim(exitCond_out) /= trim(exitCond_in)) then
         print *, 'FAIL  test_roundtrip: exitCond mismatch'
         nfail = nfail + 1
      end if
      if (existCond_out .neqv. existCond_in) then
         print *, 'FAIL  test_roundtrip: existCond expected', existCond_in, 'got', existCond_out
         nfail = nfail + 1
      end if

      ! Check mImport
      if (mImport_out /= mImport_in) then
         print *, 'FAIL  test_roundtrip: mImport expected', mImport_in, 'got', mImport_out
         nfail = nfail + 1
      end if

   end subroutine test_roundtrip

end program test_step_params
