program test_read_loads
   !! Unit tests for the read_* subroutines in indr_loads.
   !! Each test writes input text to a scratch file, calls the routine,
   !! then checks the output values against known expectations.
   use stdlib_kinds, only: dp
   use indr_loads, only: read_linear_load, read_circulating_load, &
                         read_deformation_gradient_load,           &
                         read_oedometric_load, read_pure_creep_load, &
                         read_undrained_creep
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-12_dp

   call test_linear_load_stress_ctrl()
   call test_linear_load_mixed_ctrl()
   call test_circulating_load()
   call test_deformation_gradient_load()
   call test_oedometric_load()
   call test_pure_creep_load()
   call test_undrained_creep()

   if (nfail == 0) then
      print *, 'PASS  test_read_loads'
   else
      print *, 'FAIL  test_read_loads — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   ! ---------------------------------------------------------------------------
   ! read_linear_load
   ! File format:
   !   ninc maxiter deltaTime [T deltaTemp] [: every]
   !   keyword                    (coord sys)
   !   ifstress(i) deltaLoad(i)   (x6)
   ! ---------------------------------------------------------------------------

   subroutine test_linear_load_stress_ctrl()
      integer  :: fid, ninc, maxiter, write_freq
      integer  :: ifstress(6)
      real(dp) :: deltaTime, deltaTemp, deltaLoad(9)
      character(10) :: keyword

      open(newunit=fid, status='scratch')
      write(fid, '(a)') '10 50 1.0 T 0.0'
      write(fid, '(a)') '*Roscoe'
      write(fid, '(a)') '1  100.0'
      write(fid, '(a)') '1  200.0'
      write(fid, '(a)') '1  300.0'
      write(fid, '(a)') '0  0.0'
      write(fid, '(a)') '0  0.0'
      write(fid, '(a)') '0  0.0'
      rewind(fid)

      deltaLoad = 0.0_dp
      call read_linear_load(fid, ninc, maxiter, deltaTime, deltaTemp, write_freq, keyword, ifstress, deltaLoad)
      close(fid)

      call check_int(ninc,    10,     'linear_stress_ctrl: ninc',    nfail)
      call check_int(maxiter, 50,     'linear_stress_ctrl: maxiter', nfail)
      call check_real(deltaTime, 1.0_dp, tol, 'linear_stress_ctrl: deltaTime', nfail)
      call check_real(deltaTemp, 0.0_dp, tol, 'linear_stress_ctrl: deltaTemp', nfail)
      call check_int(write_freq, 1,   'linear_stress_ctrl: write_freq', nfail)
      if (trim(keyword) /= '*Roscoe') then
         print *, 'FAIL  linear_stress_ctrl: keyword expected "*Roscoe" got "', trim(keyword), '"'
         nfail = nfail + 1
      end if
      call check_int(ifstress(1), 1,      'linear_stress_ctrl: ifstress(1)', nfail)
      call check_int(ifstress(4), 0,      'linear_stress_ctrl: ifstress(4)', nfail)
      call check_real(deltaLoad(1), 100.0_dp, tol, 'linear_stress_ctrl: deltaLoad(1)', nfail)
      call check_real(deltaLoad(2), 200.0_dp, tol, 'linear_stress_ctrl: deltaLoad(2)', nfail)
      call check_real(deltaLoad(3), 300.0_dp, tol, 'linear_stress_ctrl: deltaLoad(3)', nfail)
      call check_real(deltaLoad(4),   0.0_dp, tol, 'linear_stress_ctrl: deltaLoad(4)', nfail)
   end subroutine test_linear_load_stress_ctrl

   subroutine test_linear_load_mixed_ctrl()
      !! Verify mixed stress/strain control and write_freq from ':' field.
      integer  :: fid, ninc, maxiter, write_freq
      integer  :: ifstress(6)
      real(dp) :: deltaTime, deltaTemp, deltaLoad(9)
      character(10) :: keyword

      open(newunit=fid, status='scratch')
      write(fid, '(a)') '20 100 2.5 T 10.0 : 4'
      write(fid, '(a)') '*Cartesian'
      write(fid, '(a)') '1  50.0'
      write(fid, '(a)') '0  0.01'
      write(fid, '(a)') '0  0.0'
      write(fid, '(a)') '0  0.0'
      write(fid, '(a)') '0  0.0'
      write(fid, '(a)') '0  0.0'
      rewind(fid)

      deltaLoad = 0.0_dp
      call read_linear_load(fid, ninc, maxiter, deltaTime, deltaTemp, write_freq, keyword, ifstress, deltaLoad)
      close(fid)

      call check_int(ninc,       20,     'linear_mixed: ninc',      nfail)
      call check_int(maxiter,   100,     'linear_mixed: maxiter',   nfail)
      call check_real(deltaTime, 2.5_dp, tol, 'linear_mixed: deltaTime', nfail)
      call check_real(deltaTemp, 10.0_dp, tol, 'linear_mixed: deltaTemp', nfail)
      call check_int(write_freq,  4,     'linear_mixed: write_freq', nfail)
      call check_int(ifstress(1), 1,     'linear_mixed: ifstress(1)', nfail)
      call check_int(ifstress(2), 0,     'linear_mixed: ifstress(2)', nfail)
      call check_real(deltaLoad(1), 50.0_dp, tol, 'linear_mixed: deltaLoad(1)', nfail)
      call check_real(deltaLoad(2), 0.01_dp, tol, 'linear_mixed: deltaLoad(2)', nfail)
   end subroutine test_linear_load_mixed_ctrl

   ! ---------------------------------------------------------------------------
   ! read_circulating_load
   ! File format:
   !   ninc maxiter deltaTime [T deltaTemp] [: every]
   !   keyword
   !   ifstress(i) amp(i) phase(i) bias(i)   (x6)
   ! ---------------------------------------------------------------------------

   subroutine test_circulating_load()
      integer  :: fid, ninc, maxiter, write_freq
      integer  :: ifstress(6)
      real(dp) :: deltaTime, deltaTemp, deltaLoad(6), deltaLoadCirc(6), phase(6)
      character(10) :: keyword

      open(newunit=fid, status='scratch')
      write(fid, '(a)') '100 200 5.0 T 0.0'
      write(fid, '(a)') '*Roscoe'
      write(fid, '(a)') '1  10.0  0.0  5.0'   ! stress-ctrl: amp=10, phase=0, bias=5
      write(fid, '(a)') '0   0.0  0.0  0.0'
      write(fid, '(a)') '0   0.0  0.0  0.0'
      write(fid, '(a)') '0   0.0  0.0  0.0'
      write(fid, '(a)') '0   0.0  0.0  0.0'
      write(fid, '(a)') '0   0.0  0.0  0.0'
      rewind(fid)

      call read_circulating_load(fid, ninc, maxiter, deltaTime, deltaTemp, write_freq, &
                                 keyword, deltaLoad, ifstress, deltaLoadCirc, phase)
      close(fid)

      call check_int(ninc,     100,    'circulating: ninc',    nfail)
      call check_int(maxiter,  200,    'circulating: maxiter', nfail)
      call check_real(deltaTime, 5.0_dp, tol, 'circulating: deltaTime', nfail)
      if (trim(keyword) /= '*Roscoe') then
         print *, 'FAIL  circulating: keyword expected "*Roscoe" got "', trim(keyword), '"'
         nfail = nfail + 1
      end if
      call check_int(ifstress(1), 1,       'circulating: ifstress(1)', nfail)
      call check_int(ifstress(2), 0,       'circulating: ifstress(2)', nfail)
      call check_real(deltaLoadCirc(1), 10.0_dp, tol, 'circulating: amp(1)',   nfail)
      call check_real(phase(1),          0.0_dp, tol, 'circulating: phase(1)', nfail)
      call check_real(deltaLoad(1),       5.0_dp, tol, 'circulating: bias(1)',  nfail)
      call check_real(deltaLoad(2),       0.0_dp, tol, 'circulating: bias(2)',  nfail)
   end subroutine test_circulating_load

   ! ---------------------------------------------------------------------------
   ! read_deformation_gradient_load
   ! File format:
   !   ninc maxiter deltaTime [T deltaTemp] [: every]
   !   F11
   !   F22
   !   ... (9 values total in stretch order)
   ! ---------------------------------------------------------------------------

   subroutine test_deformation_gradient_load()
      integer  :: fid, ninc, maxiter, write_freq
      real(dp) :: deltaTime, deltaTemp, deltaLoad(9)
      character(40) :: keyword

      open(newunit=fid, status='scratch')
      write(fid, '(a)') '5 10 0.5 T 0.0'
      write(fid, '(a)') '1.1'   ! F11
      write(fid, '(a)') '1.0'   ! F22
      write(fid, '(a)') '1.0'   ! F33
      write(fid, '(a)') '0.0'   ! F12
      write(fid, '(a)') '0.0'   ! F21
      write(fid, '(a)') '0.0'   ! F13
      write(fid, '(a)') '0.0'   ! F31
      write(fid, '(a)') '0.0'   ! F23
      write(fid, '(a)') '0.0'   ! F32
      rewind(fid)

      call read_deformation_gradient_load(fid, ninc, maxiter, deltaTime, deltaTemp, write_freq, keyword, deltaLoad)
      close(fid)

      call check_int(ninc,    5,      'defgrad: ninc',    nfail)
      call check_int(maxiter, 10,     'defgrad: maxiter', nfail)
      call check_real(deltaTime, 0.5_dp, tol, 'defgrad: deltaTime', nfail)
      if (trim(keyword) /= '*Cartesian') then
         print *, 'FAIL  defgrad: keyword expected "*Cartesian" got "', trim(keyword), '"'
         nfail = nfail + 1
      end if
      call check_real(deltaLoad(1), 1.1_dp, tol, 'defgrad: deltaLoad(1)', nfail)
      call check_real(deltaLoad(2), 1.0_dp, tol, 'defgrad: deltaLoad(2)', nfail)
      call check_real(deltaLoad(4), 0.0_dp, tol, 'defgrad: deltaLoad(4)', nfail)
   end subroutine test_deformation_gradient_load

   ! ---------------------------------------------------------------------------
   ! read_oedometric_load — thin wrapper, verifies hardcoded keywords and load
   ! ---------------------------------------------------------------------------

   subroutine test_oedometric_load()
      integer  :: fid, ninc, maxiter, write_freq
      real(dp) :: deltaTime, deltaTemp, deltaLoad1
      character(40) :: kw2, kw3

      open(newunit=fid, status='scratch')
      write(fid, '(a)') '15 30 3.0 T 0.0'
      write(fid, '(a)') '0.05'
      rewind(fid)

      call read_oedometric_load(fid, ninc, maxiter, deltaTime, deltaTemp, write_freq, kw2, kw3, deltaLoad1)
      close(fid)

      call check_int(ninc,    15,      'oedometric: ninc',    nfail)
      call check_int(maxiter, 30,      'oedometric: maxiter', nfail)
      call check_real(deltaTime,  3.0_dp, tol, 'oedometric: deltaTime',  nfail)
      call check_real(deltaLoad1, 0.05_dp, tol, 'oedometric: deltaLoad1', nfail)
      if (trim(kw2) /= '*LinearLoad') then
         print *, 'FAIL  oedometric: kw2 expected "*LinearLoad" got "', trim(kw2), '"'
         nfail = nfail + 1
      end if
      if (trim(kw3) /= '*Cartesian') then
         print *, 'FAIL  oedometric: kw3 expected "*Cartesian" got "', trim(kw3), '"'
         nfail = nfail + 1
      end if
   end subroutine test_oedometric_load

   ! ---------------------------------------------------------------------------
   ! read_pure_creep_load — hardcodes ifstress=1
   ! ---------------------------------------------------------------------------

   subroutine test_pure_creep_load()
      integer  :: fid, ninc, maxiter, write_freq, ifstress(6)
      real(dp) :: deltaTime, deltaTemp
      character(40) :: kw2, kw3

      open(newunit=fid, status='scratch')
      write(fid, '(a)') '10 20 1.0 T 0.0'
      rewind(fid)

      call read_pure_creep_load(fid, ninc, maxiter, deltaTime, deltaTemp, write_freq, kw2, kw3, ifstress)
      close(fid)

      call check_int(ninc,    10,  'pure_creep: ninc',    nfail)
      if (any(ifstress /= 1)) then
         print *, 'FAIL  pure_creep: ifstress should be all 1, got', ifstress
         nfail = nfail + 1
      end if
      if (trim(kw2) /= '*LinearLoad') then
         print *, 'FAIL  pure_creep: kw2 expected "*LinearLoad" got "', trim(kw2), '"'
         nfail = nfail + 1
      end if
   end subroutine test_pure_creep_load

   ! ---------------------------------------------------------------------------
   ! read_undrained_creep — specific ifstress pattern: (0,1,1,1,1,1)
   ! ---------------------------------------------------------------------------

   subroutine test_undrained_creep()
      integer  :: fid, ninc, maxiter, write_freq, ifstress(6)
      real(dp) :: deltaTime, deltaTemp
      character(40) :: kw2, kw3

      open(newunit=fid, status='scratch')
      write(fid, '(a)') '10 20 1.0 T 0.0'
      rewind(fid)

      call read_undrained_creep(fid, ninc, maxiter, deltaTime, deltaTemp, write_freq, kw2, kw3, ifstress)
      close(fid)

      call check_int(ifstress(1), 0,  'undrained_creep: ifstress(1) should be 0', nfail)
      if (any(ifstress(2:6) /= 1)) then
         print *, 'FAIL  undrained_creep: ifstress(2:6) should be all 1, got', ifstress(2:6)
         nfail = nfail + 1
      end if
      if (trim(kw3) /= '*Roscoe') then
         print *, 'FAIL  undrained_creep: kw3 expected "*Roscoe" got "', trim(kw3), '"'
         nfail = nfail + 1
      end if
   end subroutine test_undrained_creep

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

end program test_read_loads
