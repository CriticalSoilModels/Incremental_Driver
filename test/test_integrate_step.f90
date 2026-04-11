program test_integrate_step
   !! Tests for integrate_step — proves the core computation loop works
   !! with no file I/O. Constructs config + state entirely in memory,
   !! wires up the elastic umat_runner_t, and checks results.
   use stdlib_kinds, only: dp
   use indr_step_params,  only: step_config_t, material_state_t
   use indr_model_runner, only: model_runner_t
   use indr_umat_runner,  only: umat_runner_t
   use indr_step_runner,  only: integrate_step
   use indr_umat,         only: UMAT
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-8_dp

   call test_linear_strain_controlled()
   call test_result_count()
   call test_stress_controlled()

   if (nfail == 0) then
      print *, 'PASS  test_integrate_step'
   else
      print *, 'FAIL  test_integrate_step — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   ! ---------------------------------------------------------------------------
   ! Helpers
   ! ---------------------------------------------------------------------------
   subroutine make_runner(E, nu, runner)
      real(dp), intent(in) :: E, nu
      type(umat_runner_t), intent(out) :: runner
      runner%proc   => UMAT
      allocate(runner%props(2), source=[E, nu])
      runner%nstatv  = 1
      runner%cmname  = 'ELASTIC'
   end subroutine make_runner

   subroutine make_state(state)
      type(material_state_t), intent(out) :: state
      real(dp), parameter :: identity33(3,3) = reshape( &
         [1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      state%sig    = 0.0_dp
      state%eps    = 0.0_dp
      state%time   = 0.0_dp
      state%dt     = 0.0_dp
      state%temp   = 0.0_dp
      state%F_start = identity33
      state%F_end   = identity33
      allocate(state%statev(1), source=0.0_dp)
   end subroutine make_state

   subroutine make_config(config, n_inc, deps_total)
      type(step_config_t), intent(out) :: config
      integer,  intent(in) :: n_inc
      real(dp), intent(in) :: deps_total(6)
      config%load_type       = '*LinearLoad'
      config%coord_sys       = '*Cartesian'
      config%n_inc           = n_inc
      config%max_iter        = 10
      config%ifstress        = 0          ! all strain-controlled
      config%delta_load      = 0.0_dp
      config%delta_load(1:6) = deps_total
      config%delta_load_circ = 0.0_dp
      config%phase0          = 0.0_dp
      config%delta_time      = 1.0_dp
      config%delta_temp      = 0.0_dp
      config%has_exit_cond   = .false.
      config%exit_cond       = ''
      config%import_file     = ''
      config%n_import        = 0
      config%columns_in_file = 0
      config%import_factor   = 1.0_dp
      config%dfgrd0          = reshape([1,0,0,0,1,0,0,0,1],[3,3])
      config%dfgrd1          = config%dfgrd0
      config%cMt             = 0.0_dp
      config%cMe             = 0.0_dp
      config%mbinc           = 0.0_dp
   end subroutine make_config

   ! ---------------------------------------------------------------------------

   subroutine test_linear_strain_controlled()
      !! Apply 10 equal strain increments of deps(1)=0.001 to a zero-stress
      !! elastic body. Final sig(1) should equal (lam+2mu)*0.01 exactly.
      real(dp), parameter :: E = 10000.0_dp, nu = 0.3_dp
      integer,  parameter :: n_inc = 10

      type(umat_runner_t)  :: runner
      type(material_state_t)  :: state
      type(step_config_t)     :: config
      type(material_state_t), allocatable :: results(:)

      real(dp) :: deps_total(6), lam, mu, expected_sig1

      lam = nu * E / ((1.0_dp + nu) * (1.0_dp - 2.0_dp*nu))
      mu  = E / (2.0_dp * (1.0_dp + nu))
      expected_sig1 = (lam + 2.0_dp*mu) * 0.01_dp

      deps_total    = 0.0_dp
      deps_total(1) = 0.01_dp

      call make_runner(E, nu, runner)
      call make_state(state)
      call make_config(config, n_inc, deps_total)

      call integrate_step(config, state, runner, results)

      if (abs(state%sig(1) - expected_sig1) > tol * abs(expected_sig1)) then
         print *, 'FAIL  test_linear_strain_controlled: sig(1) =', state%sig(1), &
            ' expected', expected_sig1
         nfail = nfail + 1
      end if

      ! Lateral stress from Poisson effect
      if (abs(state%sig(2) - lam * 0.01_dp) > tol * abs(lam * 0.01_dp)) then
         print *, 'FAIL  test_linear_strain_controlled: sig(2) =', state%sig(2), &
            ' expected', lam * 0.01_dp
         nfail = nfail + 1
      end if

      ! Total strain should equal the applied strain
      if (abs(state%eps(1) - 0.01_dp) > tol) then
         print *, 'FAIL  test_linear_strain_controlled: eps(1) =', state%eps(1)
         nfail = nfail + 1
      end if
   end subroutine test_linear_strain_controlled

   ! ---------------------------------------------------------------------------

   subroutine test_result_count()
      !! results(:) should contain exactly n_inc snapshots when no exit condition fires.
      real(dp), parameter :: E = 10000.0_dp, nu = 0.3_dp
      integer,  parameter :: n_inc = 10

      type(umat_runner_t)  :: runner
      type(material_state_t)  :: state
      type(step_config_t)     :: config
      type(material_state_t), allocatable :: results(:)
      real(dp) :: deps_total(6)

      deps_total    = 0.0_dp
      deps_total(1) = 0.001_dp

      call make_runner(E, nu, runner)
      call make_state(state)
      call make_config(config, n_inc, deps_total)

      call integrate_step(config, state, runner, results)

      if (size(results) /= n_inc) then
         print *, 'FAIL  test_result_count: size(results) =', size(results), 'expected', n_inc
         nfail = nfail + 1
      end if

      ! Each snapshot should have a positive dt
      if (allocated(results)) then
         if (abs(results(1)%dt - 1.0_dp/n_inc) > tol) then
            print *, 'FAIL  test_result_count: results(1)%dt =', results(1)%dt
            nfail = nfail + 1
         end if
      end if
   end subroutine test_result_count

   ! ---------------------------------------------------------------------------

   subroutine test_stress_controlled()
      !! Apply a stress-controlled step: prescribe dsig(1) = 100 kPa.
      !! For linear elastic, the resulting strain should be sig/E (plane strain approx).
      real(dp), parameter :: E = 10000.0_dp, nu = 0.3_dp
      real(dp), parameter :: target_sig1 = 100.0_dp
      integer,  parameter :: n_inc = 10

      type(umat_runner_t)  :: runner
      type(material_state_t)  :: state
      type(step_config_t)     :: config
      type(material_state_t), allocatable :: results(:)
      real(dp) :: deps_total(6)

      call make_runner(E, nu, runner)
      call make_state(state)
      call make_config(config, n_inc, deps_total)

      ! Override: stress-controlled in direction 1 only
      config%ifstress      = 0
      config%ifstress(1)   = 1
      config%delta_load    = 0.0_dp
      config%delta_load(1) = target_sig1

      call integrate_step(config, state, runner, results)

      ! Final stress should match the prescribed target
      if (abs(state%sig(1) - target_sig1) > 0.01_dp * target_sig1) then
         print *, 'FAIL  test_stress_controlled: sig(1) =', state%sig(1), &
            ' expected', target_sig1
         nfail = nfail + 1
      end if
   end subroutine test_stress_controlled

end program test_integrate_step
