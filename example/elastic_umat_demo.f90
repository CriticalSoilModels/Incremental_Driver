!! Example: run a linear strain-controlled load step using integrate_step directly.
!!
!! Demonstrates the file-free API using the incremental_driver facade module:
!!   - No test.inp, no parameters file, no initial-conditions file.
!!   - State is constructed in memory.
!!   - The user-provided UMAT is defined in elastic.f90 in this directory.
!!   - Results are returned as an array of material_state_t snapshots.
!!
!! To build and run:
!!   fpm run --example elastic_umat_demo

program elastic_umat_demo
   use incremental_driver, only: step_config_t, material_state_t, &
                                 umat_runner_t, integrate_step, STRAIN_CTRL
   use indr_umat, only: UMAT
   use stdlib_kinds, only: dp
   implicit none(type, external)

   real(dp), parameter :: E           = 10000.0_dp
   real(dp), parameter :: nu          = 0.3_dp
   real(dp), parameter :: total_eps1  = 0.01_dp   !! total axial strain applied over the step [-]
   real(dp), parameter :: I3(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])

   type(umat_runner_t)              :: runner
   type(material_state_t)           :: state
   type(step_config_t)              :: config
   type(material_state_t), allocatable :: results(:)

   real(dp) :: lam, mu, expected_sig1
   integer  :: i

   ! --- Wire up the runner ---
   runner%proc   => UMAT
   runner%cmname =  'ELASTIC'
   runner%nstatv =  1
   allocate(runner%props(2), source=[E, nu])

   ! --- Initial state: stress-free, zero strain ---
   state%sig     = 0.0_dp
   state%eps     = 0.0_dp
   state%time    = 0.0_dp
   state%dt      = 0.0_dp
   state%temp    = 0.0_dp
   state%F_start = I3
   state%F_end   = I3
   allocate(state%statev(1), source=0.0_dp)

   ! --- Step: 10 strain-controlled increments, deps(1) = 0.001 each ---
   config%load_type       = '*LinearLoad'
   config%coord_sys       = '*Cartesian'
   config%n_inc           = 10
   config%max_iter        = 1
   config%ifstress        = STRAIN_CTRL
   config%delta_load      = 0.0_dp
   config%delta_load(1)   = total_eps1
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
   config%dfgrd0          = I3
   config%dfgrd1          = I3
   config%cMt             = 0.0_dp
   config%cMe             = 0.0_dp
   config%mbinc           = 0.0_dp

   ! --- Run ---
   call integrate_step(config, state, runner, results)

   ! --- Print results ---
   write(*,'(a6,3a14)') 'inc', 'eps(1)', 'sig(1)', 'sig(2)'
   do i = 1, size(results)
      write(*,'(i6,3f14.4)') i, results(i)%eps(1), results(i)%sig(1), results(i)%sig(2)
   end do

   ! --- Analytic check ---
   lam = nu * E / ((1.0_dp + nu) * (1.0_dp - 2.0_dp*nu))
   mu  = E / (2.0_dp * (1.0_dp + nu))
   expected_sig1 = (lam + 2.0_dp*mu) * total_eps1

   write(*,*)
   write(*,'(a,f10.4)') 'Final sig(1)    : ', state%sig(1)
   write(*,'(a,f10.4)') 'Analytic sig(1) : ', expected_sig1
   if (abs(state%sig(1) - expected_sig1) < 1.0e-6_dp * abs(expected_sig1)) then
      write(*,*) 'PASS'
   else
      write(*,*) 'FAIL'
      stop 1
   end if

end program elastic_umat_demo
