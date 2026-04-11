module indr_run_model

   use stdlib_kinds, only: dp
   use stdlib_io, only: open
   use indr_parser, only: splitaLine
   use indr_types,       only: StressAlignment
   use indr_step_params, only: step_config_t, material_state_t, umat_interface
   use indr_umat_runner, only: umat_runner_t
   use indr_step_runner, only: integrate_step
   use indr_matrices, only: MRoscI, MRoscImT, MRendul, MRendulmT, MRosc, MRoscmT, MCart, MCartmT

   use indr_command_line, only: set_inputs
   use indr_file_io, only: read_parameter_file, read_init_conditions_file, set_output_name_from_test_file, &
      write_line_output_data, write_output_file_header
   use indr_alignment, only: readAlignment, tryAlignStress
   use indr_loads
   use indr_constants, only: iter_lower_limit

   implicit none(type, external)

contains
   subroutine run_model(UMAT)
      procedure(umat_interface) :: UMAT

      real(dp), parameter :: delta(3,3) = reshape([1,0,0,0,1,0,0,0,1],[3,3])

      integer :: nstatv, nprops
      integer :: test_file_id, output_file_id
      integer :: iostat
      integer :: every, ievery
      integer :: i

      ! Init-conditions temporaries (used to populate state before the step loop)
      real(dp) :: stress(6), stran(6), time(2), dtime, temp
      real(dp), allocatable :: statev(:), r_statev(:)
      character(len=15), allocatable :: statevHead(:)
      real(dp), allocatable :: props(:)
      character(len=80) :: cmname

      character(len=40) :: keywords(10), outputfilename, &
         parametersfilename, initialconditionsfilename, testfilename
      character(len=260) :: inputline(6), aLine, heading

      logical :: verbose, okSplit
      real(dp), dimension(6) :: mb   ! ObeyRestrictions RHS (read by read_obey_restrictions_load)

      integer :: ikeyword, iRepetition, nRepetitions, kStep, iStep, nSteps

      type(StressAlignment) :: align
      type(step_config_t) :: ofStep(30)   ! saved step configs for *Repetition
      type(step_config_t) :: config       ! current step configuration

      type(umat_runner_t) :: umat_runner
      type(material_state_t) :: state
      type(material_state_t), allocatable :: results(:)

      ! -----------------------------------------------------------------------
      ! [1]  File names and model parameters
      ! -----------------------------------------------------------------------
      call set_inputs(parametersfilename, initialconditionsfilename, &
         testfilename, outputfilename, verbose)

      call read_parameter_file(parametersfilename, nprops, props, cmname)

      call read_init_conditions_file(initialconditionsfilename, stress, time, stran, &
         dtime, temp, statev, r_statev, statevHead, nstatv)

      ! -----------------------------------------------------------------------
      ! [2]  Set up the UMAT runner and initial material state
      ! -----------------------------------------------------------------------
      umat_runner%proc   => UMAT
      umat_runner%props  = props
      umat_runner%cmname = cmname
      umat_runner%nstatv = nstatv

      state%sig    = stress
      state%eps    = stran
      state%statev = statev
      state%time   = time
      state%dt     = dtime
      state%temp   = temp
      state%F_start = delta
      state%F_end   = delta

      ! -----------------------------------------------------------------------
      ! [3]  Open test and output files; write header and initial state
      ! -----------------------------------------------------------------------
      test_file_id = open(testfilename)

      call set_output_name_from_test_file(test_file_id, outputfilename, heading)

      output_file_id = open(outputfilename, "w")

      call write_output_file_header(output_file_id, heading, nstatv)
      call write_line_output_data(output_file_id, state%time, state%dt, state%eps, state%sig, state%statev)

      ! -----------------------------------------------------------------------
      ! [4]  Step loop
      ! -----------------------------------------------------------------------
      kStep = 0
      do_keyword: do ikeyword = 1, 10000

         read(test_file_id, '(a)') keywords(1)

         if (iostat < 0) error stop 'Reached end of test file unexpectedly'

         keywords(1) = trim(keywords(1))
         if (keywords(1) == '*Repetition') then
            read(test_file_id, *) nSteps, nRepetitions
         else
            nRepetitions = 1
            nSteps       = 1
            keywords(2)  = keywords(1)
         end if

         do_repet: do iRepetition = 1, nRepetitions
            do_istep: do iStep = 1, nSteps
               kStep = kStep + 1

               write(*,'(A,I3,A,I3,A,F9.4,A,F9.4)') &
                  ' ikeyword = ', ikeyword, &
                  ' kstep = ', kStep, &
                  ' TEMP = ', state%temp, &
                  ' TIME = ', state%time(1)

               if (iRepetition > 1) then
                  config      = ofStep(iStep)
                  keywords(2) = config%load_type
               end if

               if (keywords(1) == '*Repetition') read(test_file_id, '(a)') keywords(2)

               call splitaLine(keywords(2), '?', keywords(2), config%exit_cond, config%has_exit_cond)
               keywords(2) = trim(keywords(2))

               ! Defaults — read routines only set the fields they use.
               config%load_type       = keywords(2)
               config%ifstress        = 0
               config%delta_load      = 0.0_dp
               config%delta_load_circ = 0.0_dp
               config%phase0          = 0.0_dp
               config%dfgrd0          = delta
               config%dfgrd1          = delta
               config%cMt             = 0.0_dp
               config%cMe             = 0.0_dp
               config%mbinc           = 0.0_dp
               config%columns_in_file = 0
               config%import_factor   = 1.0_dp
               config%n_import        = 0

               if (keywords(2) == '*DeformationGradient') then
                  call read_deformation_gradient_load(test_file_id, config, every)

               else if (keywords(2) == '*CirculatingLoad') then
                  call read_circulating_load(test_file_id, config, every)

               else if (keywords(2) == '*LinearLoad') then
                  call read_linear_load(test_file_id, config, every)

               else if (keywords(2)(1:11) == '*ImportFile') then
                  call read_file_load(test_file_id, config, every, align)

               else if (keywords(2) == '*OedometricE1') then
                  call read_oedometric_load(test_file_id, config, every)

               else if (keywords(2) == '*OedometricS1') then
                  call read_oedometric_S1_load(test_file_id, config, every)

               else if (keywords(2) == '*TriaxialE1') then
                  call read_triaxial_e1_load(test_file_id, config, every)

               else if (keywords(2) == '*TriaxialS1') then
                  call read_triaxial_s1_load(test_file_id, config, every)

               else if (keywords(2) == '*TriaxialUEq') then
                  call read_triaxial_ueq_load(test_file_id, config, every)

               else if (keywords(2) == '*TriaxialUq') then
                  call read_triaxial_uq_load(test_file_id, config, every)

               else if (keywords(2) == '*PureRelaxation') then
                  call read_pure_relaxation_load(test_file_id, config, every)

               else if (keywords(2) == '*PureCreep') then
                  call read_pure_creep_load(test_file_id, config, every)

               else if (keywords(2) == '*UndrainedCreep') then
                  call read_undrained_creep(test_file_id, config, every)

               else if (keywords(2) == '*ObeyRestrictions') then
                  call read_obey_restrictions_load(test_file_id, config, every, mb)

               else if (keywords(2) == '*PerturbationsS') then
                  call read_perturbations_S_load(test_file_id, config, every)

               else if (keywords(2) == '*PerturbationsE') then
                  call read_perturbations_E_load(test_file_id, config, every)

               else if (keywords(2) == '*RandomWalk') then
                  call read_random_walk_load(test_file_id, config, every)

               else if (keywords(2) == '*End') then
                  print *, '*End encountered in test.inp'
                  close(test_file_id)
                  close(output_file_id)
                  return
               else
                  write(*,*) 'error: unknown keywords(2) =', keywords(2)
                  error stop 'stopped by unknown keyword in test.inp'
               end if

               if (keywords(1) == '*Repetition' .and. iRepetition == 1) then
                  ofStep(iStep) = config
               end if

               ! -----------------------------------------------------------------
               ! Run the step — no file I/O inside integrate_step (except *ImportFile reads)
               ! -----------------------------------------------------------------
               call integrate_step(config, state, umat_runner, results, align)

               ! -----------------------------------------------------------------
               ! Write results respecting 'every' output frequency
               ! -----------------------------------------------------------------
               ievery = 1
               do i = 1, size(results)
                  if (ievery == 1) then
                     ! results(i)%time already holds the end-of-increment time;
                     ! pass dt=0 so write_line_output_data does not add it again.
                     call write_line_output_data(output_file_id, results(i)%time, 0.0_dp, &
                        results(i)%eps, results(i)%sig, results(i)%statev)
                  end if
                  ievery = ievery + 1
                  if (ievery > every) ievery = 1
               end do

            end do do_istep
         end do do_repet
      end do do_keyword

      close(test_file_id)
      close(output_file_id)

   end subroutine run_model
end module indr_run_model
