!! Top-level driver: parses the test file, runs each step through
!! integrate_step, and writes results. No per-increment logic lives here.
module indr_run_model

   use stdlib_kinds, only: dp
   use stdlib_io,    only: open

   use indr_step_params,  only: material_state_t, umat_interface
   use indr_umat_runner,  only: umat_runner_t
   use indr_step_runner,  only: integrate_step
   use indr_test_parser,  only: step_record_t, parse_test_file
   use indr_types,        only: StressAlignment

   use indr_command_line, only: set_inputs
   use indr_file_io,      only: read_parameter_file, read_init_conditions_file, &
      write_line_output_data, write_output_file_header, write_step_output

   implicit none(type, external)

contains

   subroutine run_model(UMAT)
      !! Orchestrate a full driver run:
      !!   1. Read parameters and initial conditions.
      !!   2. Parse test.inp upfront into steps(:) via parse_test_file.
      !!   3. For each step: integrate_step → write_step_output.
      procedure(umat_interface) :: UMAT

      real(dp), parameter :: delta(3,3) = reshape([1,0,0,0,1,0,0,0,1],[3,3])

      integer :: nstatv, nprops, kstep, isep
      integer :: test_file_id, output_file_id

      ! Initial-conditions temporaries
      real(dp)              :: stress(6), stran(6), time(2), dtime, temp
      real(dp), allocatable :: statev(:), r_statev(:), props(:)
      character(len=15), allocatable :: statevHead(:)
      character(len=80) :: cmname

      ! File names
      character(len=40) :: outputfilename, parametersfilename, &
                           initialconditionsfilename, testfilename
      logical :: verbose

      ! Parsed test-file data
      character(len=260)              :: heading, heading_comment
      type(step_record_t), allocatable :: steps(:)
      type(StressAlignment)           :: align
      integer                         :: n_steps

      ! Runner and persistent state
      type(umat_runner_t)              :: umat_runner
      type(material_state_t)           :: state
      type(material_state_t), allocatable :: results(:)

      ! -----------------------------------------------------------------------
      ! [1]  Parameters and initial conditions
      ! -----------------------------------------------------------------------
      call set_inputs(parametersfilename, initialconditionsfilename, &
         testfilename, outputfilename, verbose)

      call read_parameter_file(parametersfilename, nprops, props, cmname)

      call read_init_conditions_file(initialconditionsfilename, stress, time, stran, &
         dtime, temp, statev, r_statev, statevHead, nstatv)

      ! -----------------------------------------------------------------------
      ! [2]  Set up runner and initial material state
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
      ! [3]  Parse test.inp upfront — returns flat steps(:) with *Repetition expanded
      ! -----------------------------------------------------------------------
      test_file_id = open(testfilename)
      call parse_test_file(test_file_id, heading, steps, n_steps, align)
      close(test_file_id)

      ! Extract output filename and optional heading comment from the first line.
      ! Format: "output.txt # optional comment"  or just "output.txt"
      isep = index(heading, '#')
      if (isep == 0) then
         outputfilename  = trim(adjustl(heading))
         heading_comment = '#'
      else
         outputfilename  = trim(adjustl(heading(:isep-1)))
         heading_comment = trim(adjustl(heading(isep+1:)))
      end if

      ! -----------------------------------------------------------------------
      ! [4]  Open output file; write column headers and initial state
      ! -----------------------------------------------------------------------
      output_file_id = open(outputfilename, "w")
      call write_output_file_header(output_file_id, heading_comment, nstatv)
      call write_line_output_data(output_file_id, state%time, state%dt, &
         state%eps, state%sig, state%statev)

      ! -----------------------------------------------------------------------
      ! [5]  Step loop
      ! -----------------------------------------------------------------------
      do kstep = 1, n_steps
         write(*,'(A,I4,A,F9.4,A,F9.4)') &
            ' kstep = ', kstep, &
            '  TEMP = ', state%temp, &
            '  TIME = ', state%time(1)

         call integrate_step(steps(kstep)%config, state, umat_runner, results, align)
         call write_step_output(output_file_id, results, steps(kstep)%write_freq)
      end do

      close(output_file_id)

   end subroutine run_model

end module indr_run_model
