module indr_step_parser
   !! Reads an entire test.inp file and returns all step records.
   !!
   !! Phase 6: supports all load types except *Repetition.
   !! Phase 7: adds *Repetition block expansion.
   use stdlib_kinds,      only: dp
   use indr_types,        only: step_config_t, StressAlignment, STRAIN_CTRL
   use indr_loads,        only: read_linear_load, read_circulating_load,         &
                                read_deformation_gradient_load, read_file_load,  &
                                read_oedometric_load, read_oedometric_S1_load,   &
                                read_triaxial_e1_load, read_triaxial_s1_load,    &
                                read_triaxial_ueq_load, read_triaxial_uq_load,   &
                                read_pure_relaxation_load, read_pure_creep_load, &
                                read_undrained_creep, read_obey_restrictions_load, &
                                read_perturbations_S_load, read_perturbations_E_load, &
                                read_random_walk_load
   use indr_parser,       only: splitaLine
   use indr_alignment,    only: readAlignment
   use indr_file_io,      only: load_import_data

   implicit none
   private
   public :: step_record_t, parse_test_file

   integer, parameter :: INIT_CAP = 64  !! Initial step buffer capacity; doubles on overflow

   type step_record_t
      !! Combines a fully-populated step configuration with the output write frequency.
      type(step_config_t) :: config
      integer             :: write_freq   !! output every write_freq-th increment
   end type step_record_t

contains

   subroutine parse_test_file(file_id, heading, steps, n_steps, align)
      !! Read an entire test.inp from the already-open unit file_id.
      !! Returns all steps as an allocated array and sets heading and align.
      integer,                          intent(in)  :: file_id
      character(len=*),                 intent(out) :: heading
      type(step_record_t), allocatable, intent(out) :: steps(:)
      integer,                          intent(out) :: n_steps
      type(StressAlignment),            intent(out) :: align

      type(step_record_t), allocatable :: buf(:), rep_buf(:)
      character(len=260)  :: line
      character(len=40)   :: keyword, exit_cond, kw40
      logical             :: has_exit
      type(step_config_t) :: config
      integer             :: write_freq
      real(dp)            :: mb(6)
      real(dp), parameter :: delta(3,3) = reshape([1,0,0,0,1,0,0,0,1], [3,3])

      ! *Repetition workspace
      integer :: n_rep_steps, n_reps, i_rep, i_step

      ! Initialise outputs
      align%active = .false.
      n_steps = 0
      allocate(buf(INIT_CAP))

      ! Read heading (first line)
      read(file_id, '(a)') line
      heading = trim(line)

      ! Main parsing loop
      do
         read(file_id, '(a)', end=100) line
         line = adjustl(line)
         if (len_trim(line) == 0) cycle          ! skip blank lines
         if (line(1:1) == '!')  cycle             ! skip comment lines

         ! Extract keyword and optional exit condition after '?'
         kw40 = line(1:40)
         call splitaLine(kw40, '?', keyword, exit_cond, has_exit)
         keyword = adjustl(trim(keyword))

         ! Set defaults before dispatching to a reader
         call set_defaults(config, keyword, delta)
         ! Restore exit condition parsed from keyword line
         config%has_exit_cond = has_exit
         config%exit_cond     = exit_cond

         select case (trim(keyword))

         case ('*End')
            exit

         case ('*LinearLoad')
            call read_linear_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*CirculatingLoad')
            call read_circulating_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*DeformationGradient')
            call read_deformation_gradient_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*OedometricE1')
            call read_oedometric_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*OedometricS1')
            call read_oedometric_S1_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*TriaxialE1')
            call read_triaxial_e1_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*TriaxialS1')
            call read_triaxial_s1_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*TriaxialUEq')
            call read_triaxial_ueq_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*TriaxialUq')
            call read_triaxial_uq_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*PureRelaxation')
            call read_pure_relaxation_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*PureCreep')
            call read_pure_creep_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*UndrainedCreep')
            call read_undrained_creep(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*ObeyRestrictions')
            call read_obey_restrictions_load(file_id, config, write_freq, mb)
            call store(buf, n_steps, config, write_freq)

         case ('*PerturbationsS')
            call read_perturbations_S_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*PerturbationsE')
            call read_perturbations_E_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*RandomWalk')
            call read_random_walk_load(file_id, config, write_freq)
            call store(buf, n_steps, config, write_freq)

         case ('*Repetition')
            call read_repetition_block(file_id, align, delta, rep_buf, n_rep_steps, n_reps)
            do i_rep = 1, n_reps
               do i_step = 1, n_rep_steps
                  call store(buf, n_steps, rep_buf(i_step)%config, rep_buf(i_step)%write_freq)
               end do
            end do

         case default
            ! Check for *ImportFile prefix (may have |filename|id appended)
            if (keyword(1:11) == '*ImportFile') then
               call read_file_load(file_id, config, write_freq, align)
               ! read_file_load has align intent(in) so we call readAlignment directly
               call readAlignment(align, config%import_file)
               call load_import_data(config)
               call store(buf, n_steps, config, write_freq)
            end if
            ! Unknown keywords are silently ignored

         end select
      end do

100   continue

      allocate(steps(n_steps))
      steps = buf(1:n_steps)

   end subroutine parse_test_file

   ! ---------------------------------------------------------------------------
   ! Private helpers
   ! ---------------------------------------------------------------------------

   subroutine set_defaults(config, keyword, delta)
      !! Populate config with safe defaults before calling a load reader.
      type(step_config_t), intent(out) :: config
      character(len=40),   intent(in)  :: keyword
      real(dp),            intent(in)  :: delta(3,3)

      config%load_type       = keyword
      config%ifstress        = STRAIN_CTRL
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
      config%has_exit_cond   = .false.
      config%exit_cond       = ''
      config%import_file     = ''
      config%coord_sys       = ''
      config%delta_temp      = 0.0_dp
      config%delta_time      = 0.0_dp
      config%n_inc           = 0
      config%max_iter        = 0
   end subroutine set_defaults

   subroutine store(buf, n, config, write_freq)
      !! Append one step_record_t to the buffer, doubling capacity when full.
      type(step_record_t), allocatable, intent(inout) :: buf(:)
      integer,             intent(inout) :: n
      type(step_config_t), intent(in)    :: config
      integer,             intent(in)    :: write_freq

      type(step_record_t), allocatable :: tmp(:)

      n = n + 1
      if (n > size(buf)) then
         allocate(tmp(size(buf) * 2))
         tmp(1:n-1) = buf(1:n-1)
         call move_alloc(tmp, buf)
      end if
      buf(n)%config     = config
      buf(n)%write_freq = write_freq
   end subroutine store

   subroutine read_repetition_block(file_id, align, delta, rep_buf, n_rep_steps, n_reps)
      !! Read a *Repetition block: first line gives nSteps nRepetitions,
      !! then reads nSteps step blocks, storing them in rep_buf.
      integer,               intent(in)    :: file_id
      type(StressAlignment), intent(inout) :: align
      real(dp),              intent(in)    :: delta(3,3)
      type(step_record_t), allocatable, intent(out) :: rep_buf(:)
      integer,                          intent(out) :: n_rep_steps
      integer,                          intent(out) :: n_reps

      character(len=260)  :: line
      character(len=40)   :: keyword, exit_cond, kw40
      type(step_config_t) :: config
      integer             :: write_freq, i
      real(dp)            :: mb(6)
      logical             :: has_exit

      n_rep_steps = 0
      n_reps      = 1

      ! Read: nSteps nRepetitions
      read(file_id, *) n_rep_steps, n_reps
      allocate(rep_buf(n_rep_steps))

      ! Read exactly n_rep_steps step blocks
      do i = 1, n_rep_steps
         ! Read the keyword line for this step
         do
            read(file_id, '(a)') line
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '!') cycle
            exit
         end do

         kw40 = line(1:40)
         call splitaLine(kw40, '?', keyword, exit_cond, has_exit)
         keyword = adjustl(trim(keyword))

         call set_defaults(config, keyword, delta)
         config%has_exit_cond = has_exit
         config%exit_cond     = exit_cond

         select case (trim(keyword))
         case ('*LinearLoad')
            call read_linear_load(file_id, config, write_freq)
         case ('*CirculatingLoad')
            call read_circulating_load(file_id, config, write_freq)
         case ('*DeformationGradient')
            call read_deformation_gradient_load(file_id, config, write_freq)
         case ('*OedometricE1')
            call read_oedometric_load(file_id, config, write_freq)
         case ('*OedometricS1')
            call read_oedometric_S1_load(file_id, config, write_freq)
         case ('*TriaxialE1')
            call read_triaxial_e1_load(file_id, config, write_freq)
         case ('*TriaxialS1')
            call read_triaxial_s1_load(file_id, config, write_freq)
         case ('*TriaxialUEq')
            call read_triaxial_ueq_load(file_id, config, write_freq)
         case ('*TriaxialUq')
            call read_triaxial_uq_load(file_id, config, write_freq)
         case ('*PureRelaxation')
            call read_pure_relaxation_load(file_id, config, write_freq)
         case ('*PureCreep')
            call read_pure_creep_load(file_id, config, write_freq)
         case ('*UndrainedCreep')
            call read_undrained_creep(file_id, config, write_freq)
         case ('*ObeyRestrictions')
            call read_obey_restrictions_load(file_id, config, write_freq, mb)
         case ('*PerturbationsS')
            call read_perturbations_S_load(file_id, config, write_freq)
         case ('*PerturbationsE')
            call read_perturbations_E_load(file_id, config, write_freq)
         case ('*RandomWalk')
            call read_random_walk_load(file_id, config, write_freq)
         case default
            if (keyword(1:11) == '*ImportFile') then
               call read_file_load(file_id, config, write_freq, align)
               call readAlignment(align, config%import_file)
               call load_import_data(config)
            else
               error stop 'read_repetition_block: unknown keyword in *Repetition block'
            end if
         end select

         rep_buf(i)%config     = config
         rep_buf(i)%write_freq = write_freq
      end do
   end subroutine read_repetition_block

end module indr_step_parser
