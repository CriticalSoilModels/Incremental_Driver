module indr_test_parser
   !! Reads an entire test.inp file and returns all step records.
   !!
   !! Phase 6: supports all load types except *Repetition.
   use stdlib_kinds,      only: dp
   use indr_step_params,  only: step_config_t
   use indr_types,        only: StressAlignment
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

   implicit none
   private
   public :: step_record_t, parse_test_file

   integer, parameter :: MAX_BUF = 1000  !! Maximum number of steps that can be buffered

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

      type(step_record_t) :: buf(MAX_BUF)
      character(len=260)  :: line
      character(len=40)   :: keyword, exit_cond, kw40
      logical             :: has_exit
      type(step_config_t) :: config
      integer             :: write_freq
      real(dp)            :: mb(6)
      real(dp), parameter :: delta(3,3) = reshape([1,0,0,0,1,0,0,0,1], [3,3])

      ! Initialise outputs
      align%active = .false.
      n_steps = 0

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
            error stop 'parse_test_file: *Repetition is not supported in Phase 6'

         case default
            ! Check for *ImportFile prefix (may have |filename|id appended)
            if (keyword(1:11) == '*ImportFile') then
               call read_file_load(file_id, config, write_freq, align)
               ! read_file_load has align intent(in) so we call readAlignment directly
               call readAlignment(align, config%import_file)
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
      !! Append one step_record_t to the buffer; stops on overflow.
      type(step_record_t), intent(inout) :: buf(MAX_BUF)
      integer,             intent(inout) :: n
      type(step_config_t), intent(in)    :: config
      integer,             intent(in)    :: write_freq

      n = n + 1
      if (n > MAX_BUF) error stop 'parse_test_file: step buffer overflow (> 1000 steps)'
      buf(n)%config     = config
      buf(n)%write_freq = write_freq
   end subroutine store

end module indr_test_parser
