
 !! Module contains scripts for reading and writing information to files
module indr_file_io
   use stdlib_kinds, only: dp
   use stdlib_io, only: open, get_line
   use indr_constants, only: max_fname_len, max_mater_len, voigt_len, max_head_len
   use indr_types, only: material_state_t, step_config_t

   implicit none(type, external)
   private
   public :: read_parameter_file, read_init_conditions_file, set_output_name_from_test_file, &
      write_line_output_data, write_output_file_header, write_step_output, load_import_data

contains
   subroutine read_parameter_file(file_name, num_props, prop_vals, &
      material_name)
      !! Read the values from the parameter file
      character(len=max_fname_len), intent(in)  :: file_name
      character(len=max_mater_len), intent(out) :: material_name
      integer, intent(out) ::  num_props
      real(dp), allocatable, intent(out) :: prop_vals(:)

      integer :: file
      integer :: i, iostat

      !! Think about the problem with open and iostat some more.
      file = open(file_name)

      read(file, *) material_name
      read(file, *) num_props

      i = index(material_name, '#')
      if(i == 0) then
         material_name = trim(material_name)
      else
         material_name = material_name(:i-1)
         material_name = trim(material_name)
      endif

      allocate( prop_vals(num_props) )

      do i=1,num_props
         read(file, *) prop_vals(i)
      enddo
      
      close(file)

   end subroutine read_parameter_file

   subroutine read_init_conditions_file(file_name, stress, time, strain, &
      dtime, temperature, state_vars, init_state_vars, state_vars_head, num_state_vars)
      !! Read the inital conditions file and store the variables
      !! Fortran does an implicit save when subroutine variable names are set at
      !! declearation. making that explicit with the
      character(len=max_fname_len), intent(in) :: file_name
      real(dp), intent(out) :: stress(voigt_len)
      real(dp), intent(out) :: time(2)
      real(dp), intent(out) :: strain(voigt_len)

      real(dp), allocatable, intent(out) :: state_vars(:)
      real(dp), allocatable, intent(out) :: init_state_vars(:)
      character(len = 15), allocatable, intent(out) :: state_vars_head(:)

      real(dp), intent(out) :: dtime
      real(dp), intent(out) :: temperature
      integer , intent(out) :: num_state_vars

      ! Local variables
      integer :: i, iostat, file
      character(len = 100) :: aline
      integer :: num_components !! Number of compents in the stress and strain vectors

      ! Set the values to zero
      stress = 0.0_dp
      time   = 0.0_dp
      strain = 0.0_dp
      dtime  = 0.0_dp
      temperature   = 0.0_dp

      file = open(file_name, iostat = iostat)

      read(file,*) num_components

      do i=1,num_components
         read(file,*) stress(i)
      enddo

      read(file,'(a)') aLine          ! AN 2023  aLine may be *temperature= 30.0 or   num_state_vars
      i = index(aLine,'=')            ! AN 2023
      if(i==0) then                   ! AN 2023 '=' is absent so read num_state_vars
         read(aLine,*) num_state_vars ! AN 2023
      else                            ! AN 2023
         aLine = trim(aLine(i+1:))
         read(aLine,*) temperature    ! AN 2023
         read(1,*) num_state_vars     ! AN 2023
      endif                           ! AN 2023

      if(num_state_vars >= 1) then
         allocate(state_vars(num_state_vars)     , source = 0.0_dp)
         allocate(init_state_vars(num_state_vars), source = 0.0_dp)
         allocate(state_vars_head(num_state_vars)); state_vars_head = ' '

         do i= 1, num_state_vars
            read(file, *, end=500) init_state_vars(i)   !
            continue
         enddo

      else
         allocate( state_vars(1) , init_state_vars(1), source = 0.0_dp )             !  AN 2016 formal placeholder not really used
         allocate(state_vars_head(1)); state_vars_head = ' '
         num_state_vars = 1
      endif

500   continue !! TODO: Update this to not use this goto

      ! Store the initial values
      state_vars = init_state_vars

      close(file)
   end subroutine read_init_conditions_file

   subroutine set_output_name_from_test_file(test_file_id,&
      output_file_name, output_file_heading)
      !! Seems like you an overwrite the output file name in the test file. This subroutine
      !! gets the name from the test file
      integer, intent(out) :: test_file_id
      character(len=max_fname_len), intent(inout) :: output_file_name
      character(len=max_head_len), intent(out) :: output_file_heading

      character(len=100) :: aline
      character(len=max_fname_len) :: output_file_name_loc
      integer :: i

      ![4.1] read the outputfilename from test.inp, create/open this file and write the tablehead, heading(if any)  and the first line = initial conditions
      read(test_file_id,'(a)') aLine
      
      i = index(aLine,'#')
      if(i==0) then
         output_file_name_loc=trim(aLine)
         output_file_heading = '#'
      else
         output_file_name_loc=trim(aLine(:i-1))
         output_file_heading = trim( aLine(i+1:) )
      endif

      ! Current default file name is '--' overwrites the default if a name is passed in the test file.
      if(output_file_name == '--') output_file_name = output_file_name_loc

   end subroutine set_output_name_from_test_file

   subroutine write_output_file_header(output_file_id, output_file_heading, num_state_vars)
      integer, intent(out) :: output_file_id
      integer, intent(in) :: num_state_vars
      character(len=max_head_len), intent(out) :: output_file_heading

      character(len=10) :: time_header(2)
      character(len=10) :: strain_header(6)
      character(len=10) :: stress_header(6)
      character(len=15) :: state_var_header(num_state_vars)

      ! Get the output headers
      call get_output_headers(time_header, strain_header, stress_header, state_var_header )

      ! Write the headers
      write(output_file_id,'(a14,500a20)') time_header, strain_header ,stress_header, state_var_header

      if(output_file_heading(1:1) /= '#') write(output_file_id,*) trim(output_file_heading)

   end subroutine write_output_file_header

   subroutine get_output_headers(time_header, strain_header, stress_header, state_var_header)
      !! Writes output file headers and stores them in the *_header arrays
      character(len=10), intent(out) :: time_header(2)
      character(len=10), intent(out) :: strain_header(voigt_len)
      character(len=10), intent(out) :: stress_header(voigt_len)
      character(len=15), intent(out) :: state_var_header(:)
      
      !Local
      integer :: i
      integer :: num_state_vars
      num_state_vars = size(state_var_header)

      do i=1,2
         write(time_header(i),'(a,i0,a)')  'time(',i, ')'
      enddo

      do i=1,voigt_len
         write( strain_header(i), '(a,i0,a)' )   'stran(',i, ')'
         write( stress_header(i), '(a,i0,a)' )  'stress(',i, ')'
      enddo

      do i=1,num_state_vars
         write(state_var_header(i), '(a,i0,a)' ) 'statev(',i,')'
      enddo
   end subroutine get_output_headers

   subroutine write_line_output_data(output_file_id, time, dtime, strain, stress, state_vars)
      integer, intent(in) :: output_file_id
      real(dp), intent(in) :: time(2)
      real(dp), intent(in) :: dtime
      real(dp), intent(in) :: strain(voigt_len)
      real(dp), intent(in) :: stress(voigt_len)
      real(dp), intent(in) :: state_vars(:)

      write(output_file_id,'(500(g17.10,3h    ))') time+(/dtime,dtime/), strain, stress, state_vars

   end subroutine write_line_output_data

   subroutine load_import_data(config)
      !! Read all rows of an *ImportFile data file into config%import_data.
      !!
      !! Opens config%import_file, skips non-numeric header lines, then reads
      !! every subsequent row into a growing buffer. On return:
      !!   config%import_data(i, :) = row i of the file (1-based)
      !!   config%n_import          = number of columns read per row
      !!
      !! Row 1 is the initial state; row kinc+1 is the new state at increment kinc.
      type(step_config_t), intent(inout) :: config

      integer, parameter :: INIT_CAP = 256

      integer  :: fid, iostat, n_rows, n_cols, cap
      character(len=520) :: line
      character(len=1)   :: achar_
      real(dp), allocatable :: buf(:,:), tmp(:,:), row(:)

      n_cols = maxval(config%columns_in_file)
      config%n_import = n_cols

      fid = open(config%import_file)

      ! Skip non-numeric header lines; 'line' ends holding first data row
      do
         read(fid, '(a)', iostat=iostat) line
         if (iostat /= 0) error stop 'load_import_data: error reading file headers'
         line   = adjustl(line)
         achar_ = line(1:1)
         if (index('1234567890+-.', achar_) > 0) exit
      end do

      ! Grow buffer as rows are read
      cap = INIT_CAP
      allocate(buf(cap, n_cols), row(n_cols))
      n_rows = 0

      do
         read(line, *, iostat=iostat) row
         if (iostat /= 0) exit          ! blank or malformed line — stop
         n_rows = n_rows + 1
         if (n_rows > cap) then
            cap = cap * 2
            allocate(tmp(cap, n_cols))
            tmp(1:n_rows-1, :) = buf(1:n_rows-1, :)
            call move_alloc(tmp, buf)
         end if
         buf(n_rows, :) = row
         read(fid, '(a)', iostat=iostat) line   ! next line; EOF exits loop
         if (iostat /= 0) exit
         line = adjustl(line)
      end do

      close(fid)

      allocate(config%import_data(n_rows, n_cols))
      config%import_data = buf(1:n_rows, :)

   end subroutine load_import_data

   subroutine write_step_output(file_id, results, write_freq)
      !! Write increment snapshots to an open output file.
      !!
      !! Writes every write_freq-th entry from results(:).
      !! results(i)%time already holds the end-of-increment time;
      !! write_line_output_data is called with dtime=0 so it is not added again.
      integer,                intent(in) :: file_id
      type(material_state_t), intent(in) :: results(:)
      integer,                intent(in) :: write_freq

      integer :: i, ievery

      ievery = 1
      do i = 1, size(results)
         if (ievery == 1) then
            call write_line_output_data(file_id, results(i)%time, 0.0_dp, &
               results(i)%eps, results(i)%sig, results(i)%statev)
         end if
         ievery = ievery + 1
         if (ievery > write_freq) ievery = 1
      end do

   end subroutine write_step_output

end module indr_file_io
