module mod_command_line
   implicit none(type, external)
   private
   public :: set_inputs

contains
   subroutine get_input_vals_from_cmd(parametersfilename, initialconditionsfilename, &
      testfilename, outputfilename, verbose)
      !! Gets the names/values of the input file names and verbose setting from the commandline

      character(len=40), intent(out) :: parametersfilename
      character(len=40), intent(out) :: initialconditionsfilename
      character(len=40), intent(out) :: testfilename
      character(len=40), intent(out) :: outputfilename
      logical :: verbose

      ! Local variables
      integer :: iarg,narg, is, iargc
      integer, parameter :: argLength=40
      character(argLength) :: anArgument, argType, argValue

      narg = iargc()
      do iarg = 1,narg
         call getarg(iarg,anArgument)
         is = index(anArgument,'=')
         if(is == 0) stop 'error: a command line argument without "=" '

         argType = anArgument(:is-1)
         argValue =  anArgument(is+1:)

         select case (argType)
          case ('param')
            parametersfilename = argValue
          case ('ini')
            initialconditionsfilename = argValue
          case ('test')
            testfilename = argValue
          case ('out')
            outputfilename = argValue
          case ('verbose')
            if (argValue == 'true') verbose = .true.
            if (argValue == 'false') verbose = .false.
         end select

      enddo
   end  subroutine get_input_vals_from_cmd

   subroutine set_default_file_names(parametersfilename, initialconditionsfilename, &
      testfilename, outputfilename)
      character(len=40), intent(out) :: parametersfilename
      character(len=40), intent(out) :: initialconditionsfilename
      character(len=40), intent(out) :: testfilename
      character(len=40), intent(out) :: outputfilename

      parametersfilename = 'parameters.inp'
      initialconditionsfilename  = 'initialconditions.inp'
      testfilename = 'test.inp'
      outputfilename = '--'
   end subroutine set_default_file_names

   subroutine set_inputs(parametersfilename, initialconditionsfilename, &
      testfilename, outputfilename, verbose)
      !! Sets the input file names and verbose setting for the model
      character(len=40), intent(out) :: parametersfilename
      character(len=40), intent(out) :: initialconditionsfilename
      character(len=40), intent(out) :: testfilename
      character(len=40), intent(out) :: outputfilename
      logical          , intent(out) :: verbose


      call set_default_file_names(parametersfilename, initialconditionsfilename, &
         testfilename, outputfilename)

      ! Set the default value for verbose
      verbose = .false.

      call get_input_vals_from_cmd(parametersfilename, initialconditionsfilename, &
         testfilename, outputfilename, verbose)

   end subroutine set_inputs
end module mod_command_line
