module indr_step_params
   use stdlib_kinds, only: dp
   use indr_constants, only: max_fname_len, voigt_len

   implicit none
   private
   public :: step_config_t, material_state_t

   type step_config_t
      integer  :: n_inc, max_iter, ifstress(voigt_len), columns_in_file(7), n_import
      real(dp) :: delta_load_circ(voigt_len), phase0(voigt_len), delta_load(9), &
                  dfgrd0(3,3), dfgrd1(3,3), delta_time, import_factor(7),       &
                  delta_temp
      character(40) :: load_type, coord_sys, exit_cond, import_file
      real(dp), dimension(1:6,1:6) :: cMt, cMe
      real(dp), dimension(1:6)     :: mbinc
      logical :: has_exit_cond
   end type step_config_t

   type material_state_t
      !! Persistent model state carried across increments.
      !! Also used as a per-increment result snapshot (one entry per results(:) element).
      real(dp)              :: sig(6)        !! Cauchy stress [kPa], Voigt
      real(dp)              :: eps(6)        !! total strain [-]
      real(dp), allocatable :: statev(:)     !! model state variables
      real(dp)              :: time(2)       !! [step time, total time] [s]
      real(dp)              :: dt            !! time increment for this snapshot [s]
      real(dp)              :: temp          !! temperature [°C]
      real(dp)              :: F_start(3,3)  !! deformation gradient, start of increment [-]
      real(dp)              :: F_end(3,3)    !! deformation gradient, end of increment [-]
   end type material_state_t

end module indr_step_params
