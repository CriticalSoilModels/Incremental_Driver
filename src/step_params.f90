module indr_step_params
   use stdlib_kinds, only: dp
   use indr_constants, only: max_fname_len, voigt_len

   implicit none

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

end module indr_step_params
