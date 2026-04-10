module indr_value_checks
   use stdlib_kinds, only: dp
   use stdlib_optval, only: optval

   implicit none

contains
   function set_zero_with_tol(arr, ext_tol) result(arr_zero)
      !! prevents fortran error write 1.3E-391
      !! Check if I actually need this.
      real(dp), intent(in) :: arr(:)
      real(dp), intent(in), optional :: ext_tol
      real(dp)  :: arr_zero(size(arr))
      !Local
      real(dp) :: tol

      tol = optval(ext_tol, default = 1.0e-99_dp)

      arr_zero = arr

      where(abs(arr) < tol) arr_zero = 0.0_dp
   end function set_zero_with_tol

   subroutine check_stress_inc_size(a_dstress, u_dstress)
      !! Writes error message to the terminal of the norm squared of the a-stress increment
      !! and u-stress increment are larger than a set tolerance
      real(dp), intent(in) :: a_dstress(:)
      real(dp), intent(in) :: u_dstress(:)
      ! Local variables
      real(dp) :: a_norm_sq, u_norm_sq

      a_norm_sq = norm2(a_dstress)**2
      u_norm_sq = norm2(u_dstress)**2

      if((a_norm_sq>1.d-10 .and. u_norm_sq/a_norm_sq > 1.0d-2) .or. (a_norm_sq<1.d-10 .and. u_norm_sq > 1.0d-12)  ) then
         write(*,*) 'I cannot apply the prescribed stress components,'// &
            '||u_dstress|| too large.'                 ! check  the Rosc.stress error < toler
      endif
   end subroutine check_stress_inc_size

end module indr_value_checks
