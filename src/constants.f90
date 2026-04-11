module indr_constants
    !! Module contains constants
    use stdlib_kinds, only: dp
    use stdlib_constants, only: pi_dp
    implicit none
    private
    public :: max_fname_len, max_mater_len, max_lname_len, voigt_len, &
              max_head_len, iter_lower_limit

    integer, parameter :: max_fname_len = 40 !! Max file name length
    integer, parameter :: max_mater_len = 80 !! Max material name length
    integer, parameter :: max_lname_len = 40 !! max length of the load name
    integer, parameter :: voigt_len = 6      !! Length of the Voigt vector
    integer, parameter :: max_head_len = 260 !! max heading length
    integer, parameter :: iter_lower_limit = 5 !! Lowest allowed iterations for stress conditions
    ! real(dp), parameter:: pi = pi_dp
    
end module indr_constants