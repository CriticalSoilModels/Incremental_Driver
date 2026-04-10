module indr_types
    use stdlib_kinds, only: dp
    use indr_constants, only: voigt_len
   implicit none
!    private
!    public:: 

   type StressAlignment
      logical:: active
      character(len=40) :: ImportFileName
      integer:: kblank,nrec,kReversal, ncol
      integer,dimension(100) :: Reversal
      integer,dimension(6):: isig
      real(dp),dimension(6) :: sigFac
   end type StressAlignment
contains


end module indr_types
