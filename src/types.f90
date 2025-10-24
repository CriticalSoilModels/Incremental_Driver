module mod_types
    use stdlib_kinds, only: dp
    use mod_constants, only: voight_len
   implicit none
!    private
!    public:: 

   type descriptionOfStep
      integer:: ninc, maxiter, ifstress(voight_len),columnsInFile(7),mImport ! AN 2016
      real(dp) :: deltaLoadCirc(voight_len),phase0(voight_len),deltaLoad(9),    &
         dfgrd0(3,3), dfgrd1(3,3),deltaTime, importFactor(7),&
         deltaTemp                                            ! AN 2023 temperat
      character(40) :: keyword2, keyword3, exitCond,ImportFileName    ! AN 2016
      real(dp),dimension(1:6,1:6) :: cMt, cMe
      real(dp),dimension(1:6) :: mbinc
      logical::existCond                                              ! AN 2016
   end type  descriptionOfStep

   type StressAlignment
      logical:: active
      character(len=40) :: ImportFileName
      integer:: kblank,nrec,kReversal, ncol
      integer,dimension(100) :: Reversal
      integer,dimension(6):: isig
      real(dp),dimension(6) :: sigFac
   end type StressAlignment
contains


end module mod_types
