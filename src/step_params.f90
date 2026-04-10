module indr_step_params
   use stdlib_kinds, only: dp
   use indr_constants, only: max_fname_len, voight_len

   implicit none

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

contains

   pure function set_repetition_params(ninc, maxiter, ifstress, deltaLoadCirc, phase0, &
      deltaLoad, dfgrd0, dfgrd1, deltaTime, keywords, cMe, cMt, mbinc, deltaTemp,&
      exitCond, existCond, ImportFileName, mImport, columnsInFile, importFactor) result(ofStep)
      !! Store the repetition parameters for the last iteration
      integer, intent(in)  :: ninc
      integer, intent(in)  :: maxiter
      integer, intent(in)  :: ifstress(voight_len)
      real(dp), intent(in) :: deltaLoadCirc(6)
      real(dp), intent(in) :: phase0(6)
      real(dp), intent(in) :: deltaLoad(9)
      real(dp), intent(in) :: dfgrd0(3,3)
      real(dp), intent(in) :: dfgrd1(3,3)
      real(dp), intent(in) :: deltaTime
      character(len=10), intent(in) :: keywords(3)
      real(dp), intent(in) :: cMe(6,6)
      real(dp), intent(in) :: cMt(6,6)
      real(dp), intent(in) :: mbinc(6)
      real(dp), intent(in) :: deltaTemp
      character(len=40), intent(in) :: exitCond
      logical, intent(in) :: existCond
      character(max_fname_len), intent(in) :: ImportFileName
      integer, intent(in) :: mImport
      integer, intent(in) :: columnsInFile(7)
      real(dp), intent(in) :: importFactor(7)

      type(descriptionOfStep) :: ofStep

      ofStep%ninc          =    ninc
      ofStep%maxiter       =    maxiter
      ofStep%ifstress      =    ifstress
      ofStep%deltaLoadCirc =    deltaLoadCirc
      ofStep%phase0        =    phase0
      ofStep%deltaLoad     =    deltaLoad
      ofStep%dfgrd0        =    dfgrd0
      ofStep%dfgrd1        =    dfgrd1
      ofStep%deltaTime     =    deltaTime
      ofStep%keyword2      =    keywords(2)
      ofStep%keyword3      =    keywords(3)
      ofStep%cMe           =    cMe
      ofStep%cMt           =    cMt
      ofStep%mbinc         =    mbinc
      ofStep%deltaTemp     =    deltaTemp        ! AN 2023 temperat
      ofStep%exitCond      =    exitCond         ! AN 2016
      ofStep%existCond     =    existCond        ! AN 2016
      ofStep%ImportFileName =   ImportFileName   ! AN 2016
      ofStep%mImport        =   mImport          ! AN 2016
      ofStep%columnsInFile  =   columnsInFile    ! AN 2016    7 integers with numbers of columns  (or value = 0)
      ofStep%importFactor   =   importFactor     !! AN 2016   7 real factors to be multiplied with columns  ! jump over reading, because reading of steps is performed only on the first loop, when iRepetition==1

   end function set_repetition_params

   subroutine get_repetition_params(ofstep, ninc, maxiter, ifstress, deltaLoadCirc, phase0, &
      deltaLoad, dfgrd0, dfgrd1, deltaTime, keywords, cMe, cMt, mbinc, deltaTemp,&
      exitCond, existCond, ImportFileName, mImport, columnsInFile, importFactor) 
      !! Get the repetition parameters stored in the passed step params object
      type(descriptionOfStep), intent(in):: ofStep
      integer, intent(out)  :: ninc
      integer, intent(out)  :: maxiter
      integer, intent(out)  :: ifstress(voight_len)
      real(dp), intent(out) :: deltaLoadCirc(6)
      real(dp), intent(out) :: phase0(6)
      real(dp), intent(out) :: deltaLoad(9)
      real(dp), intent(out) :: dfgrd0(3,3)
      real(dp), intent(out) :: dfgrd1(3,3)
      real(dp), intent(out) :: deltaTime
      character(len=10), intent(out) :: keywords(3)
      real(dp), intent(out) :: cMe(6,6)
      real(dp), intent(out) :: cMt(6,6)
      real(dp), intent(out) :: mbinc(6)
      real(dp), intent(out) :: deltaTemp
      character(len=40), intent(out) :: exitCond
      logical, intent(out) :: existCond
      character(max_fname_len), intent(out) :: ImportFileName
      integer, intent(out) :: mImport
      integer, intent(out) :: columnsInFile(7)
      real(dp), intent(out) :: importFactor(7)

      ninc           = ofStep%ninc             
      maxiter        = ofStep%maxiter          
      ifstress       = ofStep%ifstress         
      deltaLoadCirc  = ofStep%deltaLoadCirc    
      phase0         = ofStep%phase0           
      deltaLoad      = ofStep%deltaLoad        
      dfgrd0         = ofStep%dfgrd0           
      dfgrd1         = ofStep%dfgrd1           
      deltaTime      = ofStep%deltaTime        
      keywords(2)    = ofStep%keyword2         
      keywords(3)    = ofStep%keyword3         
      cMe            = ofStep%cMe              
      cMt            = ofStep%cMt              
      mbinc          = ofStep%mbinc            
      deltaTemp      = ofStep%deltaTemp        ! AN 2023 temperat
      exitCond       = ofStep%exitCond         ! AN 2016
      existCond      = ofStep%existCond        ! AN 2016
      ImportFileName = ofStep%ImportFileName   ! AN 2016
      mImport        = ofStep%mImport          ! AN 2016
      columnsInFile  = ofStep%columnsInFile    ! AN 2016    7 integers with numbers of columns  (or value = 0)
      importFactor   = ofStep%importFactor     !! AN 2016   7 real factors to be multiplied with columns  ! jump over reading, because reading of steps is performed only on the first loop, when iRepetition==1
   end subroutine get_repetition_params
end module indr_step_params
