module mod_loads
   use stdlib_kinds, only: dp
   use mod_inc_driver_funcs, only: ReadStepCommons, splitaLine, parser
   use mod_alignment, only: readAlignment
   use mod_types, only: StressAlignment
   use mod_constants, only: voight_len, max_fname_len


   implicit none
   ! private
   ! public ::
contains

   subroutine read_deformation_gradient_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword, deltaLoad)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword
      real(dp), intent(out) :: deltaLoad(:)

      !Local
      integer :: i
      call ReadStepCommons(test_file_id,ninc,maxiter,&
         deltaTime,deltaTemp,write_freq)    ! AN 2016 ! AN 2023 temperat

      keyword = '*Cartesian'

      do i=1,9
         read(test_file_id,*)  deltaLoad(i)                                      !  dload means total change in the whole step here
      enddo
   end subroutine read_deformation_gradient_load


   subroutine read_circulating_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword, deltaLoad, ifstress, &
      deltaLoadCirc, phase)

      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword
      real(dp), intent(out) :: deltaLoad(:)
      integer, intent(out) :: ifstress(:)
      real(dp), intent(out) :: deltaLoadCirc(:)
      real(dp), intent(out) :: phase(:)

      integer :: i
      call ReadStepCommons(test_file_id,ninc,maxiter,&
         deltaTime,deltaTemp,write_freq)    ! AN 2016 ! AN 2023 temperat


      read(test_file_id,*) keyword  !  = Cartesian or Roscoe  or RoscoeIsomorph or Rendulic
      keyword  = trim(keyword)

      do i=1,6
         read(test_file_id,*) ifstress(i),deltaLoadCirc(i), phase(i), deltaLoad(i)   !  dload means amplitude here
      enddo

   end subroutine read_circulating_load

   subroutine read_linear_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword, ifstress, deltaLoad)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword
      integer, intent(out) :: ifstress(:)
      real(dp), intent(out) :: deltaLoad(:)

      integer :: i
      call ReadStepCommons(test_file_id,ninc,maxiter,&
         deltaTime,deltaTemp,write_freq)    ! AN 2016 ! AN 2023 temperat


      read(test_file_id,*) keyword                                           !  = Cartesian or Roscoe  or RoscoeIsomorph or Rendulic
      keyword  = trim(keyword)
      do i=1,6
         read(1,*) ifstress(i), deltaLoad(i)                           !  dload means total change in the whole step here
      enddo

   end subroutine read_linear_load


   subroutine read_file_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keywords, ifstress, &
      columnsInFile, importFactor, ImportFileName, align)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keywords(3)
      integer, intent(out) :: ifstress(voight_len)
      integer, intent(out) :: columnsInFile(7)
      real(dp), intent(out) :: importFactor(7)
      type(StressAlignment) :: align
      character(len=max_fname_len), intent(out) :: ImportFileName

      !Local variables
      character(len=40) :: load_file_name, aShortLine, leftLine, rightLine
      integer :: i
      logical :: okSplit
      character(len=10) :: keyword2
      integer :: load_file_id

      keyword2 = keywords(2)

      keywords(2) = '*ImportFile'; keyword2 = keyword2(12:)

      call splitaLine( keyword2,'|',ImportFileName, load_file_name,okSplit)

      if(.not.okSplit)    stop 'missing | in line *ImportFile'

      read(load_file_name,*) load_file_id

      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)    ! AN 2016  read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat

      read(test_file_id,'(a)') keywords(3)

      columnsInFile(:) = 0; importFactor(:) = 1

      do i=1,6
         read(test_file_id,'(a)') aShortLine
         call splitaLine(aShortLine,'*',leftLine,rightLine,okSplit )
         read(leftLine,*)  ifstress(i), columnsInFile(i)
         if(okSplit)  read(rightLine,*)  ImportFactor(i)
      enddo

      if(deltaTime <= 0) then
         read(test_file_id,'(a)') aShortLine
         call splitaLine(aShortLine,'*',leftLine,rightLine,okSplit)
         read(leftLine,*) columnsInFile(7)
         if(okSplit)  read(rightLine,*)  ImportFactor(7)
         !         read(test_file_id,*) columnsInFile(7),  ImportFac(7)
      endif !deltaTime

      !**********************************************************
      call readAlignment(align, ImportFileName )
      !**********************************************************
   end subroutine read_file_load


   subroutine read_oedometric_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword2, keyword3, deltaLoad1)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad1

      integer :: i

      keyword2 = '*LinearLoad'
      keyword3 ='*Cartesian'
      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp, write_freq)    ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
      read(test_file_id,*)   deltaLoad1

   end subroutine read_oedometric_load

   subroutine read_oedometric_S1_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword2, keyword3, deltaLoad1, ifstress1)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad1
      integer :: ifstress1

      integer :: i

      keyword2 = '*LinearLoad'
      keyword3 ='*Cartesian'
      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp, write_freq)   ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
      ifstress1 = 1
      read(test_file_id,*)   deltaLoad1

   end subroutine read_oedometric_S1_load

   subroutine read_triaxial_e1_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword2, keyword3, deltaLoad1, ifstress)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad1
      integer :: ifstress(2)

      integer :: i

      keyword2 = '*LinearLoad'
      keyword3 ='*Cartesian'
      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp, write_freq)   ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
      ifstress = 1
      read(test_file_id,*)   deltaLoad1

   end subroutine read_triaxial_e1_load

   subroutine read_triaxial_s1_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword2, keyword3, deltaLoad1, ifstress)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad1
      integer, intent(out) :: ifstress(3)

      integer :: i

      keyword2 = '*LinearLoad'
      keyword3 ='*Cartesian'
      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp, write_freq)   ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
      ifstress = 1
      read(test_file_id,*)   deltaLoad1

   end subroutine read_triaxial_s1_load

   subroutine read_triaxial_ueq_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword2, keyword3, deltaLoad2)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad2


      integer :: i

      keyword2 = '*LinearLoad'
      keyword3 ='*Roscoe'
      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp, write_freq)   ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
      read(test_file_id,*)   deltaLoad2

   end subroutine read_triaxial_ueq_load

   subroutine read_triaxial_uq_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword2, keyword3, deltaLoad2, ifstress2)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad2
      integer, intent(out) :: ifstress2

      keyword2 = '*LinearLoad'
      keyword3 ='*Roscoe'
      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp, write_freq)   ! AN 2016   read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
      read(test_file_id,*)   deltaLoad2
      ifstress2 =1

   end subroutine read_triaxial_uq_load

   subroutine read_pure_relaxation_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword2, keyword3)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword2, keyword3

      keyword2 = '*LinearLoad'
      keyword3 ='*Cartesian'
      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)    ! AN 2016  read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat
   end subroutine read_pure_relaxation_load

   subroutine read_pure_creep_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword2, keyword3, ifstress)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword2, keyword3
      integer, intent(out) :: ifstress(voight_len)

      keyword2 = '*LinearLoad'
      keyword3 ='*Cartesian'
      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)
      ifstress = 1
   end subroutine read_pure_creep_load

   subroutine read_undrained_creep(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword2, keyword3, ifstress)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword2, keyword3
      integer, intent(out) :: ifstress(voight_len)

      keyword2 = '*LinearLoad'
      keyword3 ='*Roscoe'
      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)
      ifstress(1) = 0
      ifstress(2:6) = 1
   end subroutine read_undrained_creep

   subroutine read_obey_restrictions_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword3, ifstress, cMt, cMe, mb, &
      mbinc)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword3
      integer, intent(out) :: ifstress(voight_len)
      real(dp), intent(out) :: cMt(6,6)
      real(dp), intent(out) :: cMe(6,6)
      real(dp), intent(out) :: mb(6), mbinc(6)

      !Local
      integer :: i
      character(260) :: inputline(6)

      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)      ! AN 2016    read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat

      do i=1,6
         read(test_file_id,'(a)')  inputline(i)  !=== a line of form  ''-sd1 + sd2 + 3.0*sd3 = -10  ! a comment '' is expected
         if(index(inputline(i),'=')== 0) stop 'restr without "=" '
      enddo

      call parser(inputline, cMt,cMe,mb )

      mbinc = mb/ninc
      keyword3 ='*Cartesian'
      ifstress(1:6) = 1

   end subroutine read_obey_restrictions_load

   subroutine read_perturbations_S_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword3, deltaLoad1, ifstress)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword3
      integer, intent(out) :: ifstress(voight_len)
      real(dp), intent(out) :: deltaLoad1

      !Local
      integer :: i

      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)       ! AN 2016     read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat

      deltaTemp = 0   ! AN 2023 temperat

      read(test_file_id,*) keyword3  ! = *Rendulic  or *RoscoeIsomorph

      keyword3 = trim( keyword3 )

      if(keyword3 .ne. '*Rendulic' .and. keyword3 .ne. '*RoscoeIsomorph') then
         write(*,*) 'warning: non-Isomorphic perturburbation'
      end if

      read(test_file_id,*)  deltaLoad1

      ifstress(1:6) =  1
   end subroutine read_perturbations_S_load

   subroutine read_perturbations_E_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword3, deltaLoad1, ifstress)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword3
      integer, intent(out) :: ifstress(voight_len)
      real(dp), intent(out) :: deltaLoad1

      !Local
      integer :: i

      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)       ! AN 2016     read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat

      deltaTemp = 0  ! AN 2023 temperat
      read(test_file_id,*) keyword3

      keyword3 = trim( keyword3 )

      if(keyword3 .ne. '*Rendulic' .and. keyword3 .ne. '*RoscoeIsomorph') then
         write(*,*) 'warning: Anisomorphic perturburbation'
      end if

      read(test_file_id,*) deltaLoad1

   end subroutine read_perturbations_E_load

   subroutine read_random_walk_load(test_file_id, ninc,maxiter,&
      deltaTime,deltaTemp,write_freq, keyword3, deltaLoad, ifstress)
      integer, intent(in)  :: test_file_id
      integer, intent(out) :: ninc
      integer, intent(out) :: maxiter
      real(dp), intent(out) :: deltaTime
      real(dp), intent(out) :: deltaTemp
      integer, intent(out) :: write_freq
      character(10), intent(out) :: keyword3
      integer, intent(out) :: ifstress(voight_len)
      real(dp), intent(out) :: deltaLoad(voight_len)

      !Local
      integer :: i

      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)       ! AN 2016     read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat

      deltaTemp = 0  ! AN 2023 temperat
      read(test_file_id,*) keyword3

      keyword3 = trim( keyword3 )

      if(keyword3 .ne. '*Rendulic' .and. keyword3 .ne. '*RoscoeIsomorph') then
         write(*,*) 'warning: Anisomorphic perturburbation'
      end if

      do i=1,6
         read(1,*) ifstress(i),deltaLoad(i)    !  dload means max abs value of to be multiplied by random in (-1,1)
      enddo
   end subroutine read_random_walk_load
end module mod_loads
