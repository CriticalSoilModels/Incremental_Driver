module indr_loads
   use stdlib_kinds, only: dp
   use indr_parser, only: ReadStepCommons, splitaLine, PARSER
   use indr_linalg, only: inv33
   use indr_alignment, only: readAlignment
   use indr_types, only: StressAlignment
   use indr_constants, only: voigt_len, max_fname_len, max_lname_len
   use indr_step_params, only: step_config_t

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
      character(max_lname_len), intent(out) :: keyword
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
         read(test_file_id,*) ifstress(i), deltaLoad(i)                !  dload means total change in the whole step here
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
      integer, intent(out) :: ifstress(voigt_len)
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
      character(max_lname_len), intent(inout) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad1

      integer :: i

      keyword2 = adjustl(trim('*LinearLoad'))
      keyword3 = adjustl(trim('*Cartesian'))

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
      character(max_lname_len), intent(inout) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad1
      integer :: ifstress1

      integer :: i

      keyword2 = adjustl(trim('*LinearLoad'))
      keyword3 =adjustl(trim('*Cartesian'))
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
      character(max_lname_len), intent(inout) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad1
      integer :: ifstress(2)

      integer :: i

      keyword2 = adjustl(trim('*LinearLoad'))
      keyword3 = adjustl(trim('*Cartesian'))
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
      character(max_lname_len), intent(inout) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad1
      integer, intent(out) :: ifstress(3)

      integer :: i

      keyword2 = adjustl(trim('*LinearLoad'))
      keyword3 = adjustl(trim('*Cartesian'))
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
      character(max_lname_len), intent(inout) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad2


      integer :: i

      keyword2 = adjustl(trim('*LinearLoad'))
      keyword3 = adjustl(trim('*Roscoe'))
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
      integer, intent(inout) :: write_freq
      character(max_lname_len), intent(inout) :: keyword2, keyword3
      real(dp), intent(out) :: deltaLoad2
      integer, intent(out) :: ifstress2

      keyword2 = adjustl(trim('*LinearLoad'))
      keyword3 = adjustl(trim('*Roscoe'))
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
      character(max_lname_len), intent(inout) :: keyword2, keyword3

      keyword2 = adjustl(trim('*LinearLoad'))
      keyword3 = adjustl(trim('*Cartesian'))
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
      character(max_lname_len), intent(inout) :: keyword2, keyword3
      integer, intent(out) :: ifstress(voigt_len)

      keyword2 = adjustl(trim('*LinearLoad'))
      keyword3 =adjustl(trim('*Cartesian'))
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
      character(max_lname_len), intent(inout) :: keyword2, keyword3
      integer, intent(out) :: ifstress(voigt_len)

      keyword2 = adjustl(trim('*LinearLoad'))
      keyword3 = adjustl(trim('*Roscoe'))
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
      character(max_lname_len), intent(out) :: keyword3
      integer, intent(out) :: ifstress(voigt_len)
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
      keyword3 = adjustl(trim('*Cartesian'))
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
      character(max_lname_len), intent(out) :: keyword3
      integer, intent(out) :: ifstress(voigt_len)
      real(dp), intent(out) :: deltaLoad1

      !Local
      integer :: i

      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)       ! AN 2016     read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat

      deltaTemp = 0   ! AN 2023 temperat

      read(test_file_id,*) keyword3  ! = *Rendulic  or *RoscoeIsomorph

      keyword3 = adjustl(trim(keyword3))

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
      character(max_lname_len), intent(out) :: keyword3
      integer, intent(out) :: ifstress(voigt_len)
      real(dp), intent(out) :: deltaLoad1

      !Local
      integer :: i

      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)       ! AN 2016     read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat

      deltaTemp = 0  ! AN 2023 temperat
      read(test_file_id,*) keyword3

      keyword3 = adjustl(trim(keyword3 ))

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
      character(max_lname_len), intent(out) :: keyword3
      integer, intent(out) :: ifstress(voigt_len)
      real(dp), intent(out) :: deltaLoad(voigt_len)

      !Local
      integer :: i

      call ReadStepCommons(test_file_id,ninc,maxiter,deltaTime,deltaTemp,write_freq)       ! AN 2016     read(1,*) ninc, maxiter, deltaTime ! AN 2023 temperat

      deltaTemp = 0  ! AN 2023 temperat
      read(test_file_id,*) keyword3

      keyword3 = adjustl(trim(keyword3 ))

      if(keyword3 .ne. '*Rendulic' .and. keyword3 .ne. '*RoscoeIsomorph') then
         write(*,*) 'warning: Anisomorphic perturburbation'
      end if

      do i=1,6
         read(test_file_id,*) ifstress(i),deltaLoad(i)    !  dload means max abs value, multiplied by random in (-1,1)
      enddo
   end subroutine read_random_walk_load

  

   ! Converts a step description into per-increment ddstress/dstran.
   ! Called once per increment (not once per step).
   subroutine get_increment(config, time, dtime, ddstress, dstran, dTemp, Qb33, dfgrd0, dfgrd1, drot)
      implicit none
      type(step_config_t), intent(in)    :: config
      real(dp),            intent(in)    :: time(2)
      real(dp),            intent(out)   :: dtime, ddstress(6), dstran(6), Qb33(3,3), dTemp
      real(dp),            intent(inout) :: dfgrd0(3,3), dfgrd1(3,3), drot(3,3)

      real(dp), parameter :: Pi = 3.1415926535897932385_dp
      real(dp), parameter :: delta(3,3) = reshape([1,0,0,0,1,0,0,0,1], [3,3])
      real(dp) :: Fb(3,3), Fbb(3,3), dFb(3,3), aux33(3,3), dLb(3,3), depsb(3,3), dOmegab(3,3)
      real(dp) :: wd(6), w0(6), t, arandom
      integer  :: i

      associate(load_type       => config%load_type,       &
                delta_time      => config%delta_time,      &
                n_inc           => config%n_inc,           &
                ifstress        => config%ifstress,        &
                delta_load      => config%delta_load,      &
                delta_load_circ => config%delta_load_circ, &
                phase0          => config%phase0,          &
                delta_temp      => config%delta_temp)

      dtime    = delta_time / n_inc
      dTemp    = delta_temp / n_inc
      dstran   = 0.0_dp
      ddstress = 0.0_dp
      Qb33     = delta
      drot     = delta
      dfgrd0   = delta
      dfgrd1   = delta

      if (load_type == '*LinearLoad') then
         do i = 1, 6
            if (ifstress(i) == 1) ddstress(i) = delta_load(i) / n_inc
            if (ifstress(i) == 0) dstran(i)   = delta_load(i) / n_inc
         end do
      end if

      if (load_type == '*DeformationGradient') then
         Fb = reshape([delta_load(1), delta_load(5), delta_load(7), &
                       delta_load(4), delta_load(2), delta_load(9), &
                       delta_load(6), delta_load(8), delta_load(3)], [3,3])
         Fbb    = delta + (Fb - delta) * (time(1) / delta_time)
         dfgrd0 = Fbb
         dFb    = (Fb - delta) / n_inc
         aux33  = Fbb + dFb / 2.0_dp
         dfgrd1 = Fbb + dFb
         aux33  = inv33(aux33)
         dLb    = matmul(dFb, aux33)
         depsb  = 0.5_dp * (dLb + transpose(dLb))
         dstran = [depsb(1,1), depsb(2,2), depsb(3,3), &
                   2.0_dp*depsb(1,2), 2.0_dp*depsb(1,3), 2.0_dp*depsb(2,3)]
         dOmegab = 0.5_dp * (dLb - transpose(dLb))
         aux33   = inv33(delta - 0.5_dp*dOmegab)
         Qb33    = matmul(aux33, delta + 0.5_dp*dOmegab)
         drot    = Qb33
      end if

      if (load_type == '*CirculatingLoad') then
         wd = 2.0_dp * Pi / delta_time
         w0 = phase0
         t  = time(1) + dtime / 2.0_dp
         do i = 1, 6
            if (ifstress(i) == 1) ddstress(i) = dtime*delta_load_circ(i)*wd(i)*cos(wd(i)*t + w0(i)) + delta_load(i)/n_inc
            if (ifstress(i) == 0) dstran(i)   = dtime*delta_load_circ(i)*wd(i)*cos(wd(i)*t + w0(i)) + delta_load(i)/n_inc
         end do
      end if

      if (load_type == '*PerturbationsS') then
         ddstress(1) = delta_load(1) * cos(time(1) * 2.0_dp * Pi / delta_time)
         ddstress(2) = delta_load(1) * sin(time(1) * 2.0_dp * Pi / delta_time)
      end if

      if (load_type == '*PerturbationsE') then
         dstran(1) = delta_load(1) * cos(time(1) * 2.0_dp * Pi / delta_time)
         dstran(2) = delta_load(1) * sin(time(1) * 2.0_dp * Pi / delta_time)
      end if

      if (load_type == '*RandomWalk') then
         call random_seed()
         do i = 1, 6
            call random_number(arandom)
            if (ifstress(i) == 1) ddstress(i) = 2.0_dp*(arandom - 0.5_dp)*delta_load(i)
            if (ifstress(i) == 0) dstran(i)   = 2.0_dp*(arandom - 0.5_dp)*delta_load(i)
         end do
      end if

      end associate
   end subroutine get_increment

end module indr_loads
