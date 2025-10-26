module mod_run_model

   use stdlib_kinds, only: dp
   use stdlib_io, only: open
   use mod_inc_driver_funcs, only: splitaLine, ReadStepCommons, PARSER, get_increment,USOLVER, EXITNOW

   use mod_types   , only: StressAlignment
   use mod_step_params, only: descriptionOfStep, get_repetition_params, set_repetition_params
   use mod_matrices, only: MRoscI, MRoscImt, MRendul, MRendulmT, MRosc, MRoscmT, MCart, MCartmT

   use mod_command_line, only: set_inputs
   use mod_file_io, only: read_parameter_file, read_init_conditions_file, set_output_name_from_test_file, &
      write_line_output_data, write_output_file_header
   use mod_alignment, only: readAlignment, tryAlignStress
   use mod_loads
   use mod_maps
   use mod_constants, only: iter_lower_limit
   use mod_value_checks, only: set_zero_with_tol, check_stress_inc_size

   implicit none(type, external)

   ABSTRACT INTERFACE

      subroutine umat_interface(stress,statev,ddsdde,sse,spd,scd, &
         rpl,ddsddt,drplde,drpldt, &
         stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
         ndi,nshr,ntens,nstatev,props,nprops,coords,drot,pnewdt, &
         celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

         implicit none

         character*80 cmname
         integer :: ntens,nstatev,nprops,ndi,nshr,noel, &
            npt,layer,kspt,kstep,kinc
         real(8) :: sse,spd,scd,rpl,drpldt,dtime,temp,dtemp, &
            pnewdt,celent
         real(8) :: stress(ntens),statev(nstatev), &
            ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens), &
            stran(ntens),dstran(ntens),time(2),predef(1),dpred(1), &
            props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
      end subroutine umat_interface
   END INTERFACE

contains
   subroutine run_model(UMAT)
      procedure(umat_interface), intent(in) :: UMAT

      real(dp), parameter          :: delta(3,3) = reshape([1,0,0,0,1,0,0,0,1],[3,3])
      integer, parameter           :: ntens=6, ndi=3,nshr=3,ncrds=3 ! same ntens as in SOLVER
      integer, parameter           :: noel=1 , npt=1,layer=1,kspt=1,lrebar=1
      character(len=80), parameter :: rebarn ='xxx'

      integer :: nstatv,nprops
      integer :: kinc,i
      integer :: test_file_id, output_file_id, import_file_id !! Variable to hold file id when it's opened (Will be deprecated)
      integer :: iostat

      real(dp) :: dtime,temp,dtemp,sse,spd,scd,rpl,drpldt,pnewdt,celent
      real(dp) :: stress(ntens),&
         ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),&
         stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),&
         coords(ncrds),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      character(len=80) ::cmname
      character(len=1) :: aChar                                       ! AN 2016

      character(len=40):: keywords(10), outputfilename,&
         parametersfilename,&
         initialconditionsfilename, testfilename,&
         exitCond, ImportFileName

      character(len=260) ::  inputline(6), aLine, heading
      character(len=520) :: hugeLine

      character(len=10):: timeHead(2), stranHead(6), stressHead(6)      ! AN 2016
      character(len=15), allocatable :: statevHead(:)                   ! AN 2016

      logical :: verbose
      ! logical :: EXITNOW, existCond,okSplit                             ! AN 2016 ! WaveHello: ExitNow (bool) conflicts with the function
      logical :: existCond,okSplit
      real(dp), dimension(6,6)  :: cMt , cMe
      real(dp), dimension(6)  :: mb, mbinc


      integer :: mImport, columnsInFile(7),every,ievery                 ! AN 2016
      real(dp) ::  importFactor(7)
      real(dp),dimension(20) :: oldState, newState,dState
      real(dp), allocatable :: props(:), statev(:), r_statev(:)

      real(dp),dimension(3,3):: Qb33,eps33,T33

      integer:: ifstress(ntens), maxiter, ninc,kiter, ikeyword, &
         iRepetition, nRepetitions, kStep,iStep,nSteps,ntens_in

      real(dp):: r_stress(ntens),a_dstress(ntens),u_dstress(ntens),&
         stress_Rosc(ntens),r_stress_Rosc(ntens),           &
         ddstress(ntens), c_dstran(ntens) ,                 &
         deltaLoadCirc(6),phase0(6),deltaLoad(9),           &
         dstran_Cart(6), ddsdde_bar(6,6), deltaTime,        &
         deltaTemp

      real(dp),dimension(1:6,1:6)::M,MmT      !  currrent $\cM$ and $\cM^{-T}$  for a given iStep

      type(StressAlignment) :: align
      type(descriptionOfStep) :: ofStep(30)            !  stores descriptions of up to 30 steps which are repeated


      ! [1]  Set the filenames and the verose seting
      call set_inputs(parametersfilename, initialconditionsfilename, &
         testfilename, outputfilename, verbose)

      ! Get the number of properties, the property values and the material name
      call read_parameter_file(parametersfilename, nprops, props, cmname)

      ! Set default values and read initial conditions
      call read_init_conditions_file(initialconditionsfilename, stress, time, stran,&
         dtime, temp, statev, r_statev, statevHead, nstatv)

      ! TODO: File is open for awhile doesn't get closed until later
      test_file_id = open(testfilename)

      call set_output_name_from_test_file(test_file_id, outputfilename, heading)

      ! TODO: File is open for awhile doesn't get closed until later
      output_file_id = open(outputfilename, "w")

      call write_output_file_header(output_file_id, heading, nstatv)
      call write_line_output_data(output_file_id, time, dtime, stran, stress, statev)

![4.2]  loop over keywords(1) unless keyword(1) = *Repetition  it is copied to keyword(2) which is the true type of loading
      kStep = 0  ! kStep = counter over all steps whereas  iStep = counter over steps within a *Repetition
      do_keyword: do ikeyword=1,10000

         read(test_file_id,'(a)') keywords(1)

         if (iostat < 0) error stop "Reached the end of the the test file"//testfilename

         keywords(1) = trim( keywords(1) )
         if(keywords(1) == '*Repetition') then
            read(test_file_id,*) nSteps, nRepetitions
         else
            nRepetitions=1
            nSteps=1
            keywords(2) = keywords(1)
         endif

         do_repet: do iRepetition  = 1,nRepetitions
            do_istep: do iStep = 1,nSteps
               kStep = kStep + 1

               ! ! write to the screen before the first increment
               ! write(*,'(12H ikeyword = ,i3, 8H kstep = ,i3,7H kinc = ,i5, 9H kiter = ,i2, 8H TEMP = ,f9.4, 8H TIME = ,f9.4)'), ikeyword , kstep, kinc, kiter ,TEMP,TIME(1)                    ! AN 2023 temperat
               ! write to the screen before the first increment
               ! WaveHello: Replacing this format
               write(*,'(A,I3,A,I3,A,I5,A,I2,A,F9.4,A,F9.4)') &
                  ' ikeyword = ', ikeyword, &
                  ' kstep = ', kstep, &
                  ' kinc = ', kinc, &
                  ' kiter = ', kiter, &
                  ' TEMP = ', TEMP, &
                  ' TIME = ', TIME(1)

               if(iRepetition > 1) then  ! while repeating  recall the loading parameters of the repeated step read in during the first iRepetition
                  ! jump over reading, because reading of steps is performed only on the first loop, when iRepetition==1
                  call get_repetition_params(ofstep(iStep), &
                     ninc, maxiter, ifstress, deltaLoadCirc, phase0, &
                     deltaLoad, dfgrd0, dfgrd1, deltaTime, keywords, cMe, cMt, mbinc, deltaTemp,&
                     exitCond, existCond, ImportFileName, mImport, columnsInFile, importFactor)
                  ! AN 2016   7 real factors to be multiplied with columns
               endif

               if(keywords(1) == '*Repetition') read(1,'(a)') keywords(2)          ! = *LinearLoad  or *CirculatingLoad or *ObeyRestrictions...
               ! otherwise keywords(2) = keywords(1)

               call splitaLine(keywords(2),'?', keywords(2), exitCond, existCond)  ! AN 2016 look for exit condition in keywords(2)

               keywords(2)  = trim(keywords(2)) ! Trim the test name

               ifstress(:)=0                                                      ! default strain control
               deltaLoadCirc(:)=0.0d0                                             ! default zero step increment
               phase0(:)=0.0d0                                                    ! default no phase shift
               deltaLoad(:) = 0.0d0
               deltaTemp = 0.0d0    !  AN 2023  temperature increase per step
               dfgrd0 = delta
               dfgrd1 = delta

               if(keywords(2) == '*DeformationGradient') then

                  call read_deformation_gradient_load(test_file_id, ninc, maxiter, deltaTime, &
                     deltaTemp, every, keywords(3), deltaLoad)


               else if (keywords(2)=='*CirculatingLoad') then
                  call read_circulating_load(test_file_id, ninc, maxiter, deltaTime, &
                     deltaTemp, every, keywords(3), deltaLoad, ifstress, deltaLoadCirc, &
                     phase0)

               else if(keywords(2)=='*LinearLoad') then
                  call read_linear_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, every, &
                     keywords(3), ifstress, deltaLoad)

               else if(keywords(2)(1:11) == '*ImportFile') then
                  !! TODO: Check that this isn't broken. The load file id thing is funky.
                  !! Not sure if the id is required in the main program. I don't think it's required?
                  !! This is probably broken
                  call read_file_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords, ifstress, columnsInFile, importFactor, &
                     ImportFileName, align)

               else if(keywords(2) == '*OedometricE1') then
                  call read_oedometric_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(2), keywords(3), deltaload(1))
                  keywords(3) = trim(keywords(3))
               else if(keywords(2) == '*OedometricS1') then
                  call read_oedometric_S1_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(2), keywords(3), deltaLoad(1), ifstress(1))

               else if(keywords(2) == '*TriaxialE1') then
                  call read_triaxial_e1_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(2), keywords(3), deltaLoad(1), ifstress(2:3))

               else if(keywords(2) == '*TriaxialS1') then
                  call read_triaxial_s1_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(2), keywords(3), deltaLoad(1), ifstress(1:3))

               else if(keywords(2) == '*TriaxialUEq') then
                  call read_triaxial_ueq_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(2), keywords(3), deltaLoad(2))

               else if(keywords(2) == '*TriaxialUq') then
                  call read_triaxial_uq_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(2), keywords(3), deltaLoad(2), ifstress(2))

               else if(keywords(2) == '*PureRelaxation') then
                  call read_pure_relaxation_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(2), keywords(3))

               else if(keywords(2) == '*PureCreep') then
                  call read_pure_creep_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(2), keywords(3), ifstress)

               else if(keywords(2) == '*UndrainedCreep') then
                  call read_undrained_creep(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(2), keywords(3), ifstress)

               else if(keywords(2) == '*ObeyRestrictions') then  ! ======================= *ObeyRestrictions ==================================
                  call read_obey_restrictions_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(3), ifstress, cMt, cMe, mb, mbinc)

               else if(keywords(2) == '*PerturbationsS') then
                  call read_perturbations_S_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(3), deltaLoad(1), ifstress)

               else if(keywords(2) == '*PerturbationsE') then
                  call read_perturbations_E_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(3), deltaLoad(1), ifstress)

               else if(keywords(2) == '*RandomWalk') then                         ! AN 2019
                  call read_random_walk_load(test_file_id, ninc, maxiter, deltaTime, deltaTemp, &
                     every, keywords(3), deltaLoad, ifstress)
               else
                  if(keywords(2) == '*End') stop '*End encountered in test.inp'
                  write(*,*) 'error: unknown keywords(2)=',keywords(2)
                  stop 'stopped by unknown keyword(2) in test.inp'
               end if

               keywords(3) = trim(keywords(3))

               if(keywords(1) == '*Repetition' .and. iRepetition == 1) then      !  remember the description of step for the next repetition
                  ofStep(istep) = set_repetition_params(ninc, maxiter, ifstress, deltaLoadCirc, phase0, &
                     deltaLoad, dfgrd0, dfgrd1, deltaTime, keywords, cMe, cMt, mbinc, deltaTemp, exitcond, &
                     existcond, ImportFileName, mImport, columnsInFile, importFactor )
               endif

               ! If the choosen load is stress controlled make sure at least iter_lower_limit number of iterations is done
               ! A guess at the correct strain has to be made and the stress has to be converged to
               if(any(ifstress==1)) maxiter = max(maxiter,iter_lower_limit)
               if(all(ifstress==0) .and. keywords(2) .ne. '*ObeyRestrictions') maxiter = 1  ! no iterations are necessary

               ! start the current step with zero-load call of umat() just to get the stiffness
               dstran(:)=0
               dtime=0
               dtemp=0
               kinc=0
               r_statev(:)=statev(:);  r_stress(:)=stress(:)     ! AN 21.06.2017 remember the initial state and stress
               !=== first call umat with dstrain=0 dtime=0 just for stiffness (=jacobian ddsdde)
               call  UMAT(stress,statev,ddsdde,sse,spd,scd,                       &
                  rpl,ddsddt,drplde,drpldt,                               &
                  stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
                  ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,  &
                  celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,0,kinc)   !=== some constitutive models require kStep=0 other do not

               statev(:)=r_statev(:);  stress(:)=r_stress(:)   !  AN 21.06.2017 recover stress and state  although the ZERO call of umat should not modify them

               select case( keywords(3) )
                case('*Cartesian' ) ;      M =  MCart ; MmT = MCartmT
                case('*Roscoe')     ;      M = MRosc  ; MmT = MRoscmT
                case('*RoscoeIsomorph');   M = MRoscI ; MmT = MRoscImT
                case('*Rendulic')      ;   M = MRendul; MmT = MRendulmT
                case default ;   write(*,*) 'Unknown keyword ', keywords(3)
                  stop  ' stopped by unknown keywords(3) in test.inp'
               end select

               if(keywords(2) == '*ImportFile' ) then                           ! AN 2016
                  import_file_id = open(ImportFileName)
                  do                                                             ! AN 2016
                     read(import_file_id,'(a)',iostat = iostat) hugeLine;         ! AN 2016

                     if (iostat /= 0) error stop 'Error reading ImportFile in the first non-numeric records '
                     hugeLine= adjustL(hugeLine) ; aChar = hugeLine(1:1)         ! AN 2016
                     if(index('1234567890+-.',aChar) > 0) exit                  ! preceding non-numeric lines in ImportFile will be ignored
                  enddo

                  read(hugeLine,*,iostat=iostat) oldState(1:mImport)
                  if (iostat /= 0) error stop 'Error reading ImportFile in the first numeric record '
               endif

               ievery=1
               do_kinc: do kinc=1,ninc

                  if(keywords(2) == '*ImportFile') then                             ! AN 2016
                     deltaTemp = 0                                                  ! AN 2023 temperat
                     dTemp = 0                                                      ! AN 2023 temperat
                     read(import_file_id,*,iostat=iostat)  newState(1:mImport)                         ! AN 2016
                     if(iostat > 0) then                                                  ! AN 2016
                        write(*,*) 'error Import file',ImportFileName, 'line=', kinc+1
                        error stop                                                            ! AN 2016
                     else                                                                                                         ! AN 2016
                        close(import_file_id)                                       ! AN 2016
                        write(*,*) 'finished reading file', ImportFileName          ! AN 2016
                        exit                                                        ! AN 2016
                     endif                                                           ! AN 2016

                     dState = newState(:) -  oldState(:)                             ! AN 2016
                     do i=1,6                                                        ! AN 2016
                        if  (columnsInFile(i) == 0) cycle                               ! AN 2016
                        if (ifstress(i)==1  )  ddstress(i)= dState(columnsInFile(i))* ImportFactor(i)
                        if (ifstress(i)==0)    dstran(i)  = dState(columnsInFile(i))* ImportFactor(i)
                     enddo

                     if(columnsInFile(7)/= 0) deltaTime= dState(columnsInFile(7)) * ImportFactor(i)   ! AN 2016
                     dtime = deltaTime                                             ! AN 2016
                     oldState(:) = newState(:)                                    ! AN 2016
                  endif                                                             ! AN 2016

                  if(keywords(2) /= '*ImportFile') then                            ! AN 2016

                     call get_increment(keywords, time, deltaTime, ifstress, ninc,  &    ! get inc. in terms of Rosc. variables
                        deltaLoadCirc,phase0,deltaLoad,deltaTemp,&
                        dtime, ddstress,  dstran, dTemp,  Qb33,  &    ! AB 2023 deltaTemp and dTemp added
                        dfgrd0, dfgrd1,drot )   ! to be called in each increment
                  endif


                  a_dstress(:)= 0.0d0    ! approximated Roscoe's dstress
                  r_statev(:)=statev(:)  ! remember the initial state and stress till the iteration is completed
                  r_stress(:)= stress    ! remembered Cartesian stress

                  do_kiter: do kiter=1, maxiter  !--------Equilibrium Iteration--------
                     c_dstran(:) = 0

                     if(keywords(2)== '*ObeyRestrictions'  ) then  ! ======================= ObeyRestrictions =========
                        ddsdde_bar = matmul(cMt,ddsdde) + cMe
                        u_dstress = - matmul(cMt,a_dstress)-matmul(cMe,dstran)+ mbinc
                        call  USOLVER(ddsdde_bar,c_dstran,u_dstress,ifstress,ntens)
                        dstran = dstran + c_dstran

                        call  UMAT(stress,statev,ddsdde,sse,spd,scd,                       &
                           rpl,ddsddt,drplde,drpldt,                               &
                           stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
                           ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,  &
                           celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kStep,kinc)

                        if(kiter.lt.maxiter) then                                      ! continue iteration
                           statev(:)=r_statev(:)                                       ! 1) undo the update of state (done by umat)
                           a_dstress  = stress  - r_stress                             ! 2) compute the new approximation of dstress
                           stress(:)=r_stress(:)                                       ! 3) undo the update of (stress done by umat)
                        else
                           stran(:)=stran(:)+dstran(:)                                 !  accept  the updated state and stress (Cartesian)
                        endif
                     endif  ! ==== obey-restrictions

                     if(keywords(2) /= '*ObeyRestrictions'  ) then   ! ======================= disObeyRestrictions ==========
                        u_dstress = 0.0d0
                        where (ifstress == 1)  u_dstress =ddstress -a_dstress           ! undesired Roscoe stress
                        ddsdde_bar = matmul(matmul(M,ddsdde),transpose(M))              ! Roscoe-Roscoe stiffness

                        call  USOLVER(ddsdde_bar,c_dstran,u_dstress,ifstress,ntens)     ! get Rosc. correction  c\_dstran() caused by undesired Rosc. dstress
                        where (ifstress == 1) dstran = dstran + c_dstran                ! corrected Rosc. dstran where stress-controlled
                        dstran_Cart = matmul( transpose(M),dstran )                     ! transsform Rosc. to Cartesian dstran
                        call  UMAT(stress,statev,ddsdde,sse,spd,scd,                    &
                           rpl,ddsddt,drplde,drpldt,                                    &
                           stran,dstran_Cart,time,dtime,temp,dtemp,predef,dpred,cmname, &
                           ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,       &
                           celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kStep,kinc)

                        if (kiter.lt.maxiter) then                                      ! continue iteration
                           statev(:)=r_statev(:)                                        ! 1) forget the changes of state done in umat
                           stress_Rosc = matmul(M,stress)                               !    output from umat transform to Roscoe ?
                           r_stress_Rosc = matmul(M,r_stress)
                           where (ifstress ==1) a_dstress = stress_Rosc - r_stress_Rosc ! 2) compute the new approximation of stress
                           stress(:)=r_stress(:)                                        ! 3) forget the changes of stress done in umat
                        else
                           stran(:)=stran(:)+dstran_Cart(:)                              !  accept  the updated state and stress (Cartesian)
                        endif
                     endif  ! ==== disObey-restrictions


                     if((kiter==maxiter) .and. mod(kinc,10)==0 .and. verbose ) then    ! write to screen after each increment
                        write(*, '(A,I3,A,I3,A,I5,A,I2,A,F9.4,A,F9.4)') &
                           ' ikeyword = ', ikeyword, &
                           ' kstep = ', kStep, &
                           ' kinc = ', kinc, &
                           ' kiter = ', kiter, &
                           ' TEMP = ', TEMP, &
                           ' TIME = ', TIME(1)
                     endif

                  end do do_kiter

                  call check_stress_inc_size(a_dstress, u_dstress)

                  if(keywords(2) =='*DeformationGradient' ) then                    !  rigid rotation of stress
                     T33 = map2T(stress,6)
                     T33 = matmul( matmul(Qb33,T33),transpose(Qb33))
                     stress=map2stress(T33,6)                                          ! rigid rotation of strain
                     eps33 = map2D(stran,6)
                     eps33 = matmul( matmul(Qb33,eps33),transpose(Qb33))
                     stran=map2stran(eps33,6)
                  endif

                  time = set_zero_with_tol(time)
                  stran = set_zero_with_tol(stran)
                  stress = set_zero_with_tol(stress)
                  statev = set_zero_with_tol(statev)

                  if(ievery==1) then
                     write(output_file_id,'(500(g17.10,3h    ))') time+(/dtime,dtime/), stran, stress, statev
                  endif

                  if(keywords(2) =='*PerturbationsS' .or. keywords(2) =='*PerturbationsE' ) then ! having plotted everything undo the increment
                     stran(:)=stran(:) - dstran_Cart(:)
                     statev(:)=r_statev(:)
                     stress(:)=r_stress(:)
                  endif

                  time(1)=time(1)+dtime     !  step time at the beginning of the next increment
                  time(2)=time(2)+dtime     !  total time at the beginning of the next increment
                  Temp = Temp + dTemp       !  AN 2023 total temperature at the beginning of the increment

                  if( existCond ) then                                             ! AN 2016 only if a condition exists
                     if(  EXITNOW(exitCond, stress,stran,statev,nstatv)  ) exit   ! AN 2016 depending on  exitCond go to next step
                  endif
                  ! AN 2016
                  ievery = ievery+1; if(ievery > every) ievery = 1 ! Increment ievery and reset it to 1 if applicable

                  !*****************************************************
                  if(keywords(2) == '*ImportFile' ) then
                     call  tryAlignStress(align, kinc, newState, mImport,stress,ntens)
                  endif
                  !***********************************************

               end do do_kinc
            end do do_istep
         end do do_repet
      end do do_keyword

      close(test_file_id)
      close(import_file_id)
      close(output_file_id)

   end subroutine run_model
end module mod_run_model
