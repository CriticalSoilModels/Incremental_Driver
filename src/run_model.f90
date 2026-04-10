module indr_run_model

   use stdlib_kinds, only: dp
   use stdlib_io, only: open
   use indr_parser, only: splitaLine, ReadStepCommons, PARSER, EXITNOW
   use indr_solver, only: USOLVER

   use indr_types   , only: StressAlignment
   use indr_step_params, only: step_config_t
   use indr_matrices, only: MRoscI, MRoscImt, MRendul, MRendulmT, MRosc, MRoscmT, MCart, MCartmT

   use indr_command_line, only: set_inputs
   use indr_file_io, only: read_parameter_file, read_init_conditions_file, set_output_name_from_test_file, &
      write_line_output_data, write_output_file_header
   use indr_alignment, only: readAlignment, tryAlignStress
   use indr_loads
   use indr_maps
   use indr_constants, only: iter_lower_limit
   use indr_value_checks, only: set_zero_with_tol, check_stress_inc_size

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
      procedure(umat_interface) :: UMAT

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
         initialconditionsfilename, testfilename

      character(len=260) ::  inputline(6), aLine, heading
      character(len=520) :: hugeLine

      character(len=10):: timeHead(2), stranHead(6), stressHead(6)      ! AN 2016
      character(len=15), allocatable :: statevHead(:)                   ! AN 2016

      logical :: verbose
      ! logical :: EXITNOW, existCond,okSplit                             ! AN 2016 ! WaveHello: ExitNow (bool) conflicts with the function
      logical :: okSplit
      real(dp), dimension(6) :: mb   ! constraint RHS for *ObeyRestrictions (not in config)

      integer :: every, ievery
      real(dp),dimension(20) :: oldState, newState, dState
      real(dp), allocatable :: props(:), statev(:), r_statev(:)

      real(dp),dimension(3,3):: Qb33,eps33,T33

      integer :: maxiter, kiter, ikeyword, &
         iRepetition, nRepetitions, kStep,iStep,nSteps,ntens_in

      real(dp):: r_stress(ntens),a_dstress(ntens),u_dstress(ntens),&
         stress_Rosc(ntens),r_stress_Rosc(ntens),           &
         ddstress(ntens), c_dstran(ntens),                  &
         dstran_Cart(6), ddsdde_bar(6,6)

      real(dp),dimension(1:6,1:6)::M,MmT      !  currrent $\cM$ and $\cM^{-T}$  for a given iStep

      type(StressAlignment) :: align
      type(step_config_t) :: ofStep(30)   !  stores descriptions of up to 30 steps which are repeated
      type(step_config_t) :: config       !  current step configuration, passed to get_increment


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

               if(iRepetition > 1) then  ! restore step config saved on the first iteration
                  config      = ofStep(iStep)
                  keywords(2) = config%load_type
                  keywords(3) = config%coord_sys
                  maxiter     = config%max_iter
               endif

               if(keywords(1) == '*Repetition') read(test_file_id,'(a)') keywords(2) ! = *LinearLoad  or *CirculatingLoad or *ObeyRestrictions...
               ! otherwise keywords(2) = keywords(1)

               call splitaLine(keywords(2),'?', keywords(2), config%exit_cond, config%has_exit_cond)

               keywords(2)  = trim(keywords(2)) ! Trim the test name

               ! Defaults — read routines only set the fields they use.
               config%load_type       = keywords(2)
               config%ifstress        = 0
               config%delta_load      = 0.0_dp
               config%delta_load_circ = 0.0_dp
               config%phase0          = 0.0_dp
               config%dfgrd0          = delta
               config%dfgrd1          = delta
               config%cMt             = 0.0_dp
               config%cMe             = 0.0_dp
               config%mbinc           = 0.0_dp
               config%columns_in_file = 0
               config%import_factor   = 1.0_dp
               config%n_import        = 0

               if(keywords(2) == '*DeformationGradient') then
                  call read_deformation_gradient_load(test_file_id, config, every)

               else if (keywords(2) == '*CirculatingLoad') then
                  call read_circulating_load(test_file_id, config, every)

               else if(keywords(2) == '*LinearLoad') then
                  call read_linear_load(test_file_id, config, every)

               else if(keywords(2)(1:11) == '*ImportFile') then
                  call read_file_load(test_file_id, config, every, align)

               else if(keywords(2) == '*OedometricE1') then
                  call read_oedometric_load(test_file_id, config, every)

               else if(keywords(2) == '*OedometricS1') then
                  call read_oedometric_S1_load(test_file_id, config, every)

               else if(keywords(2) == '*TriaxialE1') then
                  call read_triaxial_e1_load(test_file_id, config, every)

               else if(keywords(2) == '*TriaxialS1') then
                  call read_triaxial_s1_load(test_file_id, config, every)

               else if(keywords(2) == '*TriaxialUEq') then
                  call read_triaxial_ueq_load(test_file_id, config, every)

               else if(keywords(2) == '*TriaxialUq') then
                  call read_triaxial_uq_load(test_file_id, config, every)

               else if(keywords(2) == '*PureRelaxation') then
                  call read_pure_relaxation_load(test_file_id, config, every)

               else if(keywords(2) == '*PureCreep') then
                  call read_pure_creep_load(test_file_id, config, every)

               else if(keywords(2) == '*UndrainedCreep') then
                  call read_undrained_creep(test_file_id, config, every)

               else if(keywords(2) == '*ObeyRestrictions') then
                  call read_obey_restrictions_load(test_file_id, config, every, mb)

               else if(keywords(2) == '*PerturbationsS') then
                  call read_perturbations_S_load(test_file_id, config, every)

               else if(keywords(2) == '*PerturbationsE') then
                  call read_perturbations_E_load(test_file_id, config, every)

               else if(keywords(2) == '*RandomWalk') then
                  call read_random_walk_load(test_file_id, config, every)

               else if(keywords(2) == '*End') then
                  print *, '*End encountered in test.inp'
                  return
               else
                  write(*,*) 'error: unknown keywords(2)=',keywords(2)
                  error stop 'stopped by unknown keyword(2) in test.inp'
               end if

               ! Sync keywords(3) from config (read routines set coord_sys)
               keywords(3) = trim(config%coord_sys)

               if(keywords(1) == '*Repetition' .and. iRepetition == 1) then
                  ofStep(iStep) = config
               endif

               maxiter = config%max_iter

               ! If the choosen load is stress controlled make sure at least iter_lower_limit number of iterations is done
               ! A guess at the correct strain has to be made and the stress has to be converged to
               if(any(config%ifstress==1)) maxiter = max(maxiter,iter_lower_limit)
               if(all(config%ifstress==0) .and. keywords(2) .ne. '*ObeyRestrictions') maxiter = 1  ! no iterations are necessary

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
                  import_file_id = open(config%import_file)
                  do                                                             ! AN 2016
                     read(import_file_id,'(a)',iostat = iostat) hugeLine;         ! AN 2016

                     if (iostat /= 0) error stop 'Error reading ImportFile in the first non-numeric records '
                     hugeLine= adjustL(hugeLine) ; aChar = hugeLine(1:1)         ! AN 2016
                     if(index('1234567890+-.',aChar) > 0) exit                  ! preceding non-numeric lines in ImportFile will be ignored
                  enddo

                  read(hugeLine,*,iostat=iostat) oldState(1:config%n_import)
                  if (iostat /= 0) error stop 'Error reading ImportFile in the first numeric record '
               endif

               ievery=1
               do_kinc: do kinc=1,config%n_inc

                  if(keywords(2) == '*ImportFile') then                             ! AN 2016
                     dTemp = 0                                                      ! AN 2023 temperat
                     read(import_file_id,*,iostat=iostat)  newState(1:config%n_import)                 ! AN 2016
                     if(iostat > 0) then                                                  ! AN 2016
                        write(*,*) 'error Import file', config%import_file, 'line=', kinc+1
                        error stop                                                            ! AN 2016
                     else                                                                                                         ! AN 2016
                        close(import_file_id)                                       ! AN 2016
                        write(*,*) 'finished reading file', config%import_file      ! AN 2016
                        exit                                                        ! AN 2016
                     endif                                                           ! AN 2016

                     dState = newState(:) -  oldState(:)                             ! AN 2016
                     do i=1,6                                                        ! AN 2016
                        if  (config%columns_in_file(i) == 0) cycle                      ! AN 2016
                        if (config%ifstress(i)==1) ddstress(i)= dState(config%columns_in_file(i))* config%import_factor(i)
                        if (config%ifstress(i)==0) dstran(i)  = dState(config%columns_in_file(i))* config%import_factor(i)
                     enddo

                     if(config%columns_in_file(7)/= 0) dtime = dState(config%columns_in_file(7)) * config%import_factor(7)   ! AN 2016
                     oldState(:) = newState(:)                                    ! AN 2016
                  endif                                                             ! AN 2016

                  if(keywords(2) /= '*ImportFile') then
                     call get_increment(config, time, dtime, ddstress, dstran, dTemp, Qb33, dfgrd0, dfgrd1, drot)
                  endif


                  a_dstress(:)= 0.0d0    ! approximated Roscoe's dstress
                  r_statev(:)=statev(:)  ! remember the initial state and stress till the iteration is completed
                  r_stress(:)= stress    ! remembered Cartesian stress

                  do_kiter: do kiter=1, maxiter  !--------Equilibrium Iteration--------
                     c_dstran(:) = 0

                     if(keywords(2)== '*ObeyRestrictions'  ) then  ! ======================= ObeyRestrictions =========
                        ddsdde_bar = matmul(config%cMt,ddsdde) + config%cMe
                        u_dstress = - matmul(config%cMt,a_dstress)-matmul(config%cMe,dstran)+ config%mbinc
                        call  USOLVER(ddsdde_bar,c_dstran,u_dstress,config%ifstress,ntens)
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
                        where (config%ifstress == 1)  u_dstress =ddstress -a_dstress           ! undesired Roscoe stress
                        ddsdde_bar = matmul(matmul(M,ddsdde),transpose(M))              ! Roscoe-Roscoe stiffness

                        call  USOLVER(ddsdde_bar,c_dstran,u_dstress,config%ifstress,ntens)     ! get Rosc. correction  c\_dstran() caused by undesired Rosc. dstress
                        where (config%ifstress == 1) dstran = dstran + c_dstran                ! corrected Rosc. dstran where stress-controlled
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
                           where (config%ifstress ==1) a_dstress = stress_Rosc - r_stress_Rosc ! 2) compute the new approximation of stress
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

                  if( config%has_exit_cond ) then                                  ! AN 2016 only if a condition exists
                     if(  EXITNOW(config%exit_cond, stress,stran,statev,nstatv)  ) exit   ! AN 2016 depending on  exitCond go to next step
                  endif
                  ! AN 2016
                  ievery = ievery+1; if(ievery > every) ievery = 1 ! Increment ievery and reset it to 1 if applicable

                  !*****************************************************
                  if(keywords(2) == '*ImportFile' ) then
                     call  tryAlignStress(align, kinc, newState, config%n_import,stress,ntens)
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
end module indr_run_model
