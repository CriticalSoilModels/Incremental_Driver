! Functions that are used by incrmental driver


module mod_inc_driver_funcs
   use stdlib_kinds, only: dp
   use indr_linalg, only: inv33, spectral_decomposition_of_symmetric, &
                          app_jacobian_similarity, get_jacobian_rot
   use indr_abaqus_utils, only: SINV, ROTSIG, SPRINC, SPRIND, XIT
   use indr_parser, only: splitaLine, ReadStepCommons, PARSER, EXITNOW
   use indr_solver, only: USOLVER
   implicit none
   private
   public :: splitaLine, ReadStepCommons, PARSER, get_increment, USOLVER, EXITNOW, &
             SINV, ROTSIG, SPRINC, SPRIND, &
             inv33, spectral_decomposition_of_symmetric
contains


! ==========================================================================
!     basing on an input command with parameters converts  deltaLoad or deltaLoadCirc
!     to the canonical three lists:  dstress(), dstrain(), ifstress()
!     get\_increment is called in each increment (and not once per step )
   subroutine get_increment(keywords, time, deltaTime,ifstress,ninc,      &
      deltaLoadCirc,phase0,deltaLoad,deltaTemp, &  ! AN 2023 temperat
      dtime, ddstress,  dstran , dTemp, Qb33,   &! AN 2023 temperat
      dfgrd0, dfgrd1,drot )
      implicit none
      character(40):: keywords(10)
      integer, intent(in)  :: ifstress(6),ninc
      real(dp), intent(in) :: time(2), deltaTime, deltaLoadCirc(6),phase0(6), deltaLoad(9), deltaTemp
      real(dp), intent(out) ::  dtime, ddstress(6), dstran(6), Qb33(3,3), dTemp
      real(dp), intent(in out) ::  dfgrd0(3,3), dfgrd1(3,3), drot(3,3)


      real(dp), parameter :: Pi = 3.1415926535897932385d0
      real(dp),parameter,dimension(3,3):: delta = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
      real(dp),dimension(3,3):: Fb,Fbb, dFb,aux33,dLb,depsb,dOmegab
      real(dp):: wd(6),  & ! angular velocity (in future individual for each component)
         w0(6),  & ! initial phase shift for a component
         t         ! step time
      real(dp) :: arandom
      integer(4) :: i
      logical :: ok

      dtime =  deltaTime/ ninc
      dTemp = deltaTemp / ninc    ! AN 2023 temperat   (perturbations and random walk may need some programming)
      dstran= 0
      ddstress=0
      Qb33 = delta
      drot = delta
      dfgrd0=delta
      dfgrd1=delta

      !------------------------------------------------------
      if(keywords(2) == '*LinearLoad') then                               !  proportional loading
         do i=1,6
            if (ifstress(i)==1)   ddstress(i) = deltaLoad(i)/ ninc
            if (ifstress(i)==0)    dstran(i) = deltaLoad(i)/ ninc              ! log strain -> corresp. displac. inc. not constant
         enddo
         ! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
         ! for dfgrd0 use stran
         ! for dfgrd1 use stran-dstran
      endif
      !--------------------------------------------------
      if(keywords(2) == '*DeformationGradient') then                     ! full deformation gradient.
         ! finite rotations calculated after Hughes+Winget 1980
         Fb = reshape((/deltaLoad(1), deltaLoad(5), deltaLoad(7),  &
            deltaLoad(4), deltaLoad(2), deltaLoad(9),    &
            deltaLoad(6), deltaLoad(8), deltaLoad(3)/),  &
            (/3,3/))

         Fbb = delta + (Fb-delta)*(time(1)/deltaTime)
         dfgrd0  = Fbb
         dFb = (Fb-delta)/ninc
         aux33 =  Fbb + dFb/2.0d0
         dfgrd1   = Fbb  + dFb

         !  call matrix('inverse', aux33, 3, ok )
         aux33 = inv33(aux33)
         dLb =  matmul(dFb,aux33)
         depsb = 0.5d0*(dLb + transpose(dLb))
         dstran=(/depsb(1,1), depsb(2,2),depsb(3,3), 2.0d0*depsb(1,2),2.0d0*depsb(1,3),2.0d0*depsb(2,3)/)
         dOmegab =    0.5d0*(dLb - transpose(dLb))
         aux33 =  delta - 0.5d0*dOmegab
         !     call matrix('inverse', aux33, 3, ok )
         aux33 = inv33(aux33)
         Qb33 = matmul(aux33, (delta+0.5d0*dOmegab))
         drot=Qb33
      endif
      !------------------------------------------------------
      if(keywords(2) == '*CirculatingLoad' )then                         !  harmonic oscillation
         wd(:) = 2*Pi/deltaTime
         w0 = phase0
         t= time(1)  + dtime/2   ! step time in the middle of the increment
         do i=1,6
            if(ifstress(i)==1) ddstress(i)= dtime * deltaLoadCirc(i) * wd(i) * Cos(wd(i) * t + w0(i)) + deltaLoad(i)/ ninc
            if(ifstress(i)==0) dstran(i)  = dtime * deltaLoadCirc(i) * wd(i) * Cos(wd(i) * t + w0(i)) + deltaLoad(i)/ ninc
         enddo
         ! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
         ! for dfgrd0 use stran
         ! for dfgrd1 use stran-dstran
      endif

      !--------------------------------------------------------
      if(keywords(2) == '*PerturbationsS' )then
         ddstress(1)= deltaLoad(1)*cos( time(1)*2*Pi/deltaTime )
         ddstress(2)= deltaLoad(1)*sin( time(1)*2*Pi/deltaTime  )
         ! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
         ! for dfgrd0 use stran
         ! for dfgrd1 use stran-dstran
      endif

      !--------------------------------------------------------
      if(keywords(2) == '*PerturbationsE' )then
         dstran(1)= deltaLoad(1)*cos( time(1)*2*Pi/deltaTime )
         dstran(2)= deltaLoad(1)*sin( time(1)*2*Pi/deltaTime  )
         ! here dfgrd0 and dfgrd1   can be defined from stran assuming polar decomposition F=V.R with R=1  and V = exp(stran)
         ! for dfgrd0 use stran
         ! for dfgrd1 use stran-dstran
      endif

      if(keywords(2) == '*RandomWalk' )then
         call random_seed
         do i =1,6
            call random_number(arandom)
            if( ifstress(i)== 1) ddstress(i)= 2*(arandom-0.5d0)*deltaLoad(i)
            if( ifstress(i)== 0) dstran(i)= 2*(arandom-0.5d0)*deltaLoad(i)
         enddo
      endif

      return

   end subroutine get_increment

   ! splitaLine, ReadStepCommons, PARSER, EXITNOW moved to indr_parser (re-exported via use above)
   ! USOLVER moved to indr_solver (re-exported via use above)



   subroutine stopp(i, whyStopText)                      ! AN 2016
      USE ISO_FORTRAN_ENV  ! , ONLY : ERROR\_UNIT           ! AN 2016
      implicit none                                         ! AN 2016
      integer, intent(in) :: i                             ! AN 2016
      character(*)      :: whyStopText                    ! AN 2016
      stop  'whyStopText'                                   ! AN 2016
      WRITE(ERROR_UNIT,*)   whyStopText                  ! AN 2016
      CALL EXIT(5)                                       ! AN 2016
   end subroutine stopp                             ! AN 2016

   ! subroutine get_umat_stiffness_matrix
   !    !! Calls the users umat to get the stiffness matrix.
   !    !! Sends dummy values or resets values for variables in case the variables are updated when they shouldn't be
       

   !    dstran(:)=0
   !    dtime=0
   !    dtemp=0
   !    kinc=0
   !    r_statev(:)=statev(:);  r_stress(:)=stress(:)     ! AN 21.06.2017 remember the initial state and stress
   !    !=== first call umat with dstrain=0 dtime=0 just for stiffness (=jacobian ddsdde)
   !    call  UMAT(stress,statev,ddsdde,sse,spd,scd,                       &
   !       rpl,ddsddt,drplde,drpldt,                               &
   !       stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
   !       ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,  &
   !       celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,0,kinc)   !=== some constitutive models require kStep=0 other do not

   !    statev(:)=r_statev(:);  stress(:)=r_stress(:)   !  AN 21.06.2017 recover stress and state  although the ZERO call of umat should not modify them

   ! end subroutine get_umat_stiffness_matrix
end module mod_inc_driver_funcs
