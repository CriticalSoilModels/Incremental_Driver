! Functions that are used by incrmental driver


module mod_inc_driver_funcs
   use stdlib_kinds, only: dp
   use indr_linalg, only: inv33, spectral_decomposition_of_symmetric, &
                          app_jacobian_similarity, get_jacobian_rot
   use indr_abaqus_utils, only: SINV, ROTSIG, SPRINC, SPRIND, XIT
   use indr_parser, only: splitaLine, ReadStepCommons, PARSER, EXITNOW
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

   subroutine USOLVER(KK,u,rhs,is,ntens) ! 23.7.2008  new usolver with improvement after numerical recipes
      !  solver for unsymmetric matrix and  unknowns on both sides of equation
      !  \com    KK - stiffness  is not spoiled  within the subroutine
      !  \com    u - strain   rhs - stress
      !  \com    is(i)= 1 means rhs(i) is prescribed,
      !  \com    is(i)= 0  means u(i) is prescribed

      implicit none
      integer, intent(in)                             :: ntens
      integer, dimension(1:ntens), intent(in)         :: is
      real(dp), dimension(1:ntens,1:ntens), intent(in) :: KK
      real(dp), dimension(1:ntens), intent(inout)      ::  u,rhs
      real(dp), dimension(1:ntens)                     :: rhs1
      real(dp), allocatable                            :: rhsPrim(:), KKprim(:,:), uprim(:)
      integer ::  i,j,ii,nis
      integer,allocatable :: is1(:)

      nis = sum(is)                                                     ! \com number of prescribed stress components

      if (all( is(1:ntens)== 0) ) then
         rhs =  matmul(KK,u)
         return
      endif

      if (all(is(1:ntens) == 1)) then                                   ! \com a special case with full stress control
         u =xLittleUnsymmetricSolver(KK,rhs)
         return
      endif

      rhs1 = rhs  ! \com modify the rhs  to rhs1
      do i=1,ntens
         if (is(i) == 0) rhs1 = rhs1 - u(i)*KK(:,i)                        ! \com  modify rhs wherever strain control
      enddo

      allocate(KKprim(nis,nis), rhsprim(nis), uprim(nis), is1(nis))     ! \com re-dimension  stiffness and rhs

      ii=0
      do i=1,ntens
         if(is(i)==1) then
            ii = ii+1
            is1(ii) = i                                                    ! list with positions  of is(i) == 1
         endif
      enddo


      do i=1,nis
         rhsPrim(i) = rhs1( is1(i) )
         do j=1,nis
            KKprim(i,j) =  KK(is1(i),is1(j))
         enddo
      enddo

      if (nis ==1) uprim = rhsprim / KKprim(1,1)
      if (nis > 1) uprim =xLittleUnsymmetricSolver(KKprim,rhsprim)
      do i=1,nis
         u(is1(i)) = uprim(i)
      enddo
      do i=1,ntens
         if ( is(i) == 0 ) rhs(i) = dot_product( KK(i,:), u)             ! \com   calculate rhs where u prescribed
      enddo
      deallocate(KKprim,rhsprim,uprim,is1)


   CONTAINS  !===================================================

!  contained in USOLVER LU-decomposition from NR
      SUBROUTINE ludcmp(a,indx,d)
         IMPLICIT NONE
         real(dp), DIMENSION(:,:), INTENT(INOUT) :: a
         INTEGER, DIMENSION(:), INTENT(OUT) :: indx
         real(dp), INTENT(OUT) :: d
         real(dp), DIMENSION(size(a,1)) :: vv ,aux
         integer, dimension(1) :: imaxlocs
         real(dp), PARAMETER :: TINY=1.0d-20
         INTEGER  :: j,n,imax
         n = size(a,1)
         d=1.0
         vv=maxval(abs(a),dim=2)
         if (any(vv == 0.0)) stop 'singular matrix in ludcmp'
         vv=1.0d0/vv
         do j=1,n
            imaxlocs=maxloc(  vv(j:n)*abs( a(j:n,j) ) )
            imax=(j-1)+imaxlocs(1)
            if (j /= imax) then
               aux = a(j,:)      ! call swap(a(imax,:),a(j,:))
               a(j,:) = a(imax,:)
               a(imax,:) = aux
               d=-d
               vv(imax)=vv(j)
            end if
            indx(j)=imax
            if (a(j,j) == 0.0) a(j,j)=TINY
            a(j+1:n,j)=a(j+1:n,j)/a(j,j)
            a(j+1:n,j+1:n)=a(j+1:n,j+1:n)- spread(a(j+1:n,j),2,n-j )* spread(a(j,j+1:n),1, n-j)   ! outerprod
         end do
      END SUBROUTINE ludcmp

!  contained in USOLVER LU-back substitution from NR
      SUBROUTINE lubksb(a,indx,b)
         IMPLICIT NONE
         real(dp), DIMENSION(:,:), INTENT(IN) :: a
         INTEGER, DIMENSION(:), INTENT(IN) :: indx
         real(dp), DIMENSION(:), INTENT(INOUT) :: b
         INTEGER :: i,n,ii,ll
         real(dp) :: summ
         n=size(a,1)
         ii=0
         do i=1,n
            ll=indx(i)
            summ=b(ll)
            b(ll)=b(i)
            if (ii /= 0) then
               summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
            else if (summ /= 0.0) then
               ii=i
            end if
            b(i)=summ
         end do
         do i=n,1,-1
            b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
         end do
      END SUBROUTINE lubksb

!  contained in USOLVER  improvement of the accuracy
      SUBROUTINE mprove(a,alud,indx,b,x)
         IMPLICIT NONE
         real(dp), DIMENSION(:,:), INTENT(IN) :: a,alud
         INTEGER, DIMENSION(:), INTENT(IN) :: indx
         real(dp), DIMENSION(:), INTENT(IN) :: b
         real(dp), DIMENSION(:), INTENT(INOUT) :: x
         real(dp), DIMENSION(size(a,1)) :: r
         r=matmul(a,x)-b
         call lubksb(alud,indx,r)
         x=x-r
      END SUBROUTINE mprove

!  solver contained in USOLVER  for problems with unknowns on the left-hand side
      function xLittleUnsymmetricSolver(a,b)
         IMPLICIT NONE                !==== solves $a . x = b$ \& doesn't spoil a or b
         real(dp), DIMENSION(:), intent(inout) :: b
         real(dp), DIMENSION(:,:), intent(in) ::  a
         real(dp), DIMENSION(size(b,1)) :: x
         real(dp), DIMENSION(size(b,1),size(b,1)) ::  aa
         INTEGER, DIMENSION(1:size(b,1)) :: indx
         real(dp), DIMENSION(1:size(b,1)):: xLittleUnsymmetricSolver
         real(dp) :: d
         x(:)=b(:)
         aa(:,:)=a(:,:)
         call ludcmp(aa,indx,d)
         call lubksb(aa,indx,x)
         call mprove(a,aa,indx,b,x)
         xLittleUnsymmetricSolver= x(:)
      end function xLittleUnsymmetricSolver


   end subroutine USOLVER



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
