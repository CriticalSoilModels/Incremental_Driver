! Regression tests for the elastic UMAT (src/elastic.f90).
!
! For a linear elastic material with Lamé parameters lambda and mu:
!   lambda = nu*E / ((1+nu)*(1-2*nu))
!   mu     = E / (2*(1+nu))
!
! The stiffness tensor in Voigt notation (3 normal + 3 shear) is diagonal in blocks:
!   DDSDDE(i,i) = lambda + 2*mu  for i = 1,2,3  (normal-normal)
!   DDSDDE(i,j) = lambda          for i/=j, i,j in {1,2,3}
!   DDSDDE(k,k) = mu              for k = 4,5,6  (shear-shear)
!   all other entries = 0
!
! Expected stress increments are derived analytically and compared to UMAT output.
program test_elastic_umat
   use indr_umat, only: UMAT
   implicit none

   integer, parameter :: ntens = 6, ndi = 3, nshr = 3
   integer, parameter :: nstatev = 1, nprops = 2

   ! UMAT argument list
   real(8) :: stress(ntens), statev(nstatev)
   real(8) :: ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens)
   real(8) :: stran(ntens), dstran(ntens)
   real(8) :: time(2), dtime, temp, dtemp, predef(1), dpred(1)
   real(8) :: props(nprops), coords(3), drot(3,3)
   real(8) :: pnewdt, celent, sse, spd, scd, rpl, drpldt
   real(8) :: dfgrd0(3,3), dfgrd1(3,3)
   character(len=80) :: cmname
   integer :: noel, npt, layer, kspt, kstep, kinc

   real(8) :: E, nu, lambda, mu, eps
   real(8), parameter :: tol = 1.0d-10
   integer :: n_fail

   n_fail = 0

   ! --- Material: E = 200 000 kPa, nu = 0.3 ---
   E        = 200000.0d0
   nu       = 0.3d0
   props(1) = E
   props(2) = nu

   lambda = nu * E / ((1.0d0 + nu) * (1.0d0 - 2.0d0*nu))
   mu     = E  / (2.0d0 * (1.0d0 + nu))

   ! --- Fixed bookkeeping arguments (unused by elastic UMAT) ---
   cmname = 'LINEAR_ELASTIC'
   noel=1; npt=1; layer=1; kspt=1; kstep=1; kinc=1
   dtime=1.0d0; temp=0.0d0; dtemp=0.0d0
   time=0.0d0; predef=0.0d0; dpred=0.0d0; coords=0.0d0
   pnewdt=1.0d0; celent=1.0d0
   sse=0.0d0; spd=0.0d0; scd=0.0d0; rpl=0.0d0; drpldt=0.0d0
   drot   = 0.0d0; drot(1,1)=1.0d0;   drot(2,2)=1.0d0;   drot(3,3)=1.0d0
   dfgrd0 = 0.0d0; dfgrd0(1,1)=1.0d0; dfgrd0(2,2)=1.0d0; dfgrd0(3,3)=1.0d0
   dfgrd1 = dfgrd0
   statev = 0.0d0

   eps = 0.001d0

   ! ===========================================================================
   ! Test 1: pure volumetric strain increment
   !   dstran = [eps, eps, eps, 0, 0, 0]
   !   Each diagonal stress component = (3*lambda + 2*mu)*eps
   !   (= 3*K*eps where K = E/(3*(1-2*nu)) is the bulk modulus)
   ! ===========================================================================
   stress = 0.0d0
   stran  = 0.0d0
   dstran = [eps, eps, eps, 0.0d0, 0.0d0, 0.0d0]

   call UMAT(stress, statev, ddsdde, sse, spd, scd, &
             rpl, ddsddt, drplde, drpldt, &
             stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
             ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, &
             celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

   call assert_close(stress(1), (3.0d0*lambda + 2.0d0*mu)*eps, tol, 'vol: sig(1)', n_fail)
   call assert_close(stress(2), (3.0d0*lambda + 2.0d0*mu)*eps, tol, 'vol: sig(2)', n_fail)
   call assert_close(stress(3), (3.0d0*lambda + 2.0d0*mu)*eps, tol, 'vol: sig(3)', n_fail)
   call assert_close(stress(4), 0.0d0, tol, 'vol: sig(4)', n_fail)
   call assert_close(stress(5), 0.0d0, tol, 'vol: sig(5)', n_fail)
   call assert_close(stress(6), 0.0d0, tol, 'vol: sig(6)', n_fail)

   ! ===========================================================================
   ! Test 2: uniaxial strain increment
   !   dstran = [eps, 0, 0, 0, 0, 0]
   !   sig(1) = (lambda + 2*mu)*eps,  sig(2) = sig(3) = lambda*eps
   ! ===========================================================================
   stress = 0.0d0
   stran  = 0.0d0
   dstran = [eps, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]

   call UMAT(stress, statev, ddsdde, sse, spd, scd, &
             rpl, ddsddt, drplde, drpldt, &
             stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
             ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, &
             celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

   call assert_close(stress(1), (lambda + 2.0d0*mu)*eps, tol, 'uniaxial: sig(1)', n_fail)
   call assert_close(stress(2), lambda*eps,               tol, 'uniaxial: sig(2)', n_fail)
   call assert_close(stress(3), lambda*eps,               tol, 'uniaxial: sig(3)', n_fail)
   call assert_close(stress(4), 0.0d0,                   tol, 'uniaxial: sig(4)', n_fail)
   call assert_close(stress(5), 0.0d0,                   tol, 'uniaxial: sig(5)', n_fail)
   call assert_close(stress(6), 0.0d0,                   tol, 'uniaxial: sig(6)', n_fail)

   ! ===========================================================================
   ! Test 3: pure shear strain increment
   !   dstran = [0, 0, 0, eps, 0, 0]  (engineering shear in component 4)
   !   sig(4) = mu*eps,  all normal components = 0
   ! ===========================================================================
   stress = 0.0d0
   stran  = 0.0d0
   dstran = [0.0d0, 0.0d0, 0.0d0, eps, 0.0d0, 0.0d0]

   call UMAT(stress, statev, ddsdde, sse, spd, scd, &
             rpl, ddsddt, drplde, drpldt, &
             stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
             ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, &
             celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)

   call assert_close(stress(1), 0.0d0,   tol, 'shear: sig(1)', n_fail)
   call assert_close(stress(2), 0.0d0,   tol, 'shear: sig(2)', n_fail)
   call assert_close(stress(3), 0.0d0,   tol, 'shear: sig(3)', n_fail)
   call assert_close(stress(4), mu*eps,  tol, 'shear: sig(4)', n_fail)
   call assert_close(stress(5), 0.0d0,   tol, 'shear: sig(5)', n_fail)
   call assert_close(stress(6), 0.0d0,   tol, 'shear: sig(6)', n_fail)

   ! --- Report ---
   if (n_fail == 0) then
      print *, 'PASS  test_elastic_umat'
   else
      print *, 'FAIL  test_elastic_umat:', n_fail, 'assertion(s) failed'
      stop 1
   end if

contains

   subroutine assert_close(got, expected, tol, label, n_fail)
      real(8),      intent(in)    :: got, expected, tol
      character(*), intent(in)    :: label
      integer,      intent(inout) :: n_fail
      ! Relative tolerance scaled to max(|expected|, 1) to handle near-zero gracefully
      if (abs(got - expected) > tol * max(abs(expected), 1.0d0)) then
         print *, 'FAIL  ', trim(label), ': got', got, 'expected', expected
         n_fail = n_fail + 1
      end if
   end subroutine

end program test_elastic_umat
