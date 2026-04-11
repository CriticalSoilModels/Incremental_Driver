module indr_umat_test
   !! Linear isotropic elastic UMAT — test fixture.
   !! Provides a concrete UMAT procedure for unit tests that need one.
   use stdlib_kinds, only: dp
   implicit none
   private
   public :: UMAT

contains

   subroutine UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD,           &
         RPL, DDSDDT, DRPLDE, DRPLDT,                                &
         STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED,    &
         CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS,          &
         COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1,              &
         NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
      implicit none
      character*80  CMNAME
      integer       :: NTENS, NSTATEV, NPROPS, NDI, NSHR, NOEL, &
                       NPT, LAYER, KSPT, KSTEP, KINC
      real(8)       :: SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, DTEMP, &
                       PNEWDT, CELENT
      real(8)       :: STRESS(NTENS), STATEV(NSTATEV),                  &
                       DDSDDE(NTENS,NTENS), DDSDDT(NTENS), DRPLDE(NTENS), &
                       STRAN(NTENS), DSTRAN(NTENS), TIME(2),             &
                       PREDEF(1), DPRED(1), PROPS(NPROPS),               &
                       COORDS(3), DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)

      real(dp) :: E, nu, lam, mu
      integer  :: i, j

      E   = PROPS(1)
      nu  = PROPS(2)
      lam = nu * E / ((1.0_dp + nu) * (1.0_dp - 2.0_dp*nu))
      mu  = E / (2.0_dp * (1.0_dp + nu))

      DDSDDE = 0.0d0
      DDSDDE(1,1) = lam + 2.0_dp*mu;  DDSDDE(2,2) = lam + 2.0_dp*mu
      DDSDDE(3,3) = lam + 2.0_dp*mu
      DDSDDE(1,2) = lam;  DDSDDE(1,3) = lam;  DDSDDE(2,3) = lam
      DDSDDE(2,1) = lam;  DDSDDE(3,1) = lam;  DDSDDE(3,2) = lam
      DDSDDE(4,4) = mu;   DDSDDE(5,5) = mu;   DDSDDE(6,6) = mu

      do i = 1, NTENS
         do j = 1, NTENS
            STRESS(i) = STRESS(i) + DDSDDE(i,j) * DSTRAN(j)
         end do
      end do
   end subroutine UMAT

end module indr_umat_test
