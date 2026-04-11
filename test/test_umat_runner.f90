program test_umat_runner
   !! Tests for umat_runner_t.
   !! Uses the linear elastic UMAT from indr_umat and verifies that
   !! run_umat produces the same stress increment as calling the UMAT directly.
   use stdlib_kinds, only: dp
   use indr_types, only: material_state_t
   use indr_model_runner, only: model_runner_t
   use indr_umat_runner, only: umat_runner_t
   use indr_umat_test, only: UMAT
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-10_dp

   call test_elastic_stress_increment()
   call test_stiffness_returned()

   if (nfail == 0) then
      print *, 'PASS  test_umat_runner'
   else
      print *, 'FAIL  test_umat_runner — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   subroutine test_elastic_stress_increment()
      !! A uniaxial strain increment deps(1)=0.001 on a zero-stress state
      !! should produce a stress increment consistent with the elastic stiffness.
      integer, parameter :: ntens = 6, ndi = 3, nshr = 3
      real(dp), parameter :: E = 10000.0_dp, nu = 0.3_dp
      real(dp), parameter :: deps_ref(ntens) = [0.001_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
      real(dp), parameter :: identity33(3,3) = reshape( &
         [1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])

      type(umat_runner_t) :: runner
      type(material_state_t) :: state
      real(dp) :: ddsig_by_ddeps(ntens, ntens)

      ! Set up runner
      runner%proc  => UMAT
      allocate(runner%props(2), source=[E, nu])
      runner%nstatv = 1
      runner%cmname = 'ELASTIC'

      ! Set up state
      state%sig    = 0.0_dp
      state%eps    = 0.0_dp
      state%time   = 0.0_dp
      state%dt     = 1.0_dp
      state%temp   = 0.0_dp
      state%F_start = identity33
      state%F_end   = identity33
      allocate(state%statev(1), source=0.0_dp)

      call runner%run(state, deps_ref, ddsig_by_ddeps, ndi, nshr, ntens, 0.0_dp)

      ! For isotropic linear elastic: sig(1) = (lam + 2*mu) * deps(1)
      ! lam = nu*E / ((1+nu)*(1-2*nu)),  mu = E / (2*(1+nu))
      block
         real(dp) :: lam, mu, expected_sig1
         lam = nu * E / ((1.0_dp + nu) * (1.0_dp - 2.0_dp*nu))
         mu  = E / (2.0_dp * (1.0_dp + nu))
         expected_sig1 = (lam + 2.0_dp*mu) * deps_ref(1)

         if (abs(state%sig(1) - expected_sig1) > tol * abs(expected_sig1)) then
            print *, 'FAIL  test_elastic_stress_increment: sig(1) =', state%sig(1), &
               ' expected', expected_sig1
            nfail = nfail + 1
         end if
         ! Lateral stresses should be lam * deps(1)
         if (abs(state%sig(2) - lam * deps_ref(1)) > tol * abs(lam * deps_ref(1))) then
            print *, 'FAIL  test_elastic_stress_increment: sig(2) =', state%sig(2), &
               ' expected', lam * deps_ref(1)
            nfail = nfail + 1
         end if
      end block
   end subroutine test_elastic_stress_increment

   subroutine test_stiffness_returned()
      !! ddsig_by_ddeps(1,1) should equal lam + 2*mu for isotropic linear elastic.
      integer, parameter :: ntens = 6, ndi = 3, nshr = 3
      real(dp), parameter :: E = 10000.0_dp, nu = 0.3_dp
      real(dp), parameter :: deps_zero(ntens) = 0.0_dp
      real(dp), parameter :: identity33(3,3) = reshape( &
         [1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])

      type(umat_runner_t) :: runner
      type(material_state_t) :: state
      real(dp) :: ddsig_by_ddeps(ntens, ntens)
      real(dp) :: lam, mu, expected_D11, expected_D44

      runner%proc  => UMAT
      allocate(runner%props(2), source=[E, nu])
      runner%nstatv = 1
      runner%cmname = 'ELASTIC'

      state%sig    = 0.0_dp
      state%eps    = 0.0_dp
      state%time   = 0.0_dp
      state%dt     = 0.0_dp
      state%temp   = 0.0_dp
      state%F_start = identity33
      state%F_end   = identity33
      allocate(state%statev(1), source=0.0_dp)

      call runner%run(state, deps_zero, ddsig_by_ddeps, ndi, nshr, ntens, 0.0_dp)

      lam = nu * E / ((1.0_dp + nu) * (1.0_dp - 2.0_dp*nu))
      mu  = E / (2.0_dp * (1.0_dp + nu))
      expected_D11 = lam + 2.0_dp*mu
      expected_D44 = mu

      if (abs(ddsig_by_ddeps(1,1) - expected_D11) > tol * expected_D11) then
         print *, 'FAIL  test_stiffness_returned: D(1,1) =', ddsig_by_ddeps(1,1), &
            ' expected', expected_D11
         nfail = nfail + 1
      end if
      if (abs(ddsig_by_ddeps(4,4) - expected_D44) > tol * expected_D44) then
         print *, 'FAIL  test_stiffness_returned: D(4,4) =', ddsig_by_ddeps(4,4), &
            ' expected', expected_D44
         nfail = nfail + 1
      end if
   end subroutine test_stiffness_returned

end program test_umat_runner
