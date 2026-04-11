!! Core computation layer — file-free step integration.
!! integrate_step runs the full do_kinc loop for one load step, delegating
!! each model call to the provided runner. No output file I/O occurs here.
!! (*ImportFile load-data reading is I/O internal to the load definition,
!! not output writing, and is handled inside the do_kinc loop.)
module indr_step_runner
   use stdlib_kinds, only: dp
   use stdlib_io,    only: open
   use indr_step_params,   only: step_config_t, material_state_t
   use indr_model_runner,  only: model_runner_t
   use indr_types,         only: StressAlignment
   use indr_loads,         only: get_increment
   use indr_solver,        only: USOLVER
   use indr_parser,        only: EXITNOW
   use indr_alignment,     only: tryAlignStress
   use indr_maps,          only: map2T, map2stress, map2D, map2stran
   use indr_matrices,      only: MRoscI, MRoscImT, MRendul, MRendulmT, MRosc, MRoscmT, MCart, MCartmT
   use indr_value_checks,  only: set_zero_with_tol, check_stress_inc_size
   use indr_constants,     only: iter_lower_limit
   implicit none(type, external)
   private
   public :: integrate_step

contains

   subroutine integrate_step(config, state, runner, results, align)
      !! Integrate one load step with no output file I/O.
      !!
      !! config  -- step parameters (load type, n_inc, ifstress, coord_sys, …)
      !! state   -- model state on entry; updated in-place to end-of-step values
      !! runner  -- concrete model runner (umat_runner_t or mcss_runner_t)
      !! results -- one material_state_t snapshot per increment (allocated here)
      !! align   -- optional stress alignment for *ImportFile steps
      type(step_config_t),                 intent(in)    :: config
      type(material_state_t),              intent(inout) :: state
      class(model_runner_t),               intent(inout) :: runner
      type(material_state_t), allocatable, intent(out)   :: results(:)
      type(StressAlignment),  optional,    intent(in)    :: align

      integer, parameter :: ntens = 6, ndi = 3, nshr = 3

      ! Coordinate transform matrices for the chosen coordinate system
      real(dp) :: M(ntens, ntens), M_inv_T(ntens, ntens)

      ! Per-increment load and solver quantities
      real(dp) :: deps(ntens)             ! strain increment (Roscoe coords, or Cartesian for ObeyRestr.)
      real(dp) :: dsig(ntens)             ! prescribed stress increment
      real(dp) :: deps_Cart(ntens)        ! strain increment in Cartesian coords
      real(dp) :: deps_corr(ntens)        ! USOLVER correction to strain
      real(dp) :: dsig_approx(ntens)      ! running Newton approximation of actual stress increment
      real(dp) :: dsig_target(ntens)      ! residual stress that USOLVER must correct
      real(dp) :: ddsig_by_ddeps(ntens, ntens)      ! ∂σ/∂ε returned by the model
      real(dp) :: ddsig_by_ddeps_bar(ntens, ntens)  ! transformed tangent for USOLVER

      ! Rollback buffers (restored between Newton iterations)
      real(dp) :: sig_saved(ntens)
      real(dp), allocatable :: statev_saved(:)

      ! Roscoe-space stress for non-ObeyRestrictions Newton update
      real(dp) :: sig_Rosc(ntens), sig_Rosc_saved(ntens)

      ! Rotation from polar decomposition (DeformationGradient steps)
      real(dp) :: R_polar(3,3), drot(3,3)
      real(dp) :: T33(3,3), eps33(3,3)

      ! Time and temperature for this increment
      real(dp) :: dt, dtemp

      ! *ImportFile state
      integer  :: import_file_id, iostat
      character(len=1)   :: aChar
      character(len=520) :: hugeLine
      real(dp) :: oldState(20), newState(20), dState(20)

      ! Results buffer — pre-allocated to n_inc, trimmed after the loop
      type(material_state_t), allocatable :: results_buf(:)
      integer :: n_results

      integer :: maxiter, kiter, kinc, i
      logical :: is_obey

      ! --- Coordinate transform matrix ---
      select case (trim(config%coord_sys))
       case ('*Cartesian');        M = MCart;  M_inv_T = MCartmT
       case ('*Roscoe');           M = MRosc;  M_inv_T = MRoscmT
       case ('*RoscoeIsomorph');   M = MRoscI; M_inv_T = MRoscImT
       case ('*Rendulic');         M = MRendul; M_inv_T = MRendulmT
       case default
         write(*,*) 'integrate_step: unknown coord_sys:', trim(config%coord_sys)
         error stop 'integrate_step: unknown coord_sys'
      end select

      is_obey = (trim(config%load_type) == '*ObeyRestrictions')

      ! --- maxiter ---
      maxiter = config%max_iter
      if (any(config%ifstress == 1))                          maxiter = max(maxiter, iter_lower_limit)
      if (all(config%ifstress == 0) .and. .not. is_obey)     maxiter = 1

      ! --- Allocate rollback buffer for statev ---
      allocate(statev_saved(size(state%statev)))

      ! --- Zero-deps call for initial stiffness ---
      ! Calls the model with zero strain/time increment to get the initial
      ! tangent ddsig_by_ddeps. Stress and statev are rolled back immediately.
      sig_saved    = state%sig
      statev_saved = state%statev
      dt           = state%dt
      state%dt     = 0.0_dp
      call runner%run(state, [(0.0_dp, i=1,ntens)], ddsig_by_ddeps, ndi, nshr, ntens)
      state%sig    = sig_saved
      state%statev = statev_saved
      state%dt     = dt

      ! --- Open *ImportFile if needed ---
      if (trim(config%load_type) == '*ImportFile') then
         import_file_id = open(config%import_file)
         ! Skip leading non-numeric header lines
         do
            read(import_file_id, '(a)', iostat=iostat) hugeLine
            if (iostat /= 0) error stop 'integrate_step: error reading ImportFile headers'
            hugeLine = adjustl(hugeLine)
            aChar    = hugeLine(1:1)
            if (index('1234567890+-.', aChar) > 0) exit
         end do
         read(hugeLine, *, iostat=iostat) oldState(1:config%n_import)
         if (iostat /= 0) error stop 'integrate_step: error reading first ImportFile record'
      end if

      ! --- Allocate results buffer ---
      allocate(results_buf(config%n_inc))
      n_results = 0

      ! ======================================================================
      do_kinc: do kinc = 1, config%n_inc

         ! --- Get increment load ---
         if (trim(config%load_type) == '*ImportFile') then
            dtemp = 0.0_dp
            read(import_file_id, *, iostat=iostat) newState(1:config%n_import)
            if (iostat > 0) then
               write(*,*) 'integrate_step: error reading ImportFile, line', kinc + 1
               error stop
            else if (iostat < 0) then
               ! EOF — finished reading file
               close(import_file_id)
               write(*,*) 'integrate_step: finished reading', trim(config%import_file)
               exit do_kinc
            end if
            dState = newState - oldState
            dsig   = 0.0_dp
            deps   = 0.0_dp
            do i = 1, 6
               if (config%columns_in_file(i) == 0) cycle
               if (config%ifstress(i) == 1) dsig(i)  = dState(config%columns_in_file(i)) * config%import_factor(i)
               if (config%ifstress(i) == 0) deps(i)  = dState(config%columns_in_file(i)) * config%import_factor(i)
            end do
            if (config%columns_in_file(7) /= 0) then
               dt = dState(config%columns_in_file(7)) * config%import_factor(7)
            else
               dt = config%delta_time / config%n_inc
            end if
            oldState = newState
         else
            call get_increment(config, state%time, dt, dsig, deps, dtemp, R_polar, state%F_start, state%F_end, drot)
         end if

         ! --- Set up for Newton iteration ---
         dsig_approx  = 0.0_dp
         sig_saved    = state%sig
         statev_saved = state%statev

         ! ---- do_kiter Newton loop ----------------------------------------
         do_kiter: do kiter = 1, maxiter
            deps_corr = 0.0_dp

            if (is_obey) then
               ! ObeyRestrictions: constraints are applied directly in Cartesian space
               ddsig_by_ddeps_bar = matmul(config%cMt, ddsig_by_ddeps) + config%cMe
               dsig_target = -matmul(config%cMt, dsig_approx) - matmul(config%cMe, deps) + config%mbinc
               call USOLVER(ddsig_by_ddeps_bar, deps_corr, dsig_target, config%ifstress, ntens)
               deps      = deps + deps_corr
               deps_Cart = deps
            else
               ! Standard mixed stress/strain: work in Roscoe coords for USOLVER
               dsig_target = 0.0_dp
               where (config%ifstress == 1) dsig_target = dsig - dsig_approx
               ddsig_by_ddeps_bar = matmul(matmul(M, ddsig_by_ddeps), transpose(M))
               call USOLVER(ddsig_by_ddeps_bar, deps_corr, dsig_target, config%ifstress, ntens)
               where (config%ifstress == 1) deps = deps + deps_corr
               deps_Cart = matmul(transpose(M), deps)   ! M^{-T} * deps_Rosc → deps_Cart
            end if

            call runner%run(state, deps_Cart, ddsig_by_ddeps, ndi, nshr, ntens)

            if (kiter < maxiter) then
               ! Undo model updates — accumulate stress approximation for next iteration
               state%statev = statev_saved
               if (is_obey) then
                  dsig_approx = state%sig - sig_saved
               else
                  sig_Rosc       = matmul(M, state%sig)
                  sig_Rosc_saved = matmul(M, sig_saved)
                  where (config%ifstress == 1) dsig_approx = sig_Rosc - sig_Rosc_saved
               end if
               state%sig = sig_saved
            else
               ! Accept: advance total strain
               state%eps = state%eps + deps_Cart
            end if
         end do do_kiter

         call check_stress_inc_size(dsig_approx, dsig_target)

         ! --- *DeformationGradient: rigid rotation of stress and strain ---
         if (trim(config%load_type) == '*DeformationGradient') then
            T33         = map2T(state%sig, ntens)
            T33         = matmul(matmul(R_polar, T33), transpose(R_polar))
            state%sig   = map2stress(T33, ntens)
            eps33       = map2D(state%eps, ntens)
            eps33       = matmul(matmul(R_polar, eps33), transpose(R_polar))
            state%eps   = map2stran(eps33, ntens)
         end if

         ! --- Clamp near-zero values to avoid underflow in output ---
         state%sig    = set_zero_with_tol(state%sig)
         state%eps    = set_zero_with_tol(state%eps)
         state%statev = set_zero_with_tol(state%statev)

         ! --- Capture result snapshot (end-of-increment state) ---
         n_results = n_results + 1
         results_buf(n_results) = state
         results_buf(n_results)%time = state%time + [dt, dt]
         results_buf(n_results)%dt   = dt

         ! --- *Perturbations: undo increment so perturbations don't accumulate ---
         if (trim(config%load_type) == '*PerturbationsS' .or. &
             trim(config%load_type) == '*PerturbationsE') then
            state%eps    = state%eps - deps_Cart
            state%statev = statev_saved
            state%sig    = sig_saved
         end if

         ! --- Advance time and temperature ---
         state%time(1) = state%time(1) + dt
         state%time(2) = state%time(2) + dt
         state%temp    = state%temp + dtemp

         ! --- Exit condition ---
         if (config%has_exit_cond) then
            if (EXITNOW(config%exit_cond, state%sig, state%eps, state%statev, size(state%statev))) &
               exit do_kinc
         end if

         ! --- *ImportFile: stress alignment ---
         if (trim(config%load_type) == '*ImportFile' .and. present(align)) then
            call tryAlignStress(align, kinc, newState, config%n_import, state%sig, ntens)
         end if

      end do do_kinc
      ! ======================================================================

      ! Close ImportFile if it was opened (guard against early-exit path
      ! where close already happened on EOF)
      if (trim(config%load_type) == '*ImportFile' .and. n_results == config%n_inc) then
         close(import_file_id)
      end if

      ! Trim results to actual count (early exit reduces n_results)
      results = results_buf(1:n_results)

   end subroutine integrate_step

end module indr_step_runner
