!! Concrete model runner that wraps an Abaqus UMAT procedure.
module indr_umat_runner
   use stdlib_kinds, only: dp
   use indr_step_params, only: material_state_t, umat_interface
   use indr_model_runner, only: model_runner_t
   implicit none(type, external)
   private
   public :: umat_runner_t

   type, extends(model_runner_t) :: umat_runner_t
      !! Wraps a standard Abaqus UMAT procedure pointer.
      procedure(umat_interface), pointer, nopass :: proc => null()
      real(dp), allocatable :: props(:)    !! material properties
      character(80)         :: cmname = '' !! material name (80-char Abaqus convention)
      integer               :: nstatv = 0  !! number of state variables
   contains
      procedure :: run => run_umat
   end type umat_runner_t

contains

   subroutine run_umat(this, state, deps, ddsig_by_ddeps, ndi, nshr, ntens)
      !! Unpack state into the UMAT argument list, call this%proc, and pack
      !! the outputs (sig, statev, ddsig_by_ddeps) back into state.
      !!
      !! Fixed UMAT bookkeeping arguments (noel, npt, layer, kspt, coords,
      !! celent, predef, dpred) are set to standard dummy values since the
      !! incremental driver runs outside of any FEA mesh.
      class(umat_runner_t),  intent(inout) :: this
      type(material_state_t), intent(inout) :: state
      real(dp), intent(in)  :: deps(ntens)
      real(dp), intent(out) :: ddsig_by_ddeps(ntens, ntens)
      integer,  intent(in)  :: ndi, nshr, ntens

      ! UMAT scratch outputs — not needed by the driver
      real(dp) :: sse, spd, scd, rpl, drpldt, pnewdt
      real(dp) :: ddsddt(ntens), drplde(ntens)

      ! Fixed dummy values for FEA bookkeeping arguments
      real(dp), parameter :: coords(3)     = 0.0_dp
      real(dp), parameter :: predef(1)     = 0.0_dp
      real(dp), parameter :: dpred(1)      = 0.0_dp
      real(dp), parameter :: celent        = 1.0_dp
      ! drot is the incremental rotation from polar decomposition of F.
      ! integrate_step handles *DeformationGradient rotation separately;
      ! here we pass identity (no rigid rotation this call).
      real(dp), parameter :: drot(3,3) = reshape( &
         [1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      integer,  parameter :: noel  = 1, npt   = 1
      integer,  parameter :: layer = 1, kspt  = 1
      integer,  parameter :: kstep = 1

      integer :: nprops, kinc

      nprops = size(this%props)
      kinc   = int(state%time(1) / max(state%dt, tiny(state%dt)))

      call this%proc( &
         state%sig, state%statev, ddsig_by_ddeps, sse, spd, scd,           &
         rpl, ddsddt, drplde, drpldt,                                       &
         state%eps, deps, state%time, state%dt, state%temp, 0.0_dp,        &
         predef, dpred, this%cmname,                                        &
         ndi, nshr, ntens, this%nstatv, this%props, nprops,                 &
         coords, drot, pnewdt,                                              &
         celent, state%F_start, state%F_end, noel, npt, layer, kspt, kstep, kinc)

   end subroutine run_umat

end module indr_umat_runner
