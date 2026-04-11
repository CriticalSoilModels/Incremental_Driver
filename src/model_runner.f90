!! Abstract model runner interface.
!! Concrete implementations (umat_runner_t, mcss_runner_t) extend this type.
module indr_model_runner
   use stdlib_kinds, only: dp
   use indr_types, only: material_state_t
   implicit none(type, external)
   private
   public :: model_runner_t

   type, abstract :: model_runner_t
      !! Abstract base type for constitutive model integrators.
      !! Implementations wrap either an Abaqus UMAT procedure or an
      !! object-oriented model (e.g. euler_substep from critical-soil-models).
   contains
      procedure(model_run_i), deferred :: run
   end type model_runner_t

   abstract interface
      subroutine model_run_i(this, state, deps, ddsig_by_ddeps, ndi, nshr, ntens, dtemp)
         !! Integrate the model one increment.
         !!
         !! On entry:  state%sig, state%eps, state%statev, state%time, state%dt,
         !!            state%temp, state%F_start, state%F_end are the current values.
         !! On exit:   state%sig and state%statev are updated; other fields unchanged.
         !!
         !! deps(ntens)              -- strain increment to apply [-]
         !! ddsig_by_ddeps(ntens,ntens) -- ∂σ/∂ε returned by the model.
         !!   In practice this is the elastic stiffness (possibly stress-state
         !!   dependent), not the full consistent elastoplastic tangent.
         !!   USOLVER compensates via iteration.
         !! dtemp  -- temperature increment for this call [°C]
         import model_runner_t, material_state_t, dp
         class(model_runner_t),  intent(inout) :: this
         type(material_state_t), intent(inout) :: state
         real(dp), intent(in)  :: deps(ntens)
         real(dp), intent(out) :: ddsig_by_ddeps(ntens, ntens)
         integer,  intent(in)  :: ndi, nshr, ntens
         real(dp), intent(in)  :: dtemp
      end subroutine model_run_i
   end interface

end module indr_model_runner
