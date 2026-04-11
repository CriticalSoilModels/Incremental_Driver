!! Core computation layer — file-free step integration.
!! integrate_step runs the full do_kinc loop for one load step, delegating
!! each model call to the provided runner. No file I/O occurs here.
module indr_step_runner
   use stdlib_kinds, only: dp
   use indr_step_params, only: step_config_t, material_state_t
   use indr_model_runner, only: model_runner_t
   use indr_types, only: StressAlignment
   implicit none(type, external)
   private
   public :: integrate_step

contains

   subroutine integrate_step(config, state, runner, results, align)
      !! Integrate one load step, returning a snapshot per increment in results(:).
      !!
      !! config  -- step parameters (load type, n_inc, ifstress, etc.)
      !! state   -- model state on entry; updated in-place to end-of-step values
      !! runner  -- concrete model runner (umat_runner_t or mcss_runner_t)
      !! results -- one material_state_t snapshot per increment (allocated here)
      !! align   -- optional stress alignment (passed through to tryAlignStress)
      type(step_config_t),                 intent(in)    :: config
      type(material_state_t),              intent(inout) :: state
      class(model_runner_t),               intent(inout) :: runner
      type(material_state_t), allocatable, intent(out)   :: results(:)
      type(StressAlignment),  optional,    intent(in)    :: align

      ! Stub: allocate empty results array so callers can test size == 0
      allocate(results(0))

   end subroutine integrate_step

end module indr_step_runner
