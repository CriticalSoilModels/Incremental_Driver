module incremental_driver
   !! Public API for the Incremental Driver library.
   !!
   !! This is the only module library users need to import. It re-exports all
   !! types, procedures, and interfaces needed to drive a constitutive model:
   !!
   !!   use incremental_driver
   !!
   !! or selectively:
   !!
   !!   use incremental_driver, only: integrate_step, material_state_t, umat_runner_t

   use indr_types,        only: step_config_t, material_state_t, &
                                StressAlignment, umat_interface, &
                                STRAIN_CTRL, STRESS_CTRL
   use indr_model_runner, only: model_runner_t
   use indr_umat_runner,  only: umat_runner_t
   use indr_step_runner,  only: integrate_step
   use indr_run_model,    only: run_model
   use indr_step_parser,  only: step_record_t, parse_test_file
   use indr_file_io,      only: load_import_data, write_step_output,           &
                                read_parameter_file, read_init_conditions_file, &
                                write_output_file_header, write_line_output_data

   implicit none

   private
   public :: step_config_t, material_state_t, StressAlignment, umat_interface, &
             STRAIN_CTRL, STRESS_CTRL, &
             model_runner_t, umat_runner_t,                                     &
             integrate_step, run_model,                                         &
             step_record_t, parse_test_file,                                    &
             load_import_data, write_step_output,                               &
             read_parameter_file, read_init_conditions_file,                    &
             write_output_file_header, write_line_output_data

end module incremental_driver
