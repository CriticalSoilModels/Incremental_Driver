!! Incremental Driver — library stub.
!!
!! This repo is a library. Include it as an fpm dependency and provide your
!! own UMAT:
!!
!!   [dependencies]
!!   Incremental_Driver = { git = "https://github.com/..." }
!!
!! Then in your main program:
!!
!!   use indr_run_model,   only: run_model        ! file-based API
!!   use indr_step_runner, only: integrate_step   ! file-free API
!!   use indr_umat_runner, only: umat_runner_t
!!
!!   call run_model(your_umat)
!!
!! See example/elastic_umat_demo.f90 for a working file-free example.
program incremental_driver
   implicit none
   print *, 'Incremental Driver is a library — see the source for usage.'
end program incremental_driver
