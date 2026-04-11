module indr_types
   !! Core data types for the incremental driver.
   use stdlib_kinds, only: dp
   use indr_constants, only: max_fname_len, voigt_len

   implicit none
   private
   public :: step_config_t, material_state_t, umat_interface, StressAlignment, &
             STRAIN_CTRL, STRESS_CTRL

   integer, parameter :: STRAIN_CTRL = 0  !! ifstress flag: direction is strain-controlled [-]
   integer, parameter :: STRESS_CTRL = 1  !! ifstress flag: direction is stress-controlled [-]

   type step_config_t
      integer  :: n_inc, max_iter, ifstress(voigt_len), columns_in_file(7), n_import
      real(dp) :: delta_load_circ(voigt_len), phase0(voigt_len), delta_load(9), &
                  dfgrd0(3,3), dfgrd1(3,3), delta_time, import_factor(7),       &
                  delta_temp
      character(40) :: load_type, coord_sys, exit_cond, import_file
      real(dp), dimension(1:6,1:6) :: cMt, cMe
      real(dp), dimension(1:6)     :: mbinc
      logical :: has_exit_cond
      real(dp), allocatable :: import_data(:,:)
         !! Pre-loaded rows from the *ImportFile data file.
         !! Shape: (n_rows, n_import) where n_rows = total rows in file
         !! and n_import = maxval(columns_in_file).
         !! Row 1 is the initial state; row kinc+1 is the new state at increment kinc.
         !! Unallocated for all non-ImportFile load types.
   end type step_config_t

   type material_state_t
      !! Persistent model state carried across increments.
      !! Also used as a per-increment result snapshot (one entry per results(:) element).
      real(dp)              :: sig(6)        !! Cauchy stress [kPa], Voigt
      real(dp)              :: eps(6)        !! total strain [-]
      real(dp), allocatable :: statev(:)     !! model state variables
      real(dp)              :: time(2)       !! [step time, total time] [s]
      real(dp)              :: dt            !! time increment for this snapshot [s]
      real(dp)              :: temp          !! temperature [°C]
      real(dp)              :: F_start(3,3)  !! deformation gradient, start of increment [-]
      real(dp)              :: F_end(3,3)    !! deformation gradient, end of increment [-]
   end type material_state_t

   type StressAlignment
      !! Controls stress re-alignment at specified increments during *ImportFile steps.
      logical               :: active
      character(len=40)     :: ImportFileName
      integer               :: kblank, nrec, kReversal, ncol
      integer, dimension(100) :: Reversal
      integer, dimension(6)   :: isig
      real(dp), dimension(6)  :: sigFac
   end type StressAlignment

   abstract interface
      subroutine umat_interface(stress, statev, ddsdde, sse, spd, scd, &
            rpl, ddsddt, drplde, drpldt, &
            stran, dstran, time, dtime, temp, dtemp, predef, dpred, cmname, &
            ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, &
            celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
         implicit none
         character*80 cmname
         integer :: ntens, nstatev, nprops, ndi, nshr, noel, &
            npt, layer, kspt, kstep, kinc
         real(8) :: sse, spd, scd, rpl, drpldt, dtime, temp, dtemp, &
            pnewdt, celent
         real(8) :: stress(ntens), statev(nstatev), &
            ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens), &
            stran(ntens), dstran(ntens), time(2), predef(1), dpred(1), &
            props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      end subroutine umat_interface
   end interface

end module indr_types
