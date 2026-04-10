module indr_step_params
   use stdlib_kinds, only: dp
   use indr_constants, only: max_fname_len, voigt_len

   implicit none

   type step_config_t
      integer  :: n_inc, max_iter, ifstress(voigt_len), columns_in_file(7), n_import
      real(dp) :: delta_load_circ(voigt_len), phase0(voigt_len), delta_load(9), &
                  dfgrd0(3,3), dfgrd1(3,3), delta_time, import_factor(7),       &
                  delta_temp
      character(40) :: load_type, coord_sys, exit_cond, import_file
      real(dp), dimension(1:6,1:6) :: cMt, cMe
      real(dp), dimension(1:6)     :: mbinc
      logical :: has_exit_cond
   end type step_config_t

contains

   pure function set_repetition_params(ninc, maxiter, ifstress, deltaLoadCirc, phase0, &
      deltaLoad, dfgrd0, dfgrd1, deltaTime, keywords, cMe, cMt, mbinc, deltaTemp,      &
      exitCond, existCond, ImportFileName, mImport, columnsInFile, importFactor) result(ofStep)
      !! Pack individual step variables into a step_config_t for repetition storage.
      integer,                   intent(in) :: ninc
      integer,                   intent(in) :: maxiter
      integer,                   intent(in) :: ifstress(voigt_len)
      real(dp),                  intent(in) :: deltaLoadCirc(6)
      real(dp),                  intent(in) :: phase0(6)
      real(dp),                  intent(in) :: deltaLoad(9)
      real(dp),                  intent(in) :: dfgrd0(3,3)
      real(dp),                  intent(in) :: dfgrd1(3,3)
      real(dp),                  intent(in) :: deltaTime
      character(len=10),         intent(in) :: keywords(3)
      real(dp),                  intent(in) :: cMe(6,6)
      real(dp),                  intent(in) :: cMt(6,6)
      real(dp),                  intent(in) :: mbinc(6)
      real(dp),                  intent(in) :: deltaTemp
      character(len=40),         intent(in) :: exitCond
      logical,                   intent(in) :: existCond
      character(max_fname_len),  intent(in) :: ImportFileName
      integer,                   intent(in) :: mImport
      integer,                   intent(in) :: columnsInFile(7)
      real(dp),                  intent(in) :: importFactor(7)

      type(step_config_t) :: ofStep

      ofStep%n_inc           = ninc
      ofStep%max_iter        = maxiter
      ofStep%ifstress        = ifstress
      ofStep%delta_load_circ = deltaLoadCirc
      ofStep%phase0          = phase0
      ofStep%delta_load      = deltaLoad
      ofStep%dfgrd0          = dfgrd0
      ofStep%dfgrd1          = dfgrd1
      ofStep%delta_time      = deltaTime
      ofStep%load_type       = keywords(2)
      ofStep%coord_sys       = keywords(3)
      ofStep%cMe             = cMe
      ofStep%cMt             = cMt
      ofStep%mbinc           = mbinc
      ofStep%delta_temp      = deltaTemp
      ofStep%exit_cond       = exitCond
      ofStep%has_exit_cond   = existCond
      ofStep%import_file     = ImportFileName
      ofStep%n_import        = mImport
      ofStep%columns_in_file = columnsInFile
      ofStep%import_factor   = importFactor
   end function set_repetition_params

   subroutine get_repetition_params(ofstep, ninc, maxiter, ifstress, deltaLoadCirc, phase0, &
      deltaLoad, dfgrd0, dfgrd1, deltaTime, keywords, cMe, cMt, mbinc, deltaTemp,           &
      exitCond, existCond, ImportFileName, mImport, columnsInFile, importFactor)
      !! Unpack a step_config_t back into individual variables for re-execution.
      type(step_config_t),       intent(in)  :: ofStep
      integer,                   intent(out) :: ninc
      integer,                   intent(out) :: maxiter
      integer,                   intent(out) :: ifstress(voigt_len)
      real(dp),                  intent(out) :: deltaLoadCirc(6)
      real(dp),                  intent(out) :: phase0(6)
      real(dp),                  intent(out) :: deltaLoad(9)
      real(dp),                  intent(out) :: dfgrd0(3,3)
      real(dp),                  intent(out) :: dfgrd1(3,3)
      real(dp),                  intent(out) :: deltaTime
      character(len=10),         intent(out) :: keywords(3)
      real(dp),                  intent(out) :: cMe(6,6)
      real(dp),                  intent(out) :: cMt(6,6)
      real(dp),                  intent(out) :: mbinc(6)
      real(dp),                  intent(out) :: deltaTemp
      character(len=40),         intent(out) :: exitCond
      logical,                   intent(out) :: existCond
      character(max_fname_len),  intent(out) :: ImportFileName
      integer,                   intent(out) :: mImport
      integer,                   intent(out) :: columnsInFile(7)
      real(dp),                  intent(out) :: importFactor(7)

      ninc           = ofStep%n_inc
      maxiter        = ofStep%max_iter
      ifstress       = ofStep%ifstress
      deltaLoadCirc  = ofStep%delta_load_circ
      phase0         = ofStep%phase0
      deltaLoad      = ofStep%delta_load
      dfgrd0         = ofStep%dfgrd0
      dfgrd1         = ofStep%dfgrd1
      deltaTime      = ofStep%delta_time
      keywords(2)    = ofStep%load_type
      keywords(3)    = ofStep%coord_sys
      cMe            = ofStep%cMe
      cMt            = ofStep%cMt
      mbinc          = ofStep%mbinc
      deltaTemp      = ofStep%delta_temp
      exitCond       = ofStep%exit_cond
      existCond      = ofStep%has_exit_cond
      ImportFileName = ofStep%import_file
      mImport        = ofStep%n_import
      columnsInFile  = ofStep%columns_in_file
      importFactor   = ofStep%import_factor
   end subroutine get_repetition_params

end module indr_step_params
