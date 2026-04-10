module indr_loads
   use stdlib_kinds, only: dp
   use indr_parser, only: ReadStepCommons, splitaLine, PARSER
   use indr_linalg, only: inv33
   use indr_alignment, only: readAlignment
   use indr_types, only: StressAlignment
   use indr_constants, only: voigt_len, max_fname_len, max_lname_len
   use indr_step_params, only: step_config_t

   implicit none
   ! private
   ! public ::
contains

   subroutine read_deformation_gradient_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq
      integer :: i

      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      config%coord_sys = '*Cartesian'
      do i = 1, 9
         read(file_id, *) config%delta_load(i)
      end do
   end subroutine read_deformation_gradient_load


   subroutine read_circulating_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq
      integer :: i
      character(40) :: keyword

      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      read(file_id, *) keyword
      config%coord_sys = trim(keyword)
      do i = 1, 6
         read(file_id, *) config%ifstress(i), config%delta_load_circ(i), &
                          config%phase0(i), config%delta_load(i)
      end do
   end subroutine read_circulating_load


   subroutine read_linear_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq
      integer :: i
      character(40) :: keyword

      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      read(file_id, *) keyword
      config%coord_sys = trim(keyword)
      do i = 1, 6
         read(file_id, *) config%ifstress(i), config%delta_load(i)
      end do
   end subroutine read_linear_load


   ! Note: config%load_type must be pre-populated with the full *ImportFile|filename|id string
   ! so that the filename can be parsed out of it.
   subroutine read_file_load(file_id, config, write_freq, align)
      integer,               intent(in)    :: file_id
      type(step_config_t),   intent(inout) :: config
      integer,               intent(out)   :: write_freq
      type(StressAlignment), intent(in)    :: align

      character(len=40) :: load_file_name, aShortLine, leftLine, rightLine, keyword2
      integer :: i, load_file_id
      logical :: okSplit

      keyword2 = config%load_type(12:)   ! strip '*ImportFile' prefix
      call splitaLine(keyword2, '|', config%import_file, load_file_name, okSplit)
      if (.not. okSplit) stop 'missing | in line *ImportFile'
      read(load_file_name, *) load_file_id

      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)

      read(file_id, '(a)') aShortLine
      config%coord_sys = trim(aShortLine)

      config%columns_in_file = 0
      config%import_factor   = 1.0_dp
      do i = 1, 6
         read(file_id, '(a)') aShortLine
         call splitaLine(aShortLine, '*', leftLine, rightLine, okSplit)
         read(leftLine, *) config%ifstress(i), config%columns_in_file(i)
         if (okSplit) read(rightLine, *) config%import_factor(i)
      end do
      if (config%delta_time <= 0) then
         read(file_id, '(a)') aShortLine
         call splitaLine(aShortLine, '*', leftLine, rightLine, okSplit)
         read(leftLine, *) config%columns_in_file(7)
         if (okSplit) read(rightLine, *) config%import_factor(7)
      end if

      call readAlignment(align, config%import_file)
   end subroutine read_file_load


   subroutine read_oedometric_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq

      config%load_type = '*LinearLoad'
      config%coord_sys = '*Cartesian'
      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      read(file_id, *) config%delta_load(1)
   end subroutine read_oedometric_load


   subroutine read_oedometric_S1_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq

      config%load_type  = '*LinearLoad'
      config%coord_sys  = '*Cartesian'
      config%ifstress(1) = 1
      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      read(file_id, *) config%delta_load(1)
   end subroutine read_oedometric_S1_load


   subroutine read_triaxial_e1_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq

      config%load_type     = '*LinearLoad'
      config%coord_sys     = '*Cartesian'
      config%ifstress(1:2) = 1
      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      read(file_id, *) config%delta_load(1)
   end subroutine read_triaxial_e1_load


   subroutine read_triaxial_s1_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq

      config%load_type     = '*LinearLoad'
      config%coord_sys     = '*Cartesian'
      config%ifstress(1:3) = 1
      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      read(file_id, *) config%delta_load(1)
   end subroutine read_triaxial_s1_load


   subroutine read_triaxial_ueq_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq

      config%load_type = '*LinearLoad'
      config%coord_sys = '*Roscoe'
      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      read(file_id, *) config%delta_load(2)
   end subroutine read_triaxial_ueq_load


   subroutine read_triaxial_uq_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq

      config%load_type  = '*LinearLoad'
      config%coord_sys  = '*Roscoe'
      config%ifstress(2) = 1
      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      read(file_id, *) config%delta_load(2)
   end subroutine read_triaxial_uq_load


   subroutine read_pure_relaxation_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq

      config%load_type = '*LinearLoad'
      config%coord_sys = '*Cartesian'
      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
   end subroutine read_pure_relaxation_load


   subroutine read_pure_creep_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq

      config%load_type = '*LinearLoad'
      config%coord_sys = '*Cartesian'
      config%ifstress  = 1
      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
   end subroutine read_pure_creep_load


   subroutine read_undrained_creep(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq

      config%load_type   = '*LinearLoad'
      config%coord_sys   = '*Roscoe'
      config%ifstress(1) = 0
      config%ifstress(2:6) = 1
      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
   end subroutine read_undrained_creep


   ! mb (constraint RHS) is not stored in step_config_t — returned separately.
   subroutine read_obey_restrictions_load(file_id, config, write_freq, mb)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq
      real(dp),            intent(out)   :: mb(6)
      character(260) :: inputline(6)
      integer :: i

      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      do i = 1, 6
         read(file_id, '(a)') inputline(i)
         if (index(inputline(i), '=') == 0) stop 'restr without "=" '
      end do
      call parser(inputline, config%cMt, config%cMe, mb)
      config%mbinc     = mb / config%n_inc
      config%coord_sys = '*Cartesian'
      config%ifstress  = 1
   end subroutine read_obey_restrictions_load


   subroutine read_perturbations_S_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq
      character(40) :: keyword

      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      config%delta_temp = 0.0_dp
      read(file_id, *) keyword
      config%coord_sys = adjustl(trim(keyword))
      if (config%coord_sys /= '*Rendulic' .and. config%coord_sys /= '*RoscoeIsomorph') &
         write(*, *) 'warning: non-Isomorphic perturbation'
      read(file_id, *) config%delta_load(1)
      config%ifstress = 1
   end subroutine read_perturbations_S_load


   subroutine read_perturbations_E_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq
      character(40) :: keyword

      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      config%delta_temp = 0.0_dp
      read(file_id, *) keyword
      config%coord_sys = adjustl(trim(keyword))
      if (config%coord_sys /= '*Rendulic' .and. config%coord_sys /= '*RoscoeIsomorph') &
         write(*, *) 'warning: Anisomorphic perturbation'
      read(file_id, *) config%delta_load(1)
   end subroutine read_perturbations_E_load


   subroutine read_random_walk_load(file_id, config, write_freq)
      integer,             intent(in)    :: file_id
      type(step_config_t), intent(inout) :: config
      integer,             intent(out)   :: write_freq
      character(40) :: keyword
      integer :: i

      call ReadStepCommons(file_id, config%n_inc, config%max_iter, &
                           config%delta_time, config%delta_temp, write_freq)
      config%delta_temp = 0.0_dp
      read(file_id, *) keyword
      config%coord_sys = adjustl(trim(keyword))
      if (config%coord_sys /= '*Rendulic' .and. config%coord_sys /= '*RoscoeIsomorph') &
         write(*, *) 'warning: Anisomorphic perturbation'
      do i = 1, 6
         read(file_id, *) config%ifstress(i), config%delta_load(i)
      end do
   end subroutine read_random_walk_load


   ! Converts a step description into per-increment ddstress/dstran.
   ! Called once per increment (not once per step).
   subroutine get_increment(config, time, dtime, ddstress, dstran, dTemp, Qb33, dfgrd0, dfgrd1, drot)
      implicit none
      type(step_config_t), intent(in)    :: config
      real(dp),            intent(in)    :: time(2)
      real(dp),            intent(out)   :: dtime, ddstress(6), dstran(6), Qb33(3,3), dTemp
      real(dp),            intent(inout) :: dfgrd0(3,3), dfgrd1(3,3), drot(3,3)

      real(dp), parameter :: Pi = 3.1415926535897932385_dp
      real(dp), parameter :: delta(3,3) = reshape([1,0,0,0,1,0,0,0,1], [3,3])
      real(dp) :: Fb(3,3), Fbb(3,3), dFb(3,3), aux33(3,3), dLb(3,3), depsb(3,3), dOmegab(3,3)
      real(dp) :: wd(6), w0(6), t, arandom
      integer  :: i

      associate(load_type       => config%load_type,       &
                delta_time      => config%delta_time,      &
                n_inc           => config%n_inc,           &
                ifstress        => config%ifstress,        &
                delta_load      => config%delta_load,      &
                delta_load_circ => config%delta_load_circ, &
                phase0          => config%phase0,          &
                delta_temp      => config%delta_temp)

      dtime    = delta_time / n_inc
      dTemp    = delta_temp / n_inc
      dstran   = 0.0_dp
      ddstress = 0.0_dp
      Qb33     = delta
      drot     = delta
      dfgrd0   = delta
      dfgrd1   = delta

      if (load_type == '*LinearLoad') then
         do i = 1, 6
            if (ifstress(i) == 1) ddstress(i) = delta_load(i) / n_inc
            if (ifstress(i) == 0) dstran(i)   = delta_load(i) / n_inc
         end do
      end if

      if (load_type == '*DeformationGradient') then
         Fb = reshape([delta_load(1), delta_load(5), delta_load(7), &
                       delta_load(4), delta_load(2), delta_load(9), &
                       delta_load(6), delta_load(8), delta_load(3)], [3,3])
         Fbb    = delta + (Fb - delta) * (time(1) / delta_time)
         dfgrd0 = Fbb
         dFb    = (Fb - delta) / n_inc
         aux33  = Fbb + dFb / 2.0_dp
         dfgrd1 = Fbb + dFb
         aux33  = inv33(aux33)
         dLb    = matmul(dFb, aux33)
         depsb  = 0.5_dp * (dLb + transpose(dLb))
         dstran = [depsb(1,1), depsb(2,2), depsb(3,3), &
                   2.0_dp*depsb(1,2), 2.0_dp*depsb(1,3), 2.0_dp*depsb(2,3)]
         dOmegab = 0.5_dp * (dLb - transpose(dLb))
         aux33   = inv33(delta - 0.5_dp*dOmegab)
         Qb33    = matmul(aux33, delta + 0.5_dp*dOmegab)
         drot    = Qb33
      end if

      if (load_type == '*CirculatingLoad') then
         wd = 2.0_dp * Pi / delta_time
         w0 = phase0
         t  = time(1) + dtime / 2.0_dp
         do i = 1, 6
            if (ifstress(i) == 1) ddstress(i) = dtime*delta_load_circ(i)*wd(i)*cos(wd(i)*t + w0(i)) + delta_load(i)/n_inc
            if (ifstress(i) == 0) dstran(i)   = dtime*delta_load_circ(i)*wd(i)*cos(wd(i)*t + w0(i)) + delta_load(i)/n_inc
         end do
      end if

      if (load_type == '*PerturbationsS') then
         ddstress(1) = delta_load(1) * cos(time(1) * 2.0_dp * Pi / delta_time)
         ddstress(2) = delta_load(1) * sin(time(1) * 2.0_dp * Pi / delta_time)
      end if

      if (load_type == '*PerturbationsE') then
         dstran(1) = delta_load(1) * cos(time(1) * 2.0_dp * Pi / delta_time)
         dstran(2) = delta_load(1) * sin(time(1) * 2.0_dp * Pi / delta_time)
      end if

      if (load_type == '*RandomWalk') then
         call random_seed()
         do i = 1, 6
            call random_number(arandom)
            if (ifstress(i) == 1) ddstress(i) = 2.0_dp*(arandom - 0.5_dp)*delta_load(i)
            if (ifstress(i) == 0) dstran(i)   = 2.0_dp*(arandom - 0.5_dp)*delta_load(i)
         end do
      end if

      end associate
   end subroutine get_increment

end module indr_loads
