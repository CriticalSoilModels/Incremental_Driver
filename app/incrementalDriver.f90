! Copyright (C)  2007 ... 2023  Andrzej Niemunis
!
! incrementalDriver is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! incrementalDriver is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.


! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,USA.

! August 2007: utility routines added
! Februar 2008 :   *Repetition repaired again
! March 2008 :    *ObeyRestrictions added
! April 2008 : command line options: test= , param= , ini= , out= , verbose=  added
! November 2008: warn if  divergent equil. iter.  $\| u\_dstress \|$ too large
! February 2009: *ObeyRestrictions works with *Repetitions
! August 2010 increases output precision   +  keyword(3) error interception
! December 2015: fragments of niemunis\_tools\_lt  unsymmetric\_module incorporated into a single file incrementalDriver.f
! January 2016 if nstatev = 0 no statev( ) values will be read in but we set nstatv= 1 and statev(1)= 0.0d0
! October 2016 exit from step on inequality condition, import step loading data from a file, write every n-th state only
! November 2016 alignment of stress
! Jan 2017 exit on inequality  corrected  twice
! June 2017 undo changes in stress and state from ZERO call of  umat  (just for jacobian, with zero dstran and zero dtime )
! Dec  2017 parser disregards comments beyond \#
! Sept 2019   c\_dstran(:) = 0   in line  590  otherwise c\_dstran may be used  before being initialized.
! Sept 2019   random walk
! 2020    exit on inequality condition extended to mixed components
! Jan 2023  the initial temperature and the increase of temperature per step can be prescribed
!           (tested with thermal expansion or thermal stress in elasticity  under *ObeyRestrictions only )

!  Main program that\_calls\_umat ( performs calculation writing  to output.txt).
PROGRAM that_calls_umat   ! written by  A.Niemunis  2007 - 2023
   use stdlib_kinds, only: dp
   use stdlib_io, only: open
   use csv_module, only: csv_file
   use critical_soil_models, only: UMAT => UMAT_MCSS
   use mod_run_model, only: run_model
   ! use MOD_MCSS_ESM, only: UMAT_MohrCoulombStrainSoftening

   implicit none
   
   real(dp) :: temps(3) = [20.0_dp, 25.0_dp, 22.0_dp]
   type(csv_file) :: f
   character(len=30), dimension(:), allocatable :: header
   integer, dimension(:), allocatable :: itypes
   logical :: status_ok

   call run_model(umat)

   ! call f%read("output.txt", header_row = 1, status_ok = status_ok)
   ! call f%get_header(header, status_ok)
   ! call f%variable_types(itypes, status_ok)

   ! get some data

end program that_calls_umat

