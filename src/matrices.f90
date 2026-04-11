module indr_matrices
   !! Module contains some predifined matrices
   use stdlib_kinds, only: dp

   implicit none

   private
   public :: MRoscI, MRoscImt, MRendul, MRendulmT, MRosc, MRoscmT, MCart, MCartmT

   real(dp), parameter :: sq3 = 1.7320508075688772935_dp, &
                          sq6 = 2.4494897427831780982_dp, &
                          sq2 = 1.4142135623730950488_dp

   real(dp), parameter :: i3   = 0.3333333333333333333_dp, &
                          i2   = 0.5_dp,                   &
                          isq2 = 1/sq2,                    &
                          isq3 = 1.0_dp/sq3,               &
                          isq6 = 1.0_dp/sq6

   !  M for isomorphic Roscoe variables P,Q,Z,....
   real(dp), parameter, dimension(1:6, 1:6) :: MRoscI = reshape( &
      (/ -isq3,   -2.0_dp*isq6,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         -isq3,    isq6,         -isq2,   0.0_dp,  0.0_dp,  0.0_dp, &
         -isq3,    isq6,          isq2,   0.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,        0.0_dp,  1.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,        0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,        0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp /), &
      (/6, 6/))

   real(dp), parameter, dimension(1:6, 1:6) :: MRoscImT = MRoscI     !  latest $\cM^{-T}$ (is orthogonal)

   !  M for isomorphic Rendulic $ sigma_{11}= -T_{11}$,  $sigma_{22}  = -(T_{22} + T_{33}) / \sqrt(2) $,  $ Z= \dots$
   real(dp), parameter, dimension(1:6, 1:6) :: MRendul = reshape( &
      (/ -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp, -isq2,   -isq2,    0.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp, -isq2,    isq2,    0.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp /), &
      (/6, 6/))

   real(dp), parameter, dimension(1:6, 1:6) :: MRendulmT = MRendul   !  latest  $\cM^{-T}$   (is orthogonal)

   !  M for Roscoe variables $p,q,z,....$
   real(dp), parameter, dimension(1:6, 1:6) :: MRosc = reshape( &
      (/ -i3,     -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         -i3,      i2,     -1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         -i3,      i2,      1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp /), &
      (/6, 6/))

   !  latest  $\cM^{-T}$   (is not orthogonal)
   real(dp), parameter, dimension(1:6, 1:6) :: MRoscmT = reshape( &
      (/ -1.0_dp, -2.0_dp*i3,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         -1.0_dp,  i3,         -i2,     0.0_dp,  0.0_dp,  0.0_dp, &
         -1.0_dp,  i3,          i2,     0.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,      0.0_dp,  1.0_dp,  0.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,      0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp, &
          0.0_dp,  0.0_dp,      0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp /), &
      (/6, 6/))

   !  M for Cartesian coords $T_{11}, T_{22}, T_{33}, T_{12},.....$
   real(dp), parameter, dimension(1:6, 1:6) :: MCart = reshape( &
      (/ 1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp,  0.0_dp, &
         0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp,  1.0_dp /), &
      (/6, 6/))

   real(dp), parameter, dimension(1:6, 1:6) :: MCartmT = MCart       !  latest  $\cM^{-T}$  (is orthogonal)

end module indr_matrices
