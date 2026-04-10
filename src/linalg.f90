!! Linear algebra utilities for the incremental driver
module indr_linalg
   use stdlib_kinds, only: dp
   implicit none
   private
   public :: inv33, spectral_decomposition_of_symmetric, &
             app_jacobian_similarity, get_jacobian_rot

contains

   ! Inverts a 3x3 matrix using the analytic cofactor formula
   function inv33(a)
      implicit none
      real(dp), dimension(3,3), intent(in) :: a
      real(dp), dimension(3,3) :: inv33
      real(dp), dimension(3,3) :: b
      real(dp) :: det
      det = - a(1,3)*a(2,2)*a(3,1) + a(1,2)*a(2,3)*a(3,1) &
         + a(1,3)*a(2,1)*a(3,2) - a(1,1)*a(2,3)*a(3,2) &
         - a(1,2)*a(2,1)*a(3,3) + a(1,1)*a(2,2)*a(3,3)
      b = reshape( [-a(2,3)*a(3,2) + a(2,2)*a(3,3), a(1,3)*a(3,2) - a(1,2)*a(3,3), &
         -a(1,3)*a(2,2) + a(1,2)*a(2,3), a(2,3)*a(3,1) - a(2,1)*a(3,3), &
         -a(1,3)*a(3,1) + a(1,1)*a(3,3), a(1,3)*a(2,1) - a(1,1)*a(2,3), &
         -a(2,2)*a(3,1) + a(2,1)*a(3,2), a(1,2)*a(3,1) - a(1,1)*a(3,2), &
         -a(1,2)*a(2,1) + a(1,1)*a(2,2)], [3,3])
      inv33 = transpose(b)/det
   end function inv33

   ! Jacobi iterative eigendecomposition of a real symmetric n×n matrix.
   ! Returns eigenvalues Lam(n) and eigenvectors G(n,n) (columns).
   ! Algorithm: Kielbasinski, p. 385–386.
   subroutine spectral_decomposition_of_symmetric(A, Lam, G, n)
      implicit none
      integer, intent(in)    :: n
      real(dp), intent(in)   :: A(n,n)
      real(dp), intent(out)  :: Lam(n)
      real(dp), intent(out)  :: G(n,n)
      integer  :: iter, i, p, q
      real(dp) :: cosine, sine
      real(dp), allocatable :: pcol(:), qcol(:), x(:,:)

      allocate(pcol(n), qcol(n), x(n,n))
      x = A
      G = 0.0_dp
      do i = 1, n
         G(i,i) = 1.0_dp
      end do

      do iter = 1, 30
         call get_jacobian_rot(x, p, q, cosine, sine, n)
         call app_jacobian_similarity(x, p, q, cosine, sine, n)
         pcol = G(:,p)
         qcol = G(:,q)
         G(:,p) = pcol*cosine - qcol*sine
         G(:,q) = pcol*sine   + qcol*cosine
      end do

      do i = 1, n
         Lam(i) = x(i,i)
      end do
      deallocate(pcol, qcol, x)
   end subroutine spectral_decomposition_of_symmetric

   ! Apply one Jacobi (Givens) similarity transformation to symmetric matrix A.
   ! G_pq rotation with cosine c and sine s. Algorithm: Kielbasinski, p. 385.
   subroutine app_jacobian_similarity(A, p, q, c, s, n)
      implicit none
      integer,  intent(in)    :: p, q, n
      real(dp), intent(in)    :: c, s
      real(dp), intent(inout) :: A(n,n)
      real(dp) :: prow(n), qrow(n), App, Apq, Aqq

      if (p == q)       stop 'error: jacobian_similarity  p == q'
      if (p < 1 .or. p > n) stop 'error: jacobian_similarity p out of range'
      if (q < 1 .or. q > n) stop 'error: jacobian_similarity q out of range'

      prow = c*A(:,p) - s*A(:,q)
      qrow = s*A(:,p) + c*A(:,q)
      App  = c*c*A(p,p) - 2*c*s*A(p,q) + s*s*A(q,q)
      Aqq  = s*s*A(p,p) + 2*c*s*A(p,q) + c*c*A(q,q)
      Apq  = c*s*(A(p,p) - A(q,q)) + (c*c - s*s)*A(p,q)
      A(p,:) = prow;  A(:,p) = prow
      A(q,:) = qrow;  A(:,q) = qrow
      A(p,p) = App;   A(q,q) = Aqq
      A(p,q) = Apq;   A(q,p) = Apq
   end subroutine app_jacobian_similarity

   ! Compute optimal Jacobi rotation parameters for iterative diagonalization.
   ! Returns pivot indices p, q and rotation cosine c, sine s.
   ! Algorithm: Kielbasinski, p. 385–386.
   subroutine get_jacobian_rot(A, p, q, c, s, n)
      implicit none
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: A(n,n)
      integer,  intent(out) :: p, q
      real(dp), intent(out) :: c, s
      real(dp) :: App, Apq, Aqq, d, t, maxoff
      integer  :: i, j

      p = 0;  q = 0
      maxoff = tiny(maxoff)
      do i = 1, n-1
         do j = i+1, n
            if (abs(A(i,j)) > maxoff) then
               maxoff = abs(A(i,j))
               p = i;  q = j
            end if
         end do
      end do

      if (p > 0) then
         App = A(p,p);  Apq = A(p,q);  Aqq = A(q,q)
         d = (Aqq - App) / (2.0_dp*Apq)
         t = 1.0_dp / sign(abs(d) + sqrt(1.0_dp + d*d), d)
         c = 1.0_dp / sqrt(1.0_dp + t*t)
         s = t*c
      else
         p = 1;  q = 2;  c = 1.0_dp;  s = 0.0_dp
      end if
   end subroutine get_jacobian_rot

end module indr_linalg
