!! Mixed stress/strain Newton solver for the incremental driver.
module indr_solver
   use stdlib_kinds, only: dp
   implicit none
   private
   public :: USOLVER

contains

   ! Solve KK*u = rhs for mixed stress/strain boundary conditions.
   !   is(i) = 1  ->  rhs(i) is prescribed (stress-controlled)
   !   is(i) = 0  ->  u(i)   is prescribed (strain-controlled)
   ! KK is not modified. Both u and rhs are intent(inout):
   ! on exit, u holds all displacements and rhs all forces.
   subroutine USOLVER(KK, u, rhs, is, ntens)
      implicit none
      integer,  intent(in)    :: ntens
      integer,  intent(in)    :: is(ntens)
      real(dp), intent(in)    :: KK(ntens,ntens)
      real(dp), intent(inout) :: u(ntens), rhs(ntens)

      real(dp) :: rhs1(ntens)
      real(dp), allocatable :: rhsPrim(:), KKprim(:,:), uprim(:)
      integer,  allocatable :: is1(:)
      integer :: i, j, ii, nis

      nis = sum(is)

      if (all(is == 0)) then
         rhs = matmul(KK, u)
         return
      end if

      if (all(is == 1)) then
         u = xLittleUnsymmetricSolver(KK, rhs)
         return
      end if

      rhs1 = rhs
      do i = 1, ntens
         if (is(i) == 0) rhs1 = rhs1 - u(i)*KK(:,i)
      end do

      allocate(KKprim(nis,nis), rhsPrim(nis), uprim(nis), is1(nis))

      ii = 0
      do i = 1, ntens
         if (is(i) == 1) then
            ii = ii + 1
            is1(ii) = i
         end if
      end do

      do i = 1, nis
         rhsPrim(i) = rhs1(is1(i))
         do j = 1, nis
            KKprim(i,j) = KK(is1(i), is1(j))
         end do
      end do

      if (nis == 1) uprim = rhsPrim / KKprim(1,1)
      if (nis  > 1) uprim = xLittleUnsymmetricSolver(KKprim, rhsPrim)

      do i = 1, nis
         u(is1(i)) = uprim(i)
      end do
      do i = 1, ntens
         if (is(i) == 0) rhs(i) = dot_product(KK(i,:), u)
      end do
      deallocate(KKprim, rhsPrim, uprim, is1)

   contains

      ! LU decomposition with partial pivoting (Numerical Recipes).
      subroutine ludcmp(a, indx, d)
         implicit none
         real(dp), intent(inout) :: a(:,:)
         integer,  intent(out)   :: indx(:)
         real(dp), intent(out)   :: d
         real(dp) :: vv(size(a,1)), aux(size(a,1))
         integer, dimension(1) :: imaxlocs
         real(dp), parameter :: TINY = 1.0e-20_dp
         integer :: j, n, imax

         n = size(a,1)
         d = 1.0_dp
         vv = maxval(abs(a), dim=2)
         if (any(vv == 0.0_dp)) stop 'singular matrix in ludcmp'
         vv = 1.0_dp / vv
         do j = 1, n
            imaxlocs = maxloc(vv(j:n) * abs(a(j:n,j)))
            imax = (j-1) + imaxlocs(1)
            if (j /= imax) then
               aux = a(j,:);  a(j,:) = a(imax,:);  a(imax,:) = aux
               d = -d
               vv(imax) = vv(j)
            end if
            indx(j) = imax
            if (a(j,j) == 0.0_dp) a(j,j) = TINY
            a(j+1:n,j) = a(j+1:n,j) / a(j,j)
            a(j+1:n,j+1:n) = a(j+1:n,j+1:n) &
               - spread(a(j+1:n,j),2,n-j) * spread(a(j,j+1:n),1,n-j)
         end do
      end subroutine ludcmp

      ! LU back-substitution (Numerical Recipes).
      subroutine lubksb(a, indx, b)
         implicit none
         real(dp), intent(in)    :: a(:,:)
         integer,  intent(in)    :: indx(:)
         real(dp), intent(inout) :: b(:)
         integer  :: i, n, ii, ll
         real(dp) :: summ

         n = size(a,1);  ii = 0
         do i = 1, n
            ll = indx(i);  summ = b(ll);  b(ll) = b(i)
            if (ii /= 0) then
               summ = summ - dot_product(a(i,ii:i-1), b(ii:i-1))
            else if (summ /= 0.0_dp) then
               ii = i
            end if
            b(i) = summ
         end do
         do i = n, 1, -1
            b(i) = (b(i) - dot_product(a(i,i+1:n), b(i+1:n))) / a(i,i)
         end do
      end subroutine lubksb

      ! One step of iterative improvement for LU solution.
      subroutine mprove(a, alud, indx, b, x)
         implicit none
         real(dp), intent(in)    :: a(:,:), alud(:,:)
         integer,  intent(in)    :: indx(:)
         real(dp), intent(in)    :: b(:)
         real(dp), intent(inout) :: x(:)
         real(dp) :: r(size(a,1))
         r = matmul(a, x) - b
         call lubksb(alud, indx, r)
         x = x - r
      end subroutine mprove

      ! Solve a*x = b without modifying a or b (LU + iterative improvement).
      function xLittleUnsymmetricSolver(a, b) result(x)
         implicit none
         real(dp), intent(in)    :: a(:,:)
         real(dp), intent(inout) :: b(:)
         real(dp) :: x(size(b))
         real(dp) :: aa(size(b),size(b))
         integer  :: indx(size(b))
         real(dp) :: d

         x  = b;  aa = a
         call ludcmp(aa, indx, d)
         call lubksb(aa, indx, x)
         call mprove(a, aa, indx, b, x)
      end function xLittleUnsymmetricSolver

   end subroutine USOLVER

end module indr_solver
