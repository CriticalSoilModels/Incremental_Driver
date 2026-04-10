program test_linalg
   ! Tests for inv33 and spectral_decomposition_of_symmetric
   use stdlib_kinds, only: dp
   use mod_inc_driver_funcs, only: inv33, spectral_decomposition_of_symmetric
   implicit none

   integer :: nfail = 0
   real(dp), parameter :: tol = 1.0e-10_dp

   call test_inv33_identity()
   call test_inv33_diagonal()
   call test_inv33_roundtrip()
   call test_spectral_diagonal()
   call test_spectral_2x2_block()

   if (nfail == 0) then
      print *, 'PASS  test_linalg'
   else
      print *, 'FAIL  test_linalg — ', nfail, ' test(s) failed'
      stop 1
   end if

contains

   subroutine test_inv33_identity()
      ! inv33(I) = I
      real(dp) :: A(3,3), B(3,3)
      integer :: i, j
      A = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      B = inv33(A)
      do i = 1, 3
         do j = 1, 3
            if (abs(B(i,j) - A(i,j)) > tol) then
               print *, 'FAIL  test_inv33_identity: B(', i, ',', j, ') =', B(i,j)
               nfail = nfail + 1
            end if
         end do
      end do
   end subroutine test_inv33_identity

   subroutine test_inv33_diagonal()
      ! inv33(diag(2,4,5)) = diag(0.5, 0.25, 0.2)
      real(dp) :: A(3,3), B(3,3)
      real(dp) :: expected(3,3)
      integer :: i, j
      A = reshape([2.0_dp,0.0_dp,0.0_dp, 0.0_dp,4.0_dp,0.0_dp, 0.0_dp,0.0_dp,5.0_dp], [3,3])
      expected = reshape([0.5_dp,0.0_dp,0.0_dp, 0.0_dp,0.25_dp,0.0_dp, 0.0_dp,0.0_dp,0.2_dp], [3,3])
      B = inv33(A)
      do i = 1, 3
         do j = 1, 3
            if (abs(B(i,j) - expected(i,j)) > tol) then
               print *, 'FAIL  test_inv33_diagonal: B(', i, ',', j, ') =', B(i,j), &
                        'expected', expected(i,j)
               nfail = nfail + 1
            end if
         end do
      end do
   end subroutine test_inv33_diagonal

   subroutine test_inv33_roundtrip()
      ! A * inv33(A) = I for a general non-singular matrix
      real(dp) :: A(3,3), Ainv(3,3), prod(3,3)
      integer :: i, j
      real(dp) :: expected
      A = reshape([1.0_dp,2.0_dp,3.0_dp, 0.0_dp,4.0_dp,5.0_dp, 1.0_dp,0.0_dp,6.0_dp], [3,3])
      Ainv = inv33(A)
      prod = matmul(A, Ainv)
      do i = 1, 3
         do j = 1, 3
            expected = merge(1.0_dp, 0.0_dp, i == j)
            if (abs(prod(i,j) - expected) > tol) then
               print *, 'FAIL  test_inv33_roundtrip: A*inv(A)(', i, ',', j, ') =', prod(i,j)
               nfail = nfail + 1
            end if
         end do
      end do
   end subroutine test_inv33_roundtrip

   subroutine test_spectral_diagonal()
      ! Diagonal matrix: eigenvalues are the diagonal entries (any order)
      real(dp) :: A(3,3), Lam(3), G(3,3)
      real(dp) :: expected(3)
      logical :: found(3)
      integer :: i, j
      A = reshape([100.0_dp,0.0_dp,0.0_dp, 0.0_dp,200.0_dp,0.0_dp, 0.0_dp,0.0_dp,300.0_dp], [3,3])
      call spectral_decomposition_of_symmetric(A, Lam, G, 3)
      expected = [100.0_dp, 200.0_dp, 300.0_dp]
      found = .false.
      do i = 1, 3
         do j = 1, 3
            if (.not. found(j) .and. abs(Lam(i) - expected(j)) < tol) then
               found(j) = .true.
               exit
            end if
         end do
      end do
      if (.not. all(found)) then
         print *, 'FAIL  test_spectral_diagonal: eigenvalues =', Lam
         nfail = nfail + 1
      end if
   end subroutine test_spectral_diagonal

   subroutine test_spectral_2x2_block()
      ! 2x2 block: [[5,2],[2,8]] embedded in 3x3 with third eigenvalue = 1
      ! Analytic eigenvalues of [[5,2],[2,8]]:
      !   lambda = (13 +/- sqrt(169 - 4*36)) / 2 = (13 +/- sqrt(25)) / 2 = 9, 4
      real(dp) :: A(3,3), Lam(3), G(3,3)
      real(dp) :: expected(3)
      logical :: found(3)
      integer :: i, j
      A = reshape([5.0_dp,2.0_dp,0.0_dp, 2.0_dp,8.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
      call spectral_decomposition_of_symmetric(A, Lam, G, 3)
      expected = [4.0_dp, 9.0_dp, 1.0_dp]
      found = .false.
      do i = 1, 3
         do j = 1, 3
            if (.not. found(j) .and. abs(Lam(i) - expected(j)) < tol) then
               found(j) = .true.
               exit
            end if
         end do
      end do
      if (.not. all(found)) then
         print *, 'FAIL  test_spectral_2x2_block: eigenvalues =', Lam, &
                  'expected [4, 9, 1]'
         nfail = nfail + 1
      end if
   end subroutine test_spectral_2x2_block

end program test_linalg
