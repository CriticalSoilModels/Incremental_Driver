module mod_maps
    use stdlib_kinds, only: dp

   implicit none
   
contains
    !   contained in program\_that\_calls\_umat writes a 6x6 matrix for debugging with Mma
   subroutine write66(a)
      implicit none
      real(dp),dimension(6,6) :: a,aT
      aT = Transpose(a)
      open(12,file='nic.m',access='append')
      write(12,'(6ha66={ ,( 2h{  ,5(f15.4,2h,  ),f15.4, 3h}, ))' ) aT
      close(12)
   end subroutine write66

   !   contained in program\_that\_calls\_umat writes a 6x1 matrix  for debugging with Mma
   subroutine write6(a)
      implicit none
      real(dp), dimension(6) :: a
      open(12,file='nic.m',access='append')
      write(12,'( 5hx6={ , 5(f15.4,2h,  ),f15.4, 3h}  )' ) a
      close(12)
   end subroutine write6


   !   contained in  program\_that\_calls\_umat converts D(3,3)  to stran(6)
   function map2stran(a,ntens)
      implicit none             !===converts D(3,3)  to stran(6) with $\gamma_{12} = 2 \epsilon_{12}$ etc.
      real(dp), intent(in), dimension(1:3,1:3) :: a
      integer, intent(in) :: ntens
      real(dp),  dimension(1:ntens) :: map2stran
      real(dp), dimension(1:6) :: b
      b =[a(1,1),a(2,2),a(3,3),2*a(1,2),2*a(1,3),2*a(2,3)]
      map2stran(1:ntens)=b(1:ntens)
   end function map2stran

   !   contained in  program\_that\_calls\_umat converts strain rate from vector dstran(1:ntens) to  D(3,3)
   function map2D(a,ntens)
      implicit none
      real(dp),  dimension(1:3,1:3) :: map2D
      integer, intent(in) :: ntens
      real(dp), intent(in), dimension(:) :: a
      real(dp),dimension(1:6) :: b = 0
      b(1:ntens) = a(1:ntens)
      map2D = reshape( [b(1), b(4)/2, b(5)/2, b(4)/2,b(2),b(6)/2, b(5)/2,b(6)/2, b(3)],[3,3] )
   end function map2D

   !   contained in  program\_that\_calls\_umat converts tensor T(3,3)  to matrix stress(ntens)
   function map2stress(a,ntens)
      implicit none
      real(dp), intent(in), dimension(1:3,1:3) :: a
      integer, intent(in) :: ntens
      real(dp),  dimension(1:ntens) :: map2stress
      real(dp), dimension(1:6) :: b
      b = [a(1,1),a(2,2),a(3,3),a(1,2),a(1,3),a(2,3)]
      map2stress = b(1:ntens)
   end function map2stress

   !   contained in  program\_that\_calls\_umat converts matrix stress(1:ntens)  to tensor T(3,3)
   function map2T(a,ntens)
      implicit none
      real(dp),  dimension(1:3,1:3) :: map2T
      integer, intent(in) :: ntens
      real(dp), intent(in), dimension(:) :: a
      real(dp), dimension(1:6) :: b= 0
      b(1:ntens) = a(1:ntens)
      map2T = reshape( [b(1),b(4),b(5), b(4),b(2),b(6),  b(5),b(6),b(3) ],[3,3] )
   end function map2T

end module mod_maps