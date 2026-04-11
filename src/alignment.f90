module indr_alignment
    use stdlib_kinds, only: dp
    use indr_types, only: StressAlignment
    use indr_parser, only: splitaLine
   implicit none
   private
   public :: readAlignment, tryAlignStress

contains
   !   contained in  program\_that\_calls\_umat  reads a file with instructions for stress alignment
   subroutine  readAlignment(align, ImportFileName )
      implicit none
      character(len=40) ImportFileName, trunc, extension
      character(len=80) ReversalFileName
      logical ::  okSplit
      type(StressAlignment) :: align
      call splitaLine( ImportFileName ,'.',trunc, extension, okSplit )
      if(.not. okSplit) stop 'error  readAlignment FileName without . '
      reversalFileName = Trim(trunc) // 'rev'
      open(22, file=reversalFileName,status ='old', err=555 )
      align%active=.True.
      align%reversal(:) = 0
      read(22,*,err=556)  align%kblank,align%nrec,align%kReversal, align%ncol
      read(22,*,err=557)  align%reversal(1:align%kReversal)
      read(22,*,err=558)  align%isig(1:6)
      read(22,*,err=559)  align%sigFac(1:6)
      return
555   align%active=.False.
      return
556   stop 'error   readAlignment  cannot read kblank... '
557   stop 'error   readAlignment  cannot read reversal() '
558   stop 'error   readAlignment  cannot read sigCol() '
559   stop 'error   readAlignment  cannot read factor() '
   end subroutine  readAlignment

!   contained in  program\_that\_calls\_umat tries to align stress to values from aState(1:mImport)
   subroutine  tryAlignStress(align, kinc, aState, mImport,stress,ntens)
      implicit none
      integer:: mImport,kinc,ntens,ie
      real(dp) :: aState(mImport)
      real(dp) :: stress(ntens)
      type(StressAlignment) :: align

      if(.not. align%active) return
      if(.not. any(align%Reversal == kinc)) return

      ! only  stress components for which isig(ie) /= 0 will be aligned
      forall(ie=1:ntens, align%isig(ie) /= 0) stress(ie)= aState( align%isig(ie))*align%sigFac(ie)
      return
   end subroutine  tryAlignStress
end module indr_alignment
