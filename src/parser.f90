!! Parsing utilities: string splitting, step header reading,
!! *ObeyRestrictions parsing, and increment exit-condition evaluation.
module indr_parser
   use stdlib_kinds, only: dp
   implicit none
   private
   public :: splitaLine, ReadStepCommons, PARSER, EXITNOW

contains

!------------------------------------------------------------------------------------------
!  splitaLine: splits aLine at the first occurrence of sep.
!  Returns the portion left of sep in `left` and right of sep in `right`.
!  Sets ok=.true. if sep was found.
   subroutine splitaLine(aLine, sep, left, right, ok)
      implicit none
      character(len=40), intent(in)  :: aLine
      character(len=1),  intent(in)  :: sep
      character(len=40), intent(out) :: left, right
      logical,           intent(out) :: ok
      character(len=40) :: tmp
      integer :: iSep

      ok   = .false.
      iSep = index(aLine, sep)
      if (iSep == 0) then
         left  = trim(adjustl(aLine))
         right = '  '
      else
         ok    = .true.
         tmp   = aLine(:iSep-1)
         right = aLine(iSep+1:)
         left  = tmp
      end if
   end subroutine splitaLine

!------------------------------------------------------------------------------------------
!  ReadStepCommons: reads ninc, maxiter, deltaTime (and optional every, deltaTemp)
!  from the next line of an open test.inp file.
   subroutine ReadStepCommons(file_id, ninc, maxiter, deltaTime, deltaTemp, every)
      implicit none
      integer,  intent(in)  :: file_id
      real(dp), intent(out) :: deltaTime, deltaTemp
      integer,  intent(out) :: ninc, maxiter, every
      logical  :: okSplit
      character(len=40) :: aShortLine, leftLine, rightLine

      deltaTemp = 0.0_dp
      read(file_id, '(a)') aShortLine
      call splitaLine(aShortLine, ':', leftLine, rightLine, okSplit)
      read(leftLine, *) ninc, maxiter, deltaTime
      every = 1
      if (okSplit) read(rightLine, *) every
      if (every > ninc) every = ninc
      if (every < 1)    every = 1

      call splitaLine(aShortLine, 'T', leftLine, rightLine, okSplit)
      if (okSplit) read(rightLine, *) deltaTemp
   end subroutine ReadStepCommons

!------------------------------------------------------------------------------------------
!  EXITNOW: evaluates a condition string (e.g. "s1>100") against the current
!  stress, strain, and state variable vectors. Returns .true. when met.
   function EXITNOW(cond, stress, stran, statev, nstatv) result(res)
      implicit none
      integer, parameter :: ntens = 6, mSummands = 5
      integer,  intent(in) :: nstatv
      real(dp), intent(in) :: stress(ntens), stran(ntens), statev(nstatv)
      character(len=40), intent(in) :: cond
      logical :: res

      integer  :: i, igt, ilt, iis, imin, iplus, iminus, Nsummands, itimes
      character(len=40)  :: inp, rhs, summand(mSummands), aux
      real(dp) :: factor(mSummands), fac, x, y
      real(dp), parameter :: sq3  = 1.7320508075689_dp, &
                             sq23 = 0.81649658092773_dp

      res = .false.
      igt = index(cond, '>');  ilt = index(cond, '<');  iis = max(igt, ilt)
      if (iis == 0) goto 555
      inp = adjustl(cond(:iis));  rhs = trim(adjustl(cond(iis+1:)))

      factor(1) = 1
      if (inp(1:1) == '-') then
         factor(1) = -1;  inp = inp(2:)
      end if

      do i = 1, mSummands
         iplus  = index(inp, '+');  if (iplus  == 0) iplus  = 200
         iminus = index(inp, '-');  if (iminus == 0) iminus = 200
         igt    = index(inp, '>');  if (igt    == 0) igt    = 200
         ilt    = index(inp, '<');  if (ilt    == 0) ilt    = 200
         imin = min(iplus, iminus, igt, ilt)
         if (imin == 200) exit
         if (imin == iplus) then
            summand(i) = inp(:imin-1);  factor(i+1) = 1;   inp = inp(imin+1:)
         end if
         if (imin == iminus) then
            summand(i) = inp(:imin-1);  factor(i+1) = -1;  inp = inp(imin+1:)
         end if
         if (imin == ilt .or. imin == igt) then
            summand(i) = inp(:imin-1);  exit
         end if
      end do
      Nsummands = i

      x = 0.0_dp
      do i = 1, Nsummands
         aux    = adjustl(summand(i))
         itimes = index(aux, '*')
         if (itimes /= 0) then
            read(aux(:itimes-1), *) fac
            factor(i) = factor(i) * fac
            aux = trim(adjustl(aux(itimes+1:)))
         else
            aux = trim(aux)
         end if
         select case (aux)
          case ('s1');  x = x + factor(i)*stress(1)
          case ('s2');  x = x + factor(i)*stress(2)
          case ('s3');  x = x + factor(i)*stress(3)
          case ('s12'); x = x + factor(i)*stress(4)
          case ('s13'); x = x + factor(i)*stress(5)
          case ('s23'); x = x + factor(i)*stress(6)
          case ('v1');  x = x + factor(i)*statev(1)
          case ('v2');  x = x + factor(i)*statev(2)
          case ('v3');  x = x + factor(i)*statev(3)
          case ('v4');  x = x + factor(i)*statev(4)
          case ('v5');  x = x + factor(i)*statev(5)
          case ('v6');  x = x + factor(i)*statev(6)
          case ('v7');  x = x + factor(i)*statev(7)
          case ('v8');  x = x + factor(i)*statev(8)
          case ('v9');  x = x + factor(i)*statev(9)
          case ('p');   x = x - factor(i)*(stress(1)+stress(2)+stress(3))/3.0_dp
          case ('q');   x = x - factor(i)*(stress(1) - stress(3))
          case ('P');   x = x - factor(i)*(stress(1)+stress(2)+stress(3))/sq3
          case ('Q');   x = x - factor(i)*(stress(1) - stress(3))
          case ('e1');  x = x + factor(i)*stran(1)
          case ('e2');  x = x + factor(i)*stran(2)
          case ('e3');  x = x + factor(i)*stran(3)
          case ('g12'); x = x + factor(i)*stran(4)
          case ('g13'); x = x + factor(i)*stran(5)
          case ('g23'); x = x + factor(i)*stran(6)
          case ('ev');  x = x - factor(i)*(stran(1)+stran(2)+stran(3))
          case ('eq');  x = x - 2.0_dp*(stran(1) - stran(3))/3.0_dp
          case ('eP');  x = x - factor(i)*(stran(1)+stran(2)+stran(3))/sq3
          case ('eQ');  x = x - factor(i)*sq23*(stran(1) - stran(3))
          case DEFAULT; goto 555
         end select
      end do

      read(rhs, *) y
      igt = index(cond, '>');  ilt = index(cond, '<')
      if (igt /= 0) res = (x > y)
      if (ilt /= 0) res = (x < y)
      return
555   write(*,*) 'inp syntax error: ', cond, ' exit condition ignored'
      res = .false.
   end function EXITNOW

!------------------------------------------------------------------------------------------
!  PARSER: reads *ObeyRestrictions lines and fills the Mt, Me, mb matrices.
   subroutine PARSER(inputline, Mt, Me, mb)
      implicit none
      character(260),    intent(in)  :: inputline(6)
      real(dp),          intent(out) :: Mt(6,6), Me(6,6), mb(6)

      character(len=260) :: inp, aux, aux3
      character(40)      :: summand(13)
      integer  :: iis, i, iplus, iminus, iequal, imin, iex, itimes, Irestr, ihash, Nsummands
      real(dp) :: factor(13), fac

      Mt = 0.0_dp;  Me = 0.0_dp;  mb = 0.0_dp

      do Irestr = 1, 6
         inp   = trim(adjustl(inputline(Irestr)))
         ihash = index(inp, '#')
         if (ihash /= 0) inp = inp(:ihash-1)
         iis = index(inp, '=')
         if (iis == 0) stop 'parser error: no = in restriction'

         factor(1) = 1
         if (inp(1:1) == '-') then
            factor(1) = -1;  inp = inp(2:)
         end if

         do i = 1, 13
            iplus  = index(inp, '+');  if (iplus  == 0) iplus  = 200
            iminus = index(inp, '-');  if (iminus == 0) iminus = 200
            iequal = index(inp, '=');  if (iequal == 0) iequal = 200
            imin = min(iplus, iminus, iequal)
            if (imin == 200) stop 'parser err: no +,-,= in restric'
            if (imin == iplus) then
               summand(i) = inp(:imin-1);  factor(i+1) = 1;   inp = inp(imin+1:)
            end if
            if (imin == iminus) then
               summand(i) = inp(:imin-1);  factor(i+1) = -1;  inp = inp(imin+1:)
            end if
            if (imin == iequal) then
               summand(i) = inp(:imin-1)
               inp = inp(imin+1:)
               iminus = index(inp, '-');  if (iminus == 0) iminus = 200
               iex    = index(inp, '!');  if (iex    == 0) iex    = len(inp) + 1
               if (iminus == 200) then
                  factor(i+1) = 1;   summand(i+1) = inp(:iex-1)
               else
                  factor(i+1) = -1;  summand(i+1) = inp(iminus+1:iex-1)
               end if
               exit
            end if
         end do
         Nsummands = i + 1

         do i = 1, Nsummands - 1
            aux    = adjustl(summand(i))
            itimes = index(aux, '*')
            if (itimes /= 0) then
               read(aux(:itimes-1), *) fac
               factor(i) = factor(i) * fac
               aux = adjustl(aux(itimes+1:))
            end if
            aux3 = aux(1:3)
            select case (aux3)
             case ('sd1'); Mt(Irestr,1) = factor(i)
             case ('sd2'); Mt(Irestr,2) = factor(i)
             case ('sd3'); Mt(Irestr,3) = factor(i)
             case ('sd4'); Mt(Irestr,4) = factor(i)
             case ('sd5'); Mt(Irestr,5) = factor(i)
             case ('sd6'); Mt(Irestr,6) = factor(i)
             case ('ed1'); Me(Irestr,1) = factor(i)
             case ('ed2'); Me(Irestr,2) = factor(i)
             case ('ed3'); Me(Irestr,3) = factor(i)
             case ('ed4'); Me(Irestr,4) = factor(i)
             case ('ed5'); Me(Irestr,5) = factor(i)
             case ('ed6'); Me(Irestr,6) = factor(i)
            end select
         end do
         read(summand(Nsummands), *) mb(Irestr)
         mb(Irestr) = mb(Irestr) * factor(Nsummands)
      end do
   end subroutine PARSER

end module indr_parser
