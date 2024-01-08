subroutine minmax1(a,n,amin,amax)
  implicit none
  !f2py intent(hidden) :: n
  !f2py intent(out) :: amin,amax
  !f2py intent(in) :: a
  integer n
  real a(n),amin,amax
  integer i

  amin = a(1)
  amax = a(1)
  do i=2, n
     if(a(i) > amax)then
        amax = a(i)
     elseif(a(i) < amin) then
        amin = a(i)
     endif
  enddo
end subroutine minmax1

subroutine minmax2(a,n,amin,amax)
  implicit none
  !f2py intent(hidden) :: n
  !f2py intent(out) :: amin,amax
  !f2py intent(in) :: a
  integer n
  real a(n),amin,amax
  amin = minval(a)
  amax = maxval(a)
end subroutine minmax2
