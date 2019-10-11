module mod_math

contains
!********************************************************************
!    spline (numerical recipes)
!********************************************************************
subroutine spline(x, y, n, y2)
  implicit none
  integer   i, k, n, nmax
  parameter  (nmax=5000)
  real*8    yp1, ypn, x(n), y(n), y2(n), p, qn, sig, un, u(nmax)

  y2(1) = 0.
  u(1)  = 0.
  do i=2, n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.
     y2(i)=(sig-1.)/p
     u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/ &
           (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo

  qn=0.
  un=0.
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

  do k=n-1, 1, -1
     y2(k)=y2(k)*y2(k+1)+u(k)
  enddo

  return
  end subroutine spline


!********************************************************************
!     spline (numerical recipes)
!********************************************************************
subroutine splint(xa,ya,y2a,n,x,y,khi,klo)
    implicit none
    integer n,k,khi,klo
    real*8 x,y,xa(n),y2a(n),ya(n), a,b,h

     !if ((khi.eq.0) .and. (klo.eq.0)) then
    klo=1
    khi=n
 1  if (khi-klo.gt.1) then
       k=(khi+klo)/2
       if(xa(k).gt.x)then
          khi=k
       else
          klo=k
       endif
       goto 1
    endif

    h=xa(khi)-xa(klo)
    if (h.eq.0.) stop 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
    return
end subroutine splint

end module mod_math
