module fourier
use parameters
implicit none
integer, parameter, private :: dp=kind(0.0d0)
real(kind=dp), parameter :: pi=2.d0*asin(1.d0), vk=2.0d0
!
contains
!
!========================================
subroutine basis_f(cl,ea,ew,ev,ep,eb,ich)
!========================================
real(kind=dp), intent(in) :: cl
integer, intent(in) :: ich
real(kind=dp), dimension(:), intent(out) :: ea,ew,ev,ep,eb
real(kind=dp) :: vk2,tol,z
integer :: i
vk2=vk*vk
tol=1.d-8
z=0.d0 ! a parameter used to ensure the sequence of eigenvalues is regular
do i=1,size(ev)
    z=newton_raph(i,cl,tol,z) ! z is updated
    ew(i)=z                                                 ! frequency
    if(mod(i,2)==1) then
        ea(i)=sqrt(2.d0/(1.d0+sin(ew(i)*cl)/(ew(i)*cl)))    ! amplitude
        ep(i)=0.d0                                          ! phase
    else
        ea(i)=sqrt(2.d0/(1.d0-sin(ew(i)*cl)/(ew(i)*cl)))
        ep(i)=pi/2.d0
    end if
    ev(i)=2.d0*vk/(vk2+ew(i)*ew(i))/cl                      ! eigenvalue
    eb(i)=ea(i)/ew(i)/cl                                    ! integral amplitude scaled by density
    write(ich,'(i3,5f10.4)') i,ev(i),ew(i),ep(i),ea(i),eb(i)
end do
end subroutine basis_f
!
!======================================
real(kind=dp) function efunc(n,c2,c3,x)
!======================================
integer, intent(in) :: n
real(kind=dp), intent(in) :: c2,c3,x
if(mod(n,2)==1) then ! cosine solution
    efunc=x*sin(c2*x)-c3*cos(c2*x)
else ! sine solution
    efunc=x*cos(c2*x)+c3*sin(c2*x)
end if
end function efunc
!
!======================================
real(kind=dp) function edash(n,c2,c3,x)
!======================================
integer, intent(in) :: n
real(kind=dp), intent(in) :: c2,c3,x
if(mod(n,2)==1) then 
    edash=c2*x*cos(c2*x)+(1.d0+c3*c2)*sin(c2*x)
else
    edash=(1.0d0+c3*c2)*cos(c2*x)-c2*x*sin(c2*x)
end if
end function edash
!
!=============================================
real(kind=dp) function newton_raph(n,cl,tol,x)
!=============================================
integer, intent(in) :: n
real(kind=dp), intent(in) :: cl,tol
real(kind=dp), intent(inout) :: x
real(kind=dp) :: c2,c3,dx
c2=cl/2.0d0
c3=vk
if(n==1) then
    x=pi/2.0/cl
else
    x=x+pi/cl
end if
dx=huge(x)
do while (abs(dx)>tol)
    dx=efunc(n,c2,c3,x)/edash(n,c2,c3,x)
    x=x-dx
end do
newton_raph=x
end function newton_raph
!
!=====================================
function net_value(xnet,spec,ea,ew,ep)
!=====================================
real(kind=dp), dimension(:), intent(in) :: xnet,spec,ea,ew,ep
real(kind=dp), dimension(size(xnet)) :: net_value
integer :: i
net_value(:)=0.d0
do i=1,size(spec)
    net_value(:)=net_value(:)+spec(i)*ea(i)*cos(ew(i)*xnet(:)+ep(i))
end do
end function net_value
!
end module fourier
