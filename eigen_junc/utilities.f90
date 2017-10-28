module utilities
use parameters
use sort_kit
implicit none
integer, parameter, private :: dp=kind(0.0d0)
!
contains
!    
!=================================
logical function equal_gamete(a,b)
!=================================
type(gamete) :: a,b
integer :: i
if(a%nj/=b%nj) then
    equal_gamete=.false.
else
    equal_gamete=.true.
    junctions: do i=1,a%nj
        if(equal_junction(a%junc(i),b%junc(i))) cycle
        equal_gamete=.false.
        exit junctions
    end do junctions
end if
end function equal_gamete
!
!===================================
logical function equal_junction(a,b)
!===================================
type(junction) :: a,b
if((a%anc==b%anc).and.(abs(a%pos-b%pos)<1.0d-6)) then
    equal_junction=.true.
else
    equal_junction=.false.
end if
end function equal_junction
!
!===================
subroutine matings()
!===================
integer :: i
real(kind=dp) :: u
do i=1,ni ! for each newborn individual
    select case (msys)
    case('di') ! dioecious
            call random_number(u) ! to choose a sire
            gp(i)%ped(1)=ceiling(u*real(ni/2))
            gp(i)%htype(1)=pass_gamete(gp(i)%ped(1))
            call random_number(u) ! to choose a dam
            gp(i)%ped(2)=(ni/2)+ceiling(u*real(ni/2)) ! choose a dam
            gp(i)%htype(2)=pass_gamete(gp(i)%ped(2)) ! sample a gamete from dam
        case('ms') ! monoecious selfing ... true random
            call random_number(u) ! to choose a sire
            gp(i)%ped(1)=ceiling(u*real(ni))
            gp(i)%htype(1)=pass_gamete(gp(i)%ped(1))
            call random_number(u) ! to choose a dam
            gp(i)%ped(2)=ceiling(u*real(ni)) ! choose a dam
            gp(i)%htype(2)=pass_gamete(gp(i)%ped(2)) ! sample a gamete from dam
        case('mx') ! monoecious but no selfing
            call random_number(u) ! to choose a sire
            gp(i)%ped(1)=ceiling(u*real(ni))
            gp(i)%htype(1)=pass_gamete(gp(i)%ped(1))
            call random_number(u) ! to choose a dam
            gp(i)%ped(2)=ceiling(u*real(ni-1)) ! choose a dam
            if(gp(i)%ped(2)>=gp(i)%ped(1)) gp(i)%ped(2)=gp(i)%ped(2)+1 ! shift up to avoid and account for selfing
            gp(i)%htype(2)=pass_gamete(gp(i)%ped(2)) ! sample a gamete from dam
        case default
            print *, 'error nominating mating system'; stop
    end select
end do
end subroutine matings
!
!================================
integer function crossovers(clen)
! clen is chromosome length in Morgans
!================================
real(kind=dp), intent(in) :: clen
integer :: nc
real(kind=dp) :: u,x,y
call random_number(u)
x=exp(-clen)
y=x
nc=0
do
  if(x>u) exit ! no crossovers
  nc=nc+1
  y=y*clen/real(nc)
  x=x+y
end do
crossovers=nc
end function crossovers
!
!================================
subroutine map_crossings(clen,xc)
!================================
real(kind=dp), intent(in) :: clen
real(kind=dp), dimension(:), intent(out) :: xc
real(kind=dp) :: u
integer :: i,j
call random_number(u); xc(1)=u
if(size(xc)>1) then
    next: do j=2,size(xc)
        call random_number(u)
        sort: do i=1,j-1
            if(u<xc(i)) then
                xc(i+1:j)=xc(i:j-1)
                xc(i)=u
                cycle next
            end if
            xc(j)=u
        end do sort
    end do next
end if
xc(:)=clen*xc(:) ! scale to chromosome length
end subroutine map_crossings
!
!========================================
type(gamete) function pass_gamete(parent)
!========================================
integer, intent(in) :: parent
type(gamete) :: gsel,galt,goff
integer :: m,jg,flag
real(kind=dp), dimension(:), allocatable :: xc
real(kind=dp) :: u
goff=gnull
flag=0
if(pg_test) then 
    call random_number(u) ! to determine whether this gamete is reported
    if(u<2.d-3) then
        flag=1
        write(12,'(/a,5i6)') ' ... test ... ',ir,ig,parent,xgp(parent)%htype(:)%nj
        call xgp(parent)%htype(1)%output(12)
        call xgp(parent)%htype(2)%output(12)
    end if
end if
m=crossovers(cl) ! find number of crossovers
poisson_test(min(m,10))=poisson_test(min(m,10))+1.0d0   ! accumulate for Poisson check
if(flag==1) write(12,*) ' ... crossovers ... ',m   ! for gamete check
!
call random_number(u) ! to determine gamete
jg=1; if(u<0.5d0) jg=2
gsel=xgp(parent)%htype(jg); galt=xgp(parent)%htype(3-jg)
if(m==0) then ! easy
    pass_gamete=gsel
    if(flag==1) call gsel%output(12)
else ! harder
    allocate(xc(m))
    call map_crossings(cl,xc) ! to obtain positions
    if(flag==1) write(12,'(10f8.4)') xc(:)
    call construct_gamete(gsel,galt,xc,goff)
    if(.not.goff%ordered()) call goff_error(goff,gsel,galt,xc,'a')
    pass_gamete=goff
    deallocate(xc)
end if
if(flag==1) call goff%output(12)
end function pass_gamete
!
!=============================================
subroutine construct_gamete(gsel,galt,xc,goff)
!=============================================
type(gamete), intent(inout) :: gsel,galt
real(kind=dp), dimension(:), intent(in) :: xc ! a list of crossovers
type(gamete), intent(inout) :: goff
integer :: jd,joff,jsel,jalt 
real(kind=dp) :: xp
goff%junc(1)%pos=0.0d0
goff%junc(1)%anc=gsel%junc(1)%anc
jd=1; xp=xc(jd) ! jd counts the crossovers
joff=1; jsel=2; jalt=2 ! set up so that jsel and jalt indicate positions ahead of joff
gamete_scan: do
    sel_cycle: do
        if(gsel%junc(jsel)%pos<xp) then ! junction needs to be transferred to offspring
            joff=joff+1
            goff%junc(joff)=gsel%junc(jsel) ! transfers junction
            jsel=jsel+1 ! ready to transfer next junction
            if(jsel>gsel%nj) exit gamete_scan
        else ! a crossover needs to be processed
            exit sel_cycle
        end if
    end do sel_cycle
    alt_cycle: do ! to keep alt in step
        if(galt%junc(jalt)%pos<xp) then ! the junction is bypassed
            jalt=jalt+1
            if(jalt>galt%nj) call goff_error(goff,gsel,galt,xc,'b')
        else ! junction still in play
            exit alt_cycle
        end if
    end do alt_cycle
    ! process the crossover
    joff=joff+1 ! adds to junction count of offspring
    goff%junc(joff)%pos=xc(jd) ! saves position
    goff%junc(joff)%anc=galt%junc(jalt-1)%anc ! picks up the alternative gamete ancestor
    call swap(gsel,galt) 
    call swap(jsel,jalt)
    ! prepare to move on
    if(jd==size(xc)) then ! no more insertions
        xp=maxval(gsel%junc(:)%pos)+0.01d0 ! so that all further junctions on the current gsel are transferred
    else ! more crossovers to process
        jd=jd+1
        xp=xc(jd)
    end if
end do gamete_scan
goff%nj=joff
call goff%reduce()
end subroutine construct_gamete
!
!=============================================
subroutine goff_error(goff,gsel,galt,xc,etype)
!=============================================
type(gamete), intent(in) :: goff,gsel,galt
real(kind=dp), dimension(:), intent(in) :: xc
character(len=*) :: etype
write(11,*) " ... pass_gamete problem _",etype,"_ ... ",size(xc),xc
call gsel%output(11)
call galt%output(11)
call goff%output(11)
stop
end subroutine goff_error
!
!======================
subroutine setup_base()
!======================
integer :: i,j
real(kind=dp) :: u
! parents are already set to 0
forall (i=1:ni) gp(i)%htype(:)%nj=2
forall (i=1:ni) gp(i)%htype(:)%junc(2)%pos=cl ! junc(1) is already set to 0.0d0
if(ibda) then ! only 1 allele is tracked
    forall (i=1:ni) gp(i)%htype(:)%junc(1)%anc=1 
    call random_number(u); i=ceiling(u*real(ni)) ! individual i
    call random_number(u); j=ceiling(u*2.) ! gamete j
    gp(i)%htype(j)%junc(1)%anc=2
else
    j=0
    do i=1,ni
        j=j+1; gp(i)%htype(1)%junc(1)%anc=j
        j=j+1; gp(i)%htype(2)%junc(1)%anc=j
    end do
end if
end subroutine setup_base
!
!=================================
subroutine check_poisson(clen,ich)
!=================================
integer :: ich,i
real(kind=dp) :: clen,y
write(11,'(a/11f7.4)') 'check Poisson ... ',poisson_test(0:10)/sum(poisson_test)
y=exp(-clen)
poisson_test(0)=y
do i=1,9
    y=y*clen/real(i)
    poisson_test(i)=y
end do
poisson_test(10)=1.d0-sum(poisson_test(0:9))
write(11,'(a/11f7.4)') 'true Poisson  ... ',poisson_test(0:10)
end subroutine check_poisson
!
end module utilities