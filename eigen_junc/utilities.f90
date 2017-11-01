module utilities
use parameters
use sort_kit
use fourier
implicit none
integer, parameter, private :: dp=kind(0.0d0)
!
contains
!
!===================
subroutine matings()            ! change to subroutine
!===================
integer :: i
real(kind=dp) :: u
do i=1,ni ! for each newborn individual
    gp(i)%id=i
    select case (msys)
    case('di') ! dioecious
        call random_number(u) 
        gp(i)%ped(1)=ceiling(u*real(ni/2)) ! choose a sire
        call pass_gamete(gp(i)%ped(1),gp(i)%htype(1),gp(i)%vmend(1))
        call random_number(u) ! to choose a dam
        gp(i)%ped(2)=(ni/2)+ceiling(u*real(ni/2)) ! choose a dam
        call pass_gamete(gp(i)%ped(2),gp(i)%htype(2),gp(i)%vmend(2))
    case('ms') ! monoecious selfing ... true random
        call random_number(u) ! to choose a sire
        gp(i)%ped(1)=ceiling(u*real(ni))
        call pass_gamete(gp(i)%ped(1),gp(i)%htype(1),gp(i)%vmend(1))
        call random_number(u) ! to choose a dam
        gp(i)%ped(2)=ceiling(u*real(ni)) ! choose a dam
        call pass_gamete(gp(i)%ped(2),gp(i)%htype(2),gp(i)%vmend(2))
    case('mx') ! monoecious but no selfing
        call random_number(u) ! to choose a sire
        gp(i)%ped(1)=ceiling(u*real(ni))
        call pass_gamete(gp(i)%ped(1),gp(i)%htype(1),gp(i)%vmend(1))
        call random_number(u) ! to choose a dam
        gp(i)%ped(2)=ceiling(u*real(ni-1)) ! choose a dam
        if(gp(i)%ped(2)>=gp(i)%ped(1)) gp(i)%ped(2)=gp(i)%ped(2)+1 ! shift up to avoid and account for selfing
        call pass_gamete(gp(i)%ped(2),gp(i)%htype(2),gp(i)%vmend(2))
    case default
        stop ' ... stopping: error nominating mating system!'
    end select
end do
end subroutine matings
!
!=======================================
subroutine pass_gamete(parent,goff,gmen)
!=======================================
integer, intent(in) :: parent
type(gamete), intent(out) :: goff,gmen
type(gamete) :: gsel,galt
integer :: nc,jd,jg,joff,jmen,jsel,jalt,jact
real(kind=dp), dimension(:), allocatable :: xc
real(kind=dp) :: u
! set up gametes
call goff%gnull(cl,.true.)
if(ibda) call gmen%gnull(cl,.false.)
! set up recombinations
nc=crossovers(cl) ! find number of crossovers
poisson_test(min(nc,10))=poisson_test(min(nc,10))+1.0d0 ! accumulate for Poisson check
allocate(xc(nc+1)) 
call map_crossings(cl,xc) ! to obtain positions
! determine starting gamete for x=0
call random_number(u)
jg=1; if(u<0.5d0) jg=2
gsel=xgp(parent)%htype(jg); galt=xgp(parent)%htype(3-jg)
! initiate goff & gmen if required
goff%jct(1)%anc=gsel%jct(1)%anc
if(ibda) gmen%jct(1)%anc=wh(gsel%jct(1)%anc)-wh(galt%jct(1)%anc)
! initiate transfer
jsel=2; jalt=2; jd=1 ! set up so that jsel and jalt indicate positions ahead of joff, jd counts the breaks
gamete_scan: do     
    if(jsel>=gsel%nj.and.jalt>=galt%nj.and.jd>=size(xc)) exit gamete_scan ! only at end the counters will all be at maximum values ... each array has precisely one value >= cl
    jact=sum(minloc((/gsel%jct(jsel)%pos,galt%jct(jalt)%pos,xc(jd)/)))
    select case (jact)
    case (1)
        call goff%gappend(gsel%jct(jsel)) ! transfers junction
        if(equal_junction(gsel%jct(jsel),galt%jct(jalt))) jalt=jalt+1 ! keep jalt ahead
        if(ibda) call ms_calc(gmen,gsel%jct(jsel)%pos,gsel%jct(jsel)%anc,galt%jct(jalt-1)%anc) 
        jsel=jsel+1 ! move on to next junction  
    case (2) ! junction in alternative, update gmen, not goff
        if(equal_junction(gsel%jct(jsel),galt%jct(jalt))) jsel=jsel+1 ! keep jsel ahead
        if(ibda) call ms_calc(gmen,galt%jct(jalt)%pos,gsel%jct(jsel-1)%anc,galt%jct(jalt)%anc) 
        jalt=jalt+1 ! move on to next alternative junction 
    case (3) ! crossover
        call swap(gsel,galt)
        call swap(jsel,jalt)
        ! junction formed only if different, i.e. if(gsel%jct(jsel-1)%anc/=galt%jct(jalt-1)%anc) but this is ignored
        call goff%gappend( junction(xc(jd),gsel%jct(jsel-1)%anc) )
        if(ibda) call ms_calc(gmen,xc(jd),gsel%jct(jsel-1)%anc,galt%jct(jalt-1)%anc)  ! note both jsel and jalt are ahead
        jd=jd+1
    end select
end do gamete_scan
call goff%greduce()
if(.not.goff%gordered()) call gam_prob(goff,gsel,galt,xc,'a','y')
call goff%gfc(ew,ep,eb,wh)
pnet(:)=pnet(:)+goff%gprofile(xnet,wh)/real(2*ni)
eave(:)=eave(:)+goff%fc(:)/real(2*ni)
if(ibda) then 
    call gmen%greduce()
    call gmen%gfc(ew,ep,eb,wh)
end if 
! monitoring report    
if(pg_test) then 
    call random_number(u) ! to determine whether this gamete is reported
    if(u<ptest) then
        write(12,'(/a,7i6)') ' ... test ... ',ir,ig,parent,xgp(parent)%htype(:)%nj
        call xgp(parent)%htype(1)%goutput(12)
        call xgp(parent)%htype(2)%goutput(12)
        write(12,'(a,i4/10f8.4)') ' ... crossovers ... ',nc,xc(:)  ! for gamete check
        call goff%goutput(12)
        if(ibda) call gmen%goutput(12)
    end if
end if  
! tidy up
deallocate(xc)   
end subroutine pass_gamete
!
!=================================================
subroutine gam_prob(gsel,galt,goff,xc,etype,astop)
!=================================================
type(gamete), intent(in) :: goff,gsel,galt
real(kind=dp), dimension(:), intent(in) :: xc
character(len=*) :: etype,astop
write(11,'(a,a,a,i4,10f6.4)') " ... pass_gamete problem _",etype,"_ ... ",size(xc),xc
call gsel%goutput(11)
call galt%goutput(11)
call goff%goutput(11)
if(astop=='y'.or.astop=='Y') stop
end subroutine gam_prob
!
!======================
subroutine setup_base()
!======================
integer :: i,j
real(kind=dp) :: u
do i=1,ni
    gp(i)%id=i
    do j=1,2
        gp(i)%htype(j)%jct(1)%anc=2*(i-1)+j
        call gp(i)%htype(j)%gfc(ew,ep,eb,wh)
        pnet(:)=pnet(:)+gp(i)%htype(j)%gprofile(xnet,wh)/real(2*ni)
        eave(:)=eave(:)+gp(i)%htype(j)%fc(:)/real(2*ni)
    end do
end do
end subroutine setup_base
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
  if(x>u) exit
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
if(size(xc)>1) then
    call random_number(u); xc(1)=u
    next: do j=2,size(xc)-1
        call random_number(u)
        sort: do i=1,j-1 ! simple sort as array is short
            if(u<xc(i)) then
                xc(i+1:j)=xc(i:j-1)
                xc(i)=u
                cycle next
            end if
            xc(j)=u
        end do sort
    end do next
end if
xc(size(xc))=1.1d0 ! add a long stop
xc(:)=xc(:)-0.5d0 ! centre the distribution
xc(:)=clen*xc(:) ! scale to chromosome length
end subroutine map_crossings   
! 
!=================================
subroutine check_poisson(clen,ich)
!=================================
integer, intent(in) :: ich
real(kind=dp), intent(in) :: clen
integer :: i
real(kind=dp) :: y
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
!=================================
logical function equal_gamete(a,b)
!=================================
type(gamete), intent(in) :: a,b
integer :: i
if(a%nj/=b%nj) then
    equal_gamete=.false.
else
    equal_gamete=.true.
    junctions: do i=1,a%nj
        if(equal_junction(a%jct(i),b%jct(i))) cycle
        equal_gamete=.false.
        exit junctions
    end do junctions
end if
end function equal_gamete
!
!===================================
logical function equal_junction(a,b)
!===================================
type(junction), intent(in) :: a,b
if((a%anc==b%anc).and.(abs(a%pos-b%pos)<1.0d-6)) then
    equal_junction=.true.
else
    equal_junction=.false.
end if
end function equal_junction
!
!==============================
subroutine ms_calc(gm,xp,js,ja)
!==============================
type(gamete), intent(inout) :: gm
real(kind=dp) :: xp
integer, intent(in) :: js,ja
integer :: ms
ms=wh(js)-wh(ja)
call gm%gappend( junction(xp,ms) )
end subroutine ms_calc
!
end module utilities