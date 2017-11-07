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
subroutine matings()
!===================
integer :: i
real(kind=dp) :: u
do i=1,ni ! for each newborn individual
    gp(i)%id=i
    select case (msys)
    case('di') ! dioecious
        call random_number(u) 
        gp(i)%ped(1)=ceiling(u*real(ni/2)) ! choose a sire
        call pass_gamete(gp(i)%ped(1),gp(i)%htype(1))
        call random_number(u) ! to choose a dam
        gp(i)%ped(2)=(ni/2)+ceiling(u*real(ni/2)) ! choose a dam
        call pass_gamete(gp(i)%ped(2),gp(i)%htype(2))
    case('ms') ! monoecious selfing ... true random
        call random_number(u) ! to choose a sire
        gp(i)%ped(1)=ceiling(u*real(ni))
        call pass_gamete(gp(i)%ped(1),gp(i)%htype(1))
        call random_number(u) ! to choose a dam
        gp(i)%ped(2)=ceiling(u*real(ni)) ! choose a dam
        call pass_gamete(gp(i)%ped(2),gp(i)%htype(2))
    case('mx') ! monoecious but no selfing
        call random_number(u) ! to choose a sire
        gp(i)%ped(1)=ceiling(u*real(ni))
        call pass_gamete(gp(i)%ped(1),gp(i)%htype(1))
        call random_number(u) ! to choose a dam
        gp(i)%ped(2)=ceiling(u*real(ni-1)) ! choose a dam
        if(gp(i)%ped(2)>=gp(i)%ped(1)) gp(i)%ped(2)=gp(i)%ped(2)+1 ! shift up to avoid and account for selfing
        call pass_gamete(gp(i)%ped(2),gp(i)%htype(2))
    case default
        stop ' ... stopping: error nominating mating system!'
    end select
end do
end subroutine matings
!
!==================================
subroutine pass_gamete(parent,goff)
!==================================
integer, intent(in) :: parent
type(gamete), intent(out) :: goff
type(gamete) :: gsel,galt
integer :: nc,jd,jg,joff,jsel,jalt,jact
real(kind=dp), dimension(:), allocatable :: xc
real(kind=dp) :: u
! set up gamete & recombinations
goff=gnull
nc=crossovers(cl) ! find number of crossovers
poisson_test(min(nc,10))=poisson_test(min(nc,10))+1.0d0 ! accumulate for Poisson check
allocate(xc(nc+1)) 
call map_crossings(cl,xc) ! to obtain positions
! determine starting gamete for x=0
call random_number(u)
jg=1; if(u<0.5d0) jg=2
gsel=xgp(parent)%htype(jg)
galt=xgp(parent)%htype(3-jg)
! initiate goff
goff%jct(1)%anc=gsel%jct(1)%anc
goff%key(gsel%jct(1)%anc)=1
! initiate transfer
jsel=2; jalt=2; jd=1 ! set up so that jsel and jalt indicate positions ahead of joff, jd counts the breaks
gamete_scan: do     
    if(jsel>=gsel%nj.and.jalt>=galt%nj.and.jd>=size(xc)) exit gamete_scan ! only at end the counters will all be at maximum values ... each array has precisely one value >= cl
    jact=minloc((/gsel%jct(jsel)%pos,galt%jct(jalt)%pos,xc(jd)/),dim=1)
    select case (jact)
    case (1)
        call goff%add(gsel%jct(jsel)) ! transfers junction
        goff%key(gsel%jct(jsel)%anc)=1
        if(equal_junction(gsel%jct(jsel),galt%jct(jalt))) jalt=jalt+1 ! keep jalt ahead
        jsel=jsel+1 ! move on to next junction  
    case (2) ! junction in alternative, update gmen, not goff
        if(equal_junction(gsel%jct(jsel),galt%jct(jalt))) jsel=jsel+1 ! keep jsel ahead
        jalt=jalt+1 ! move on to next alternative junction 
    case (3) ! crossover
        call swap(gsel,galt)
        call swap(jsel,jalt)
        ! junction formed only if different, i.e. if(gsel%jct(jsel-1)%anc/=galt%jct(jalt-1)%anc) but this is ignored
        call goff%add( junction(xc(jd),gsel%jct(jsel-1)%anc) )
        goff%key(gsel%jct(jsel-1)%anc)=1
        jd=jd+1
    end select
end do gamete_scan
! tidy & check gamete junctions
call goff%reduce()
if(.not.goff%ordered()) call gam_prob(goff,gsel,galt,xc,'a','y')
! complete gamete
call goff%wrap(ew,ep,eb,xnet)
! monitor
if(pg_test) then
    call random_number(u)
    if(u<ptest) then
        write(12,'(/a)') ' ... new test ... '
        call xgp(parent)%show(12)
        write(12,'(/i4,<nc+1>f8.4)') nc,xc(:)
        call goff%show(12)
    end if
end if
if(fc_test) then
    if(maxval( abs(check_fc(:)-sum(goff%fc(:,:),1)) ) >1.d-6) then
        write(11,'(/a,3i5)') ' ... fc_test fail detected ...'
        call goff%show(11)
    end if
end if
if(xn_test) then
    if( dot_product( (sum(goff%xn(:,:),1)-1),(sum(goff%xn(:,:),1)-1) ) >1.d-6) then
        write(11,'(/a,3i5)') ' ... xn_test fail detected ...'
        call goff%show(11)
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
write(11,'(3a,i4,10f6.4)') " ... pass_gamete problem _",etype,"_ ... ",size(xc),xc
call gsel%show(11)
call galt%show(11)
call goff%show(11)
if(astop=='y'.or.astop=='Y') stop
end subroutine gam_prob
!
!======================
subroutine setup_base()
!======================
integer :: i,j,k
real(kind=dp) :: u
do i=1,ni
    gp(i)%id=i
    do j=1,2
        k=2*(i-1)+j
        gp(i)%htype(j)%jct(1)%anc=k
        gp(i)%htype(j)%key(k)=1
        call gp(i)%htype(j)%wrap(ew,ep,eb,xnet)
    end do
end do
if(fc_test) check_fc(:)=gp(1)%htype(1)%fc(1,:) ! it always exists and is always flat
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
!============================
subroutine check_ms(id,pa,ma)
!============================
type(info), intent(in) :: id,pa,ma
integer, dimension(:), allocatable :: ic
real(kind=dp), dimension(:), allocatable :: xms,yms,xtot,ytot
integer :: ja,jb,ii
allocate(ic(2*ni),xms(ne),xtot(ne),yms(nn),ytot(nn))
xtot(:)=0.d0
ytot(:)=0.d0
ic(:)=0
ic(:)=ic(:)+id%htype(1)%key(:)+id%htype(2)%key(:)
ic(:)=ic(:)+pa%htype(1)%key(:)+pa%htype(2)%key(:)
ic(:)=ic(:)+ma%htype(1)%key(:)+ma%htype(2)%key(:)
do ja=1,2*ni
    if(ic(ja)==0) cycle
    xms(:)=0.d0
    yms(:)=0.d0
    if(id%htype(1)%key(ja)>0) then
        jb=id%htype(1)%key(ja)
        xms(:)=xms(:)+0.5d0*id%htype(1)%fc(jb,:)
        yms(:)=yms(:)+0.5d0*real(id%htype(1)%xn(jb,:))
    end if
    if(id%htype(2)%key(ja)>0) then
        jb=id%htype(2)%key(ja)
        xms(:)=xms(:)+0.5d0*id%htype(2)%fc(jb,:)
        yms(:)=yms(:)+0.5d0*real(id%htype(2)%xn(jb,:))
    end if
    if(pa%htype(1)%key(ja)>0) then
        jb=pa%htype(1)%key(ja)
        xms(:)=xms(:)-0.25d0*pa%htype(1)%fc(jb,:)
        yms(:)=yms(:)-0.25d0*real(pa%htype(1)%xn(jb,:))
    end if
    if(pa%htype(2)%key(ja)>0) then
        jb=pa%htype(2)%key(ja)
        xms(:)=xms(:)-0.25d0*pa%htype(2)%fc(jb,:)
        yms(:)=yms(:)-0.25d0*real(pa%htype(2)%xn(jb,:))
    end if
    if(ma%htype(1)%key(ja)>0) then
        jb=ma%htype(1)%key(ja)
        xms(:)=xms(:)-0.25d0*ma%htype(1)%fc(jb,:)
        yms(:)=yms(:)-0.25d0*real(ma%htype(1)%xn(jb,:))
    end if
    if(ma%htype(2)%key(ja)>0) then
        jb=ma%htype(2)%key(ja)
        xms(:)=xms(:)-0.25d0*ma%htype(2)%fc(jb,:)
        yms(:)=yms(:)-0.25d0*real(ma%htype(2)%xn(jb,:))
    end if  
    write(19,'(/a,i4)') ' ... ancestor ... ',ja
    ii=min(15,size(xms))
    write(19,'(15f8.4)') xms(1:ii)
    ii=min(15,size(yms))
    write(19,'(15f8.4)') yms(1:ii)
    xtot(:)=xtot(:)+xms(:)
    ytot(:)=ytot(:)+yms(:)
end do
if(dot_product(xtot,xtot)>1.d-6) then
    write(19,'(a)') ' ... xtot: dot product fail ... '
    write(11,'(a)') ' ... xtot: dot product fail ... see testing file '
end if
if(dot_product(ytot,ytot)>1.d-6) then
    write(19,'(a)') ' ... ytot: dot product fail ... '
    write(11,'(a)') ' ... ytot: dot product fail ... see testing file '
end if
end subroutine check_ms
!
!====================
subroutine get_msfc()
!====================
real(kind=dp), dimension(:,:), allocatable :: fms
type(info), pointer :: dad=>null(),mum=>null(),me=>null()
type(gamete), dimension(:), pointer :: gdum=>null()
integer :: i,j,k,id
allocate(fms(2*ni,ne))
do id=1,ni
    fms(:,:)=0.d0
    me=>gp(id)
    dad=>xgp(me%ped(1))
    mum=>xgp(me%ped(2))
    !
    gdum=>me%htype
    do i=1,2
        do j=1,gdum(i)%na
            fms(gdum(i)%mix(j),:)=fms(gdum(i)%mix(j),:)+0.5d0*gdum(i)%fc(j,:)
        end do
    end do
    !
    gdum=>dad%htype
    do i=1,2
        do j=1,gdum(i)%na
            fms(gdum(i)%mix(j),:)=fms(gdum(i)%mix(j),:)-0.25d0*gdum(i)%fc(j,:)
        end do
    end do
    !
    gdum=>mum%htype
    do i=1,2
        do j=1,gdum(i)%na
            fms(gdum(i)%mix(j),:)=fms(gdum(i)%mix(j),:)-0.25d0*gdum(i)%fc(j,:)
        end do
    end do
    !       
    efc(:)=efc(:)+sum(fms(:,:),1)/real(ni)
    do i=1,2*ni
        vfc(:,:)=vfc(:,:)+matmul(fms(i:i,:),transpose(fms(i:i,:)))/real(ni)
    end do
end do
!
nullify(dad,mum,me,gdum)
deallocate(fms)
end subroutine get_msfc
!
!=====================
subroutine get_msnet()
!=====================
real(kind=dp), dimension(:,:), allocatable :: fms
type(info), pointer :: dad=>null(),mum=>null(),me=>null()
type(gamete), dimension(:), pointer :: gdum=>null()
integer :: i,j,k,id
allocate(fms(2*ni,nn))
do id=1,ni
    fms(:,:)=0.d0
    me=>gp(id)
    dad=>xgp(me%ped(1))
    mum=>xgp(me%ped(2))
    !
    gdum=>me%htype
    do i=1,2
        j=2; k=1 ! set off to find junctions ahead of the next net assignment 
        do while ((k<=nn).and.(j<gdum(i)%nj))     
            if(xnet(k)<gdum(i)%jct(j)%pos) then
                fms(gdum(i)%jct(j-1)%anc,k)=fms(gdum(i)%jct(j-1)%anc,k)+0.5d0
                k=k+1 ! move on
            else
                j=j+1 ! move on
                if(j==gdum(i)%nj) fms(gdum(i)%jct(j-1)%anc,k:nn)=fms(gdum(i)%jct(j-1)%anc,k:nn)+0.5d0
            endif
        end do
    end do
    !
    gdum=>dad%htype
    do i=1,2
        j=2; k=1 ! set off to find junctions ahead of the next net assignment 
        do while ((k<=nn).and.(j<gdum(i)%nj))      
            if(xnet(k)<gdum(i)%jct(j)%pos) then
                fms(gdum(i)%jct(j-1)%anc,k)=fms(gdum(i)%jct(j-1)%anc,k)-0.25d0
                k=k+1 ! move on
            else
                j=j+1 ! move on
                if(j==gdum(i)%nj) fms(gdum(i)%jct(j-1)%anc,k:nn)=fms(gdum(i)%jct(j-1)%anc,k:nn)-0.25d0
            endif
        end do
    end do
    !
    gdum=>mum%htype
    do i=1,2
        j=2; k=1 ! set off to find junctions ahead of the next net assignment 
        do while ((k<=nn).and.(j<gdum(i)%nj))     
            if(xnet(k)<gdum(i)%jct(j)%pos) then
                fms(gdum(i)%jct(j-1)%anc,k)=fms(gdum(i)%jct(j-1)%anc,k)-0.25d0
                k=k+1 ! move on
            else
                j=j+1 ! move on keeping junctions ahead
                if(j==gdum(i)%nj) fms(gdum(i)%jct(j-1)%anc,k:nn)=fms(gdum(i)%jct(j-1)%anc,k:nn)-0.25d0
            end if
        end do
    end do
    !
    efc(:)=efc(:)+sum(fms(:,:),1)/real(ni)
    do i=1,2*ni
        vfc(:,:)=vfc(:,:)+matmul(fms(i:i,:),transpose(fms(i:i,:)))/real(ni)
    end do
end do
nullify(dad,mum,me,gdum)
deallocate(fms)
end subroutine get_msnet
!
end module utilities