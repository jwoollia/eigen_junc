module parameters
implicit none
integer, parameter, private :: dp=kind(0.0d0)
integer, parameter, public :: maxj=100,maxa=30
integer :: ng,ig,nr,ir ! totals and counters for generations, replicates
integer :: ni,nsel ! number of newborn and number selected, both equally divided among males and females
integer :: ne ! number of eigenvalues
integer :: nn ! number of positions in a net
integer :: isd ! channel number of sdout
integer :: nobj, noba ! empirical maxima
real(kind=dp) :: cl,ptest
real(kind=dp), dimension(:), allocatable :: ew,ea,ev,ep,eb,eave ! eigenfunction stuff
real(kind=dp), dimension(:), allocatable :: xnet,pnet,fnet
real(kind=dp), dimension(:), allocatable :: efc,enet,check_fc
real(kind=dp), dimension(:,:), allocatable :: vfc,vnet
real(kind=dp), dimension(0:10) :: poisson_test ! testing Poisson sampling
character(len=2) :: msys ! mating system
logical :: pg_test ! .true. initiates sampling of constructed gametes for checking
logical :: fc_test ! .true. initiates testing fc construction
logical :: xn_test ! .true. initiates testing xn construction
logical :: ms_test ! .true. initiates testing of Mendelian sampling
!
type junction
    real(kind=dp) :: pos ! location
    integer :: anc ! base generation allele
end type junction
!
type, extends(junction) :: gamete
    type(junction), dimension(maxj) :: jct
    integer :: nj ! number of junctions, >= 2 as it includes the ends
    integer :: na ! number of generation 0 ancestors contributing
    integer, dimension(:), allocatable :: key
    integer, dimension(maxa) :: mix ! ancestor list
    real(kind=dp), dimension(:,:), allocatable :: fc
    integer, dimension(:,:), allocatable :: xn
contains
    procedure, pass :: make => make_gamete          ! initial set up for gamete
    procedure, pass :: wrap => wrap_gamete          ! finalising gamete data
    procedure, pass :: add  => add_junction         ! inserts internal junction before chromosome end
    procedure, pass :: reduce  => reduce_gamete     ! keep junctions as change-points
    procedure, pass :: show    => write_gamete      ! self-explanatory
    procedure, pass :: allele  => get_allele        ! finds ancestor for given position
    procedure, pass :: ordered => check_pos         ! logical check on junction order
    procedure, pass :: profile => profile           ! derives an ibd allele profile by positional net
end type gamete
!
type, extends(gamete) :: info
    integer :: id
    integer, dimension(2) :: ped ! dad 1, mum 2
    type(gamete), dimension(2) :: htype ! gametes of an individual from each parent ped(1) and ped(2)
    ! real(kind=dp) :: bv,ptype
contains
    procedure, pass :: form => form_info            ! initial set up for info
    procedure, pass :: show => write_info
end type info
!
type(gamete) :: gnull
type(info) :: inull
type(info), dimension(:), allocatable, target :: gp,xgp ! holds data on individuals within replicate and generation
!
interface swap
    module procedure sp_swap, dp_swap, int_swap, gam_swap
end interface swap
!
contains   
!
!==================================
integer function get_allele(self,x)     ! allele
!==================================
class(gamete) :: self
real(kind=dp), intent(in) :: x
integer :: k
junctions: do k=2,self%nj
    if(self%jct(k)%pos>x) then
        get_allele=self%jct(k-1)%anc
        exit junctions
    end if
end do junctions
end function get_allele
!
!=============================
subroutine reduce_gamete(self)          ! reduce
!=============================
class(gamete), intent(inout) :: self
integer :: j
j=1 ! preparing to look for redundant junctions
reduce: do 
    if(self%jct(j)%anc==self%jct(j+1)%anc) then ! then j+1 has disappeared
        self%jct(j+1:self%nj-1)=self%jct(j+2:self%nj) ! note j+2 must exist as gamete end is huge() and is constant
        self%jct(self%nj)%pos=0.d0; self%jct(self%nj)%anc=0 ! zero the removed end for clarity
        self%nj=self%nj-1 ! reduce the numbers
    else ! continue along chromosome
        j=j+1 ! move on
        if(j==self%nj) exit reduce
    end if
end do reduce
end subroutine reduce_gamete
!
!===============================
logical function check_pos(self)        ! ordered
!===============================
class(gamete), intent(in) :: self
check_pos=(minval(self%jct(2:self%nj)%pos-self%jct(1:self%nj-1)%pos)>0.d0)
end function check_pos
!
!================================
subroutine write_gamete(self,jch)       ! show
!================================
class(gamete), intent(in) :: self
integer, intent(in), optional :: jch
integer :: i,j,ich
! get channel
if(present(jch)) then
    ich=jch
else
    inquire(file='sdout',number=ich)
end if
! counts
write(ich,'(/2i4)') self%nj,self%na
write(ich,'(30i4)') self%mix(:)
! gametes
i=(self%nj-1)/15
do j=1,i
    write(ich,'(15f8.4)') self%jct(15*j-14:15*j)%pos
    write(ich,'(15(4x,i4))') self%jct(15*j-14:15*j)%anc
end do
write(ich,'(15f8.4)') self%jct(15*i+1:self%nj)%pos
write(ich,'(15(4x,i4))') self%jct(15*i+1:self%nj)%anc
! transforms
i=min(size(self%fc,2),15) ! limit depth of transform output
do j=1,self%na
    write(ich,'(15f8.4)') self%fc(j,1:i)
end do
! nets
i=min(size(self%xn,2),30) ! limit net output
do j=1,self%na
    write(ich,'(30i4)') self%xn(j,1:i)
end do
end subroutine write_gamete
!
!=================================
subroutine add_junction(self,jadd)   ! add
!=================================
class(gamete), intent(inout) :: self
type(junction), intent(in) :: jadd
integer :: jend
jend=self%nj
self%nj=self%nj+1
if(self%nj>maxj) stop '... junctions exceed gamete capacity, stopping: increase maxj!' 
self%jct(self%nj)=self%jct(jend)
self%jct(jend)=jadd
end subroutine add_junction
!
!=========================================
subroutine make_gamete(self,clen,ni,ne,nn)       ! form
!=========================================
class(gamete), intent(inout) :: self
real(kind=dp), intent(in) :: clen
integer, intent(in) :: ni,ne,nn
self%nj=2
self%jct(1)=junction(-clen/2.d0,0)  ! can use to obtain chromosome length from within the gamete
self%jct(2)=junction(clen/2.d0,huge(0))
self%na=0
self%mix(:)=0
allocate(self%key(2*ni),self%fc(size(self%mix),ne),self%xn(size(self%mix),nn))
self%key(:)=0
self%fc(:,:)=0.d0
self%xn(:,:)=0
end subroutine make_gamete
!
!=========================================
subroutine wrap_gamete(self,ew,ep,eb,xnet)       ! wrap
!=========================================
class(gamete), intent(inout) :: self
real(kind=dp), dimension(:), intent(in) :: ew,ep,eb,xnet
integer :: j,ja,jb,je,jn
! list ancestors contributing
self%na=sum(self%key)
if(self%na>size(self%mix)) stop '... ancestors exceed list capacity, stopping: increase maxa!' 
jb=0
do ja=1,size(self%key)
    if(self%key(ja)==0) cycle
    jb=jb+1
    self%mix(jb)=ja
    self%key(ja)=jb
end do
! get fourier coefficients
junctions: do j=1,self%nj-1
    ja=self%jct(j)%anc ! must be in the key
    jb=self%key(ja)
    transform: do je=1,size(ew)
        self%fc(jb,je)=self%fc(jb,je)+eb(je)*sin(ew(je)*self%jct(j+1)%pos+ep(je))
        self%fc(jb,je)=self%fc(jb,je)-eb(je)*sin(ew(je)*self%jct(j)%pos+ep(je))
    end do transform
end do junctions
! get net values
j=1; jn=1
do while ((jn<=size(xnet)).and.(j<self%nj))     
    if(xnet(jn)<self%jct(j)%pos) then    
        ja=self%jct(j-1)%anc
        jb=self%key(ja)
        self%xn(jb,jn)=1
        jn=jn+1 ! move on
    else
        j=j+1 ! move on
        if(j==self%nj) then
            ja=self%jct(j-1)%anc
            jb=self%key(ja)
            self%xn(jb,jn:)=1
        end if
    endif
end do
end subroutine wrap_gamete
!
!==========================
function profile(self,xnet) ! xnet assumed to be mid-points of intervals to avoid ends
!==========================
class(gamete), intent(in) :: self
real(kind=dp), dimension(:), intent(in) :: xnet
real(kind=dp), dimension(size(xnet),size(self%key)) :: profile
integer :: i,j
profile(:,:)=0.d0
i=1; j=2
scan: do
    fill: do while (xnet(i)<self%jct(j)%pos)
        profile(i,self%jct(j-1)%anc)=1.d0
        i=i+1
        if(i>size(xnet)) exit scan
    end do fill
    j=j+1
    if(j>self%nj) exit scan ! should not be needed, but catches bad xnet
end do scan        
end function profile
!
!===============================
subroutine form_info(self,gnull)       ! form
!===============================
class(info), intent(inout) :: self
type(gamete), intent(in) :: gnull
self%htype(1)=gnull
self%htype(2)=gnull
self%id=0
self%ped(:)=0
end subroutine form_info
!
!==============================
subroutine write_info(self,jch)    ! show
!==============================
class(info), intent(in) :: self
integer, intent(in), optional :: jch
integer :: ich
! get channel
if(present(jch)) then
    ich=jch
else
    inquire(file='sdout',number=ich)
end if
write(ich,'(//7i4)') self%id,self%ped
call write_gamete(self%htype(1),ich)
call write_gamete(self%htype(2),ich)
end subroutine write_info
!
!======================
subroutine sp_swap(x,y)
!======================
real, intent(inout) :: x,y
real :: xdum
xdum=x
x=y
y=xdum
end subroutine sp_swap
!
!======================
subroutine dp_swap(x,y)
!======================
real(kind=dp), intent(inout) :: x,y
real(kind=dp) :: xdum
xdum=x
x=y
y=xdum
end subroutine dp_swap
!
!======================
subroutine int_swap(i,j)
!======================
integer, intent(inout) :: i,j
integer :: idum
idum=i
i=j
j=idum
end subroutine int_swap
!
!=========================
subroutine gam_swap(gi,gj)
!=========================
type(gamete), intent(inout) :: gi,gj
type(gamete) :: gdum
gdum=gi
gi=gj
gj=gdum
end subroutine gam_swap
!
end module parameters