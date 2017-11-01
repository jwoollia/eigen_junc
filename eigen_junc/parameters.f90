module parameters
implicit none
integer, parameter, private :: dp=kind(0.0d0)
integer, parameter, public :: maxj=100
integer :: ng,ig,nr,ir ! totals and counters for generations, replicates
integer :: ni,nsel ! number of newborn and number selected, both equally divided among males and females
integer :: ne ! number of eigenvalues
integer :: nn ! number of positions in a net
integer, dimension(:), allocatable :: wh ! look up array mapping single alleles to ancestor 
real(kind=dp) :: cl,ptest
real(kind=dp), dimension(:), allocatable :: ew,ea,ev,ep,eb,eave ! eigenfunction stuff
real(kind=dp), dimension(:), allocatable :: xnet,pnet,fnet
real(kind=dp), dimension(0:10) :: poisson_test ! testing Poisson sampling
character(len=2) :: msys ! mating system
logical :: pg_test ! .true. initiates sampling of constructed gametes for checking
logical :: ibda ! .true. labels only ancestral gamete 2, remainder set to 1
!
type junction
    real(kind=dp) :: pos ! location
    integer :: anc ! base generation allele
contains
    procedure, pass :: jdefine => put_pos_anc ! puts location and ancestral allele into junction 
end type junction
!
type, extends(junction) :: gamete
    type(junction), dimension(maxj) :: jct
    integer :: nj ! number of junctions, >= 2 as it includes the ends
    logical :: gtype ! .true. if haplotype, .false. if Mendelian sampling
    real(kind=dp) :: clen ! length of chromosome in Morgans, simplifies self%jct(self%nj)%pos
    real(kind=dp), dimension(30) :: fc ! fourier coefficients
contains
    procedure, pass :: gnull    => null_gamete      ! initial set up for gamete
    procedure, pass :: gappend  => append_junction  ! inserts internal junction before chromosome end
    procedure, pass :: greduce  => reduce_gamete    ! keep junctions as change-points
    procedure, pass :: goutput  => write_gamete     ! self-explanatory
    procedure, pass :: gallele  => get_allele       ! finds ancestor gor given position
    procedure, pass :: gordered => check_pos        ! logical check on junction order
    procedure, pass :: gfc      => get_fc           ! obtain fourier coefficients
    procedure, pass :: gprofile => profile          ! derives an ibd allele profile
end type gamete
!
type, extends(gamete) :: info
    type(gamete), dimension(2) :: htype ! gametes of an individual from each parent
    type(gamete), dimension(2) :: vmend ! Mendelian sampling from each parent
    integer :: id
    integer, dimension(2) :: ped ! dad 1, mum 2
    ! real(kind=dp) :: bv,ptype
contains
    procedure, pass :: inull => null_info ! initial set up for info
    procedure, pass :: ioutput => write_info
end type info
!
type(info), dimension(:), allocatable :: gp,xgp ! holds data on individuals within replicate and generation
!
interface swap
    module procedure sp_swap, dp_swap, int_swap, gam_swap
end interface swap
!
contains
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
!==================================
integer function get_allele(self,x)     ! gallele
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
subroutine reduce_gamete(self)          ! greduce
!=============================
class(gamete), intent(inout) :: self
integer :: j
j=1 ! preparing to look for redundant junctions
reduce: do 
    if(self%jct(j)%anc==self%jct(j+1)%anc) then ! then j+1 has disappeared
        self%jct(j+1:self%nj-1)=self%jct(j+2:self%nj) ! note j+2 must exist as gamete end is huge() and is constant
        self%jct(self%nj)%pos=0.d0; self%jct(self%nj)%anc=0 ! zero the end for clarity
        self%nj=self%nj-1
    else ! continue along chromosome
        j=j+1 ! move on
        if(j==self%nj) exit reduce
    end if
end do reduce
end subroutine reduce_gamete
!
!===============================
logical function check_pos(self)        ! gordered
!===============================
class(gamete), intent(in) :: self
check_pos=(minval(self%jct(2:self%nj)%pos-self%jct(1:self%nj-1)%pos)>0.d0)
end function check_pos
!
!================================
subroutine write_gamete(self,ich)       ! goutput
!================================
class(gamete), intent(in) :: self
integer, optional :: ich
integer :: i,j,ie
ie=min(15,size(self%fc))
i=(self%nj-1)/15
do j=1,i
    if(present(ich)) then
        write(ich,'(15f8.4)') self%jct(15*j-14:15*j)%pos
        write(ich,'(15(4x,i4))') self%jct(15*j-14:15*j)%anc
    else
        print '(15f8.4)', self%jct(15*j-14:15*j)%pos
        print '(15(4x,i4))', self%jct(15*j-14:15*j)%anc  
    end if
end do
if(present(ich)) then
    write(ich,'(15f8.4)') self%jct(15*i+1:self%nj)%pos
    write(ich,'(15(4x,i4))') self%jct(15*i+1:self%nj)%anc
else 
    print '(15f8.4)', self%jct(15*i+1:self%nj)%pos
    print '(15(4x,i4))', self%jct(15*i+1:self%nj)%anc  
end if
if(present(ich)) then
    if(ie>1) write(ich,'(15f8.4)') self%fc(1:ie)
else
    if(ie>1) print '(15f8.4)', self%fc(1:ie)
end if
end subroutine write_gamete
!
!====================================
subroutine append_junction(self,jadd)   ! gappend
!====================================
class(gamete), intent(inout) :: self
type(junction), intent(in) :: jadd
integer :: jend
jend=self%nj
self%nj=self%nj+1
if(self%nj>maxj) stop '... junctions exceed gamete capacity, stopping: increase maxp!' 
self%jct(self%nj)=self%jct(jend)
self%jct(jend)=jadd
end subroutine append_junction
!
!======================================
subroutine null_gamete(self,clen,atype)       ! gnull
!======================================
class(gamete), intent(inout) :: self
real(kind=dp), intent(in) :: clen
logical, intent(in) :: atype
self%nj=2
self%clen=clen
call self%jct(1)%jdefine(-clen/2.d0,0)
call self%jct(2)%jdefine(clen/2.d0,huge(0))
self%gtype=atype
self%fc(:)=0.d0
end subroutine null_gamete
!
!==================================
subroutine get_fc(self,ew,ep,eb,wh)           ! gfc
!==================================
class(gamete), intent(inout) :: self 
real(kind=dp), dimension(:), intent(in) :: ew,ep,eb
integer, dimension(:), intent(in) :: wh
integer :: i,j,ky
! self%fc are zero'd already
ky=0
do i=1,min(size(self%fc),size(ew))
    integral: do j=1,self%nj-1
        ky=self%jct(j)%anc
        if(self%gtype) ky=wh(ky)
        if(ky==0) cycle
        self%fc(i)=self%fc(i)+real(ky)*eb(i)*sin(ew(i)*self%jct(j+1)%pos+ep(i))
        self%fc(i)=self%fc(i)-real(ky)*eb(i)*sin(ew(i)*self%jct(j)%pos+ep(i))
    end do integral
end do 
end subroutine get_fc
!
!=============================
function profile(self,xnet,wh) ! xnet assumed to be mid-points of intervals to avoid ends
!=============================
class(gamete), intent(in) :: self
real(kind=dp), dimension(:), intent(in) :: xnet
integer, dimension(:), intent(in) :: wh
real(kind=dp), dimension(size(xnet)) :: profile
integer :: i,j
profile(:)=0.d0
i=1; j=2
scan: do
    fill: do while (xnet(i)<self%jct(j)%pos)
        profile(i)=real(wh(self%jct(j-1)%anc))
        i=i+1
        if(i>size(xnet)) exit scan
    end do fill
    j=j+1
    if(j>self%nj) exit scan ! should not be needed, but catches bad xnet
end do scan        
end function profile
!
!==============================
subroutine null_info(self,clen)       ! inull
!==============================
class(info), intent(inout) :: self
real(kind=dp), intent(in) :: clen
integer :: i
do i=1,2
    call self%htype(i)%gnull(clen,.true.)
    call self%vmend(i)%gnull(clen,.false.)
end do
self%id=0
self%ped(:)=0
end subroutine null_info
!
!===================================
subroutine write_info(self,ich,ibda)    ! ioutput
!===================================
class(info), intent(in) :: self
integer, intent(in) :: ich
logical, intent(in) :: ibda
write(ich,'(7i4)') self%id,self%ped,self%htype(:)%nj,self%vmend(:)%nj
call write_gamete(self%htype(1),ich)
call write_gamete(self%htype(2),ich)
if(ibda) then
    call write_gamete(self%vmend(1),ich)
    call write_gamete(self%vmend(2),ich)
end if
end subroutine write_info
!
!===================================
subroutine put_pos_anc(self,pos,anc)    ! jdefine
!===================================
class(junction), intent(inout) :: self
real(kind=dp) :: pos
integer :: anc
self%pos=pos
self%anc=anc
end subroutine put_pos_anc
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
end module parameters
