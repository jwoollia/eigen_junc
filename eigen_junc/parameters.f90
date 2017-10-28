module parameters
implicit none
integer, parameter, private :: dp=kind(0.0d0)
!
type junction
    real(kind=dp) :: pos ! location
    integer :: anc
end type junction
!
type gamete
    type(junction), dimension(100) :: junc
    integer :: nj ! number of junctions, >= 2 as it includes the ends
contains
    procedure, pass :: reduce => reduce_gamete
    procedure, pass :: output => write_gamete
    procedure, pass :: allele => get_allele
    procedure, pass :: ordered => check_pos
end type gamete
!
type info
    type(gamete), dimension(2) :: htype ! gametes of an individual
    integer, dimension(2) :: ped ! dad 1, mum 2
    ! real(kind=dp) :: bv,ptype
end type info
!
integer :: ng,ig,nr,ir ! totals and counters for generations, replicates
integer :: ni,nsel ! number of newborn and number selected, both equally divided among males and females
type(info), dimension(:), allocatable :: gp,xgp ! holds data on individuals within replicate and generation
real(kind=dp), dimension(0:10) :: poisson_test ! testing Poisson sampling
real(kind=dp) :: cl
character(len=2) :: msys ! mating system
logical :: pg_test ! .true. initiates sampling of constructed gametes for checking
logical :: ibda ! .true. labels only ancestral gamete 1, remainder set to 0 
type(junction) :: jnull 
type(gamete) :: gnull
type(info) :: pnull
!
interface swap
    module procedure sp_swap, dp_swap, int_swap, gam_swap
end interface swap
!
contains
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
!==================================
integer function get_allele(self,x)
!==================================
class(gamete) :: self
real(kind=dp), intent(in) :: x
integer :: k
junctions: do k=2,self%nj
    if(self%junc(k)%pos>x) then
        get_allele=self%junc(k-1)%anc
        exit junctions
    end if
end do junctions
end function get_allele
!
!=============================
subroutine reduce_gamete(self)
!=============================
class(gamete), intent(inout) :: self
integer :: j
j=1 ! preparing to look for redundant junctions
reduce: do
    if(self%junc(j)%anc==self%junc(j+1)%anc) then ! then j+1 has disappeared
        self%junc(j+1:self%nj-1)=self%junc(j+2:self%nj) ! note j+2 must exist as anc=0 for first time only at the end
        self%junc(self%nj)%pos=0.d0; self%junc(self%nj)%anc=0 ! zero the end for clarity
        self%nj=self%nj-1
    else ! continue along chromosome
        j=j+1 ! move on
        if(j==self%nj) exit reduce
    end if
end do reduce
end subroutine reduce_gamete
!
logical function check_pos(self)
class(gamete), intent(in) :: self
check_pos=(minval(self%junc(2:self%nj)%pos-self%junc(1:self%nj-1)%pos)>0.d0)
end function check_pos
!
!================================
subroutine write_gamete(self,ich)
!================================
class(gamete), intent(in) :: self
integer, optional :: ich
integer :: i,j
i=(self%nj-1)/10
do j=1,i
    if(present(ich)) then
        write(ich,'(10f8.4)') self%junc(10*j-9:10*j)%pos
        write(ich,'(10i8)') self%junc(10*j-9:10*j)%anc
    else
        print '(10f8.4)', self%junc(10*j-9:10*j)%pos
        print '(10i8)', self%junc(10*j-9:10*j)%anc  
    end if
end do
if(present(ich)) then
    write(ich,'(10f8.4)') self%junc(10*i+1:self%nj)%pos
    write(ich,'(10i8)') self%junc(10*i+1:self%nj)%anc
else 
    print '(10f8.4)', self%junc(10*i+1:self%nj)%pos
    print '(10i8)', self%junc(10*i+1:self%nj)%anc 
end if
end subroutine write_gamete
!
end module parameters
