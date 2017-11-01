!     last change:  jaw   3 aug 2006    9:59 am
module jaw_seeds
implicit none
integer, parameter, private :: dp=kind(0.0d0)
!
contains
!
subroutine seed_set(sunit,runit)
integer, intent(in) :: sunit ! seed unit number, dedicated to seed storage
integer, intent(in), optional :: runit ! report channel number, so seed can be written into output files
real(kind=dp) :: u
integer :: seed_size,ioerr,k
integer, dimension(:), allocatable :: iseed
logical :: iexist,iopen
!
if(present(runit)) then ! check some details
  inquire(runit,opened=iopen)
  if(.not.iopen) print *, 'subroutine seed_set: report channel number is unopened'
  if(runit==sunit) then
    print *, 'subroutine seed_set: seed channel number is same as report channel number'
    write(runit,*) 'subroutine seed_set: seed channel number is same as report channel number'
  end if
end if
!
call random_seed(size=seed_size)
allocate(iseed(seed_size))
inquire(file='jawseed.txt', exist=iexist)
if(iexist) then ! read seed and prepare for next
  open(sunit,file='jawseed.txt',status='old',action='readwrite')
  read(sunit,*,iostat=ioerr) iseed
  if(ioerr<0) then
    print *, 'subroutine seed_set: error reading seed from seed channel'
    stop
  end if
  call random_seed(put=iseed)
  rewind(sunit)
else ! find a seed and open
  call random_seed(get=iseed)
  open(sunit,file='jawseed.txt',status='new',action='readwrite')
end if
if(present(runit)) write(runit,*) '... random number seed ... ',iseed
do k=1,size(iseed) ! for next seed
  call random_number(u)
  iseed(k)=floor(1000000.0d0*u)
end do
write(sunit,*) iseed
close(sunit)
deallocate(iseed)
end subroutine seed_set
end module
