module jaw_seeds
implicit none
integer, parameter, private :: dp=kind(0.0d0)
!
contains
!
!=========================    
subroutine seed_set(runit)
!=========================
integer, intent(in), optional :: runit  ! report channel number, so seed can be written into output files
integer, dimension(:), allocatable :: iseed
real(kind=dp), dimension(:), allocatable :: u
integer :: seed_size,sunit,ierr,k
logical :: iexist,iopen
!
call random_seed(size=seed_size) ! find the set up of the random seeds
allocate(iseed(seed_size),u(seed_size)) 
inquire(file='jaw_seed.txt',exist=iexist) ! whether the standard seed file exists
if(iexist) then ! read seed
    open(newunit=sunit,file='jaw_seed.txt',status='old',action='readwrite') ! open the standard unit
    read(sunit,*,iostat=ierr) iseed ! reading the seed
    if(ierr<0) stop 'subroutine seed_set: error reading seed from seed channel' ! stopping if errors
    call random_seed(put=iseed) ! setting the seed
    rewind(sunit) ! preparing for next seed
else ! if there is no file 
    call random_seed(get=iseed) ! get a seed
    open(newunit=sunit,file='jaw_seed.txt',status='new',action='readwrite') ! create the standard unit for next seed
end if
!
if(present(runit)) then ! reporting of current seed is expected
    inquire(runit,opened=iopen) ! whether the reporting unit is open
    if(.not.iopen) then ! not much can be done
        print *, 'subroutine seed_set: report channel number is unopened, seed not reported'
    else ! report the seed
        write(runit,*) '... random number seed ... ',iseed
    end if
end if 
! now get next seed and write to standard file
call random_number(u(:))
iseed(:)=floor(1000000.0d0*u(:))
write(sunit,*) iseed 
! tidy up
close(sunit)
deallocate(iseed)
end subroutine seed_set
!
end module
    