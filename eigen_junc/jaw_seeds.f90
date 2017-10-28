!     Last change:  JAW   3 Aug 2006    9:59 am
module jaw_seeds
implicit none
INTEGER, PARAMETER, PRIVATE :: dp=KIND(0.0d0)
!
contains
!
subroutine seed_set(sunit,runit)
INTEGER, INTENT(IN) :: sunit ! seed unit number, dedicated to seed storage
INTEGER, INTENT(IN), OPTIONAL :: runit ! report channel number, so seed can be written into output files
REAL(KIND=dp) :: u
INTEGER :: seed_size,ioerr,k
INTEGER, DIMENSION(:), ALLOCATABLE :: iseed
LOGICAL :: iexist,iopen
!
IF(PRESENT(runit)) then ! check some details
  INQUIRE(runit,OPENED=iopen)
  IF(.not.iopen) PRINT *, 'SUBROUTINE seed_set: report channel number is unopened'
  IF(runit==sunit) then
    PRINT *, 'SUBROUTINE seed_set: seed channel number is same as report channel number'
    WRITE(runit,*) 'SUBROUTINE seed_set: seed channel number is same as report channel number'
  END if
END if
!
call RANDOM_SEED(SIZE=seed_size)
ALLOCATE(iseed(seed_size))
INQUIRE(FILE='jawseed.txt', EXIST=iexist)
IF(iexist) then ! read seed and prepare for next
  OPEN(sunit,FILE='jawseed.txt',STATUS='OLD',ACTION='READWRITE')
  READ(sunit,*,IOSTAT=ioerr) iseed
  IF(ioerr<0) then
    PRINT *, 'SUBROUTINE seed_set: error reading seed from seed channel'
    STOP
  END if
  call RANDOM_SEED(put=iseed)
  REWIND(sunit)
else ! find a seed and open
  call RANDOM_SEED(get=iseed)
  OPEN(sunit,FILE='jawseed.txt',STATUS='NEW',ACTION='READWRITE')
END if
IF(PRESENT(runit)) WRITE(runit,*) ' ... random number seed ... ',iseed
do k=1,SIZE(iseed) ! for next seed
  call RANDOM_NUMBER(u)
  iseed(k)=FLOOR(1000000.0d0*u)
END do
WRITE(sunit,*) iseed
CLOSE(sunit)
DEALLOCATE(iseed)
END subroutine seed_set
end module
