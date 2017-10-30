program eigen_junc
use parameters
use jaw_seeds
use utilities
implicit none
integer, parameter :: dp=kind(0.0d0)
integer :: i,ierr
logical :: iex
character(len=4) :: arep
character(len=10) :: rtime,now
!
inquire(file='junc_in.txt',exist=iex); 
if(.not.iex) stop ' ... no parameter file! '
call date_and_time(time=rtime)
open(11,file='junc_log_'//rtime(1:6)//'.txt')
call seed_set(10,11)
!
open(14,file='junc_in.txt')
read(14,*,iostat=ierr) nr; write(11,'(a,i8)')  ' ... replicates    ... ',nr
read(14,*,iostat=ierr) ng; write(11,'(a,i8)')  ' ... generations   ... ',ng
read(14,*,iostat=ierr) ni; write(11,'(a,i8)')  ' ... individuals   ... ',ni
read(14,*,iostat=ierr) cl; write(11,'(a,f6.2)')  ' ... length (M)    ... ',cl
read(14,*,iostat=ierr) msys; write(11,'(a,a)') ' ... mating system ... ',msys
if(msys=='di'.and.mod(ni,2)>0) stop ' ... make ni even for the dioecious option'
read(14,*,iostat=ierr) pg_test
if(pg_test) then
    open(12,file='junc_test_'//rtime(1:6)//'.txt')
    write(11,'(a)') ' ... gametes are being sampled for checking construction'
end if
read(14,*,iostat=ierr) ibda
if(ibda) write(11,'(a)') ' ... tracking single gamete from base'
!
allocate(gp(ni),xgp(ni))
if(ibda) allocate(wh(2*ni))
replicates: do ir=1,nr
    do i=1,ni 
        call xgp(i)%inull(cl) ! set up xgp
    end do
    write(arep,'(i4.4)') ir
    open(15,file='junc_rep'//trim(arep)//'_'//rtime(1:6)//'.txt')
    generations: do ig=0,ng
        do i=1,ni 
            call gp(i)%inull(cl) ! set up gp
        end do
        if(ig==0) then ! it is the base
            call setup_base()
        else
            call matings()
        end if
        xgp=gp
        ! report
        write(15,'(a,i4)') ' ... generation ... ',ig
        do i=1,ni
            call gp(i)%ioutput(15,ibda)
        end do
    end do generations
    call date_and_time(time=now)
    print *, now,' ... finished replicate ... ',ir
    write(11,'(a,a,i5)') now,' ... finished replicate ... ',ir
    close(15)
end do replicates
!
call check_poisson(cl,11)
deallocate(gp,xgp)
read(*,*)
!
end program eigen_junc
