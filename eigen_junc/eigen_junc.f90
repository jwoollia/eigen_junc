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
if(.not.iex) then 
    print *, ' ... no parameter file! ... '; call exit
else
    call date_and_time(time=rtime)
    open(11,file='junc_log_'//rtime(1:6)//'.txt')
    open(14,file='junc_in.txt')
end if
call seed_set(10,11)
!
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
jnull%pos=0.0d0; jnull%anc=0
gnull%nj=0; gnull%junc(:)=jnull
pnull%htype(:)=gnull; pnull%ped(:)=0
!
replicates: do ir=1,nr
    xgp(:)=pnull
    write(arep,'(i4.4)') ir
    open(15,file='junc_rep'//trim(arep)//'_'//rtime(1:6)//'.txt')
    generations: do ig=0,ng
        gp(:)=pnull
        if(ig==0) then ! it is the base
            call setup_base()
        else
            call matings()
        end if
        xgp=gp
        if(max(maxval(gp(:)%htype(1)%nj),maxval(gp(:)%htype(2)%nj))>=size(gnull%junc)) then
            print *, 'junctions per gamete exceeded gamete capacity ... increase capacity'; stop
        end if
        do i=1,ni
            write(15,'(i3,3i4)') ig,i,gp(i)%htype(:)%nj
            call write_gamete(gp(i)%htype(1),15)
            call write_gamete(gp(i)%htype(2),15)
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
