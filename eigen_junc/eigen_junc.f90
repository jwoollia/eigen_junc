program eigen_junc
use parameters
use jaw_seeds
use utilities
use fourier
implicit none
integer, parameter :: dp=kind(0.0d0)
integer :: i,ierr
logical :: iex
character(len=4) :: arep
character(len=10) :: rtime,now
real(kind=dp) :: u
!
inquire(file='junc_in.txt',exist=iex) 
if(.not.iex) stop ' ... no parameter file! '
call date_and_time(time=rtime)
open(11,file='junc_log_'//rtime(1:6)//'.txt')
call seed_set(11)
do i=1,1000 ! burn in
    call random_number(u)
end do
!
open(14,file='junc_in.txt')
read(14,*,iostat=ierr) nr
    write(11,'(a,i8)')   ' ... replicates     ... ',nr
read(14,*,iostat=ierr) ng
    write(11,'(a,i8)')   ' ... generations    ... ',ng
read(14,*,iostat=ierr) ni
    write(11,'(a,i8)')   ' ... individuals    ... ',ni
    allocate(gp(ni),xgp(ni),wh(2*ni))
read(14,*,iostat=ierr) cl
    write(11,'(a,f8.2)') ' ... length (M)     ... ',cl
read(14,*,iostat=ierr) msys
    write(11,'(a,a)')  ' ... mating system  ... ',msys
    if(msys=='di'.and.mod(ni,2)>0) stop ' ... make ni even for the dioecious option'
read(14,*,iostat=ierr) pg_test
    if(pg_test) then
        open(12,file='junc_test_'//rtime(1:6)//'.txt')
        write(11,'(/a)') ' ... gametes are being sampled for checking construction'
        ptest=10./real(2*ni*ng*nr) ! to get ~ 10
        write(11,'(a,f12.8)') ' ... sampling probability ... ',ptest
    end if
read(14,*,iostat=ierr) ibda
    if(ibda) then
        write(11,'(/a)') ' ... tracking single gamete from base'
        wh(:)=0 
        call random_number(u)
        i=ceiling(u*real(2*ni)) ! ibd allele i
        wh(i)=1
        write(11,'(a,i8)') ' ... tracked ibd allele is ... ',i
    end if
read(14,*,iostat=ierr) ne
    write(11,'(/a,i8)')   ' ... eigenfunctions ... ',ne
    if(ne<=0.or.ne>30) stop ' ... # eigenfunctions must be positive and <30 (a setting in type(gamete)'
    allocate(ea(ne),eb(ne),ep(ne),ev(ne),ew(ne),eave(ne))
    write(11,'(a)') ' ... eigenvalues, frequencies, phases, amplitudes and integration amplitudies'
    call basis_f(cl,ea,ew,ev,ep,eb,11)
read(14,*,iostat=ierr) nn
    write(11,'(/a,i8)') ' ... number of internal net positions ... ',nn
    allocate(xnet(nn),pnet(nn),fnet(nn))
    forall(i=1:nn) xnet(i)=cl*( (real(i)-0.5d0)/real(nn)-0.5d0)
    write(11,'(15f8.4)') xnet
!
replicates: do ir=1,nr
    do i=1,ni 
        call xgp(i)%inull(cl) ! set up xgp
    end do
    write(arep,'(i4.4)') ir
    open(15,file='junc_rep'//trim(arep)//'_'//rtime(1:6)//'.txt')
    if(ibda) open(17,file='junc_plot'//trim(arep)//'_'//rtime(1:6)//'.txt')
    generations: do ig=0,ng
        do i=1,ni 
            call gp(i)%inull(cl) ! set up gp
        end do
        if(ibda) then
            pnet(:)=0.d0
            fnet(:)=0.d0
            eave(:)=0.d0
        end if
        if(ig==0) then ! it is the base
            call setup_base()
        else
            call matings()
        end if
        ! report
        write(15,'(/a,i4)') ' ... generation ... ',ig
        do i=1,ni
            call gp(i)%ioutput(15,ibda)
        end do
        if(ibda) then
            fnet(:)=net_value(xnet,eave,ea,ew,ep)
            write(15,'(a)') ' ... pnet ... '
            write(15,'(15f8.4)') pnet(:)
            write(15,'(a)') ' ... fnet ... '
            write(15,'(15f8.4)') fnet(:)
            write(15,'(a)') ' ... eave ... '
            write(15,'(15f8.4)') eave(:)
            write(17,'(/a,i4)') ' ... generation ... ',ig
            write(17,'(3f8.4)') (xnet(i),pnet(i),fnet(i), i=1,size(xnet))
        end if
        ! move on
        xgp=gp
    end do generations
    call date_and_time(time=now)
    print *, now,' ... finished replicate ... ',ir
    write(11,'(/a,a,i5)') now,' ... finished replicate ... ',ir
    close(15)
end do replicates
!
call check_poisson(cl,11)
deallocate(gp,xgp)
read(*,*)
!
end program eigen_junc
