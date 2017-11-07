program eigen_junc
use parameters
use jaw_seeds
use utilities
use fourier
implicit none
integer, parameter :: dp=kind(0.0d0)
integer :: i,ierr,jp,jm
logical :: iex
character(len=4) :: arep
character(len=10) :: rtime,now
real(kind=dp) :: u
!
inquire(file='sdout',number=isd)
inquire(file='junc_in.txt',exist=iex) 
if(.not.iex) stop ' ... no parameter file! '
call date_and_time(time=rtime)
open(11,file='j'//rtime(1:6)//'_log.txt')
call seed_set(11)
do i=1,1000 ! burn in
    call random_number(u)
end do
open(14,file='junc_in.txt')
! 11 ... log file
! 12 ... testing gametes
! 13 ... unlucky
! 14 ... input file
! 15 ... output per replicate
! 16 ... pedigree per replicate
! 17 ... positional net output for plotting per replicate
! 18 ... means and variances of ms
! 19 ... testing mendelian sampling
read(14,*,iostat=ierr) nr
    write(11,'(a,i8)')   ' ... replicates     ... ',nr
read(14,*,iostat=ierr) ng
    write(11,'(a,i8)')   ' ... generations    ... ',ng
read(14,*,iostat=ierr) ni
    write(11,'(a,i8)')   ' ... individuals    ... ',ni
    if(ni>50) stop ' ... increase parameter maxa'
    allocate(gp(ni),xgp(ni))
read(14,*,iostat=ierr) cl
    write(11,'(a,f8.2)') ' ... length (M)     ... ',cl
read(14,*,iostat=ierr) msys
    write(11,'(a,a)')  ' ... mating system  ... ',msys
    if(msys=='di'.and.mod(ni,2)>0) stop ' ... make ni even for the dioecious option'
read(14,*,iostat=ierr) pg_test
    if(pg_test) then
        open(12,file='j'//rtime(1:6)//'_gtest.txt')
        write(11,'(/a)') ' ... gametes are being sampled for checking construction'
        ptest=10./real(2*ni*ng*nr) ! to get ~ 10
        write(11,'(a,f12.8)') ' ... sampling probability ... ',ptest
    end if
read(14,*,iostat=ierr) ne
    write(11,'(/a,i8)')   ' ... eigenfunctions ... ',ne
    if(ne<=0) stop ' ... # eigenfunctions must be positive'
    allocate(ea(ne),eb(ne),ep(ne),ev(ne),ew(ne),eave(ne))
    write(11,'(a)') ' ... eigenvalues, frequencies, phases, amplitudes and integration amplitudies'
    call basis_f(cl,ea,ew,ev,ep,eb,11)
    open(18,file='j'//rtime(1:6)//'_ms.txt')
    allocate(efc(ne),vfc(ne,ne))
read(14,*,iostat=ierr) nn
    write(11,'(/a,i8)') ' ... number of internal net positions ... ',nn
    allocate(xnet(nn),enet(nn),vnet(nn,nn))
    forall(i=1:nn) xnet(i)=cl*((real(i)-0.5d0)/real(nn)-0.5d0)
    write(11,'(15f8.4)') xnet
read(14,*,iostat=ierr) fc_test
    if(fc_test) then
        write(11,'(/a)') ' ... fourier coefficients are checked after gamete construction'
        allocate(check_fc(ne))
    end if
read(14,*,iostat=ierr) xn_test
    if(xn_test) then
        write(11,'(/a)') ' ... net values are checked after gamete construction'
    end if
read(14,*,iostat=ierr) ms_test
    if(ms_test) then
        open(19,file='j'//rtime(1:6)//'_mtest.txt')
        write(11,'(/a)') ' ... Mendelian sampling is being checked'
        write(11,'(a,f12.8)') ' ... sampling probability ... ',ptest
    end if
!
call gnull%make(cl,ni,ne,nn)
call inull%form(gnull)
nobj=0
noba=0
!
replicates: do ir=1,nr
    write(arep,'(i4.4)') ir
    open(15,file='j'//rtime(1:6)//'_'//trim(arep)//'_rep.txt')
    open(16,file='j'//rtime(1:6)//'_'//trim(arep)//'_ped.txt')
    open(17,file='j'//rtime(1:6)//'_'//trim(arep)//'_plot.txt')
    efc(:)=0.d0
    vfc(:,:)=0.d0
    xgp(:)=inull  ! set up xgp
    generations: do ig=0,ng
        gp(:)=inull
        if(ig==0) then ! it is the base
            call setup_base()
        else
            call matings()
        end if
        ! report & monitor
        write(15,'(//a,i4//)') ' ... generation ... ',ig
        do i=1,ni
            call gp(i)%show(15) 
            write(16,'(i2,i3.3,2(1x,i2,i3.3))') ig,gp(i)%id,max(0,ig-1),gp(i)%ped(1),max(0,ig-1),gp(i)%ped(1)
            if(ms_test.and.ig>0) then
                call random_number(u)
                if(u<ptest) then
                    write(19,'(/a)') ' ... new test ... new test ... new test ... '
                    call gp(i)%show(19)
                    jp=gp(i)%ped(1); jm=gp(i)%ped(2)
                    call xgp(jp)%show(19)
                    call xgp(jm)%show(19)
                    call check_ms(gp(i),xgp(jp),xgp(jm))
                end if
            end if
        end do
        if(ig>0) then
            call get_msfc()
            write(18,'(2a,2i4,15f8.5)') 'f','m',ir,ig,efc(:)
            write(18,'(2a,2i4,15f8.5,14(/9x,15f8.4))') 'f','v',ir,ig,vfc(:,:)
            call get_msnet()
            write(18,'(2a,2i4,15f8.5)') 'n','m',ir,ig,enet(:)
            write(18,'(2a,2i4,15f8.5,14(/9x,15f8.4))') 'n','v',ir,ig,vnet(:,:)
        end if
        nobj=max(nobj,maxval(gp(:)%htype(1)%nj),maxval(gp(:)%htype(2)%nj))
        noba=max(noba,maxval(gp(:)%htype(1)%na),maxval(gp(:)%htype(2)%na))
        ! move on
        xgp(:)=gp(:)
    end do generations
    call date_and_time(time=now)
    print *, now,' ... finished replicate ... ',ir
    write(11,'(/2a,i5)') now,' ... finished replicate ... ',ir
    close(15)
end do replicates
!
write(11,'(a,2i4)') ' ... maximum number of junctions and ancestors observed ... ',nobj,noba
call check_poisson(cl,11)
deallocate(gp,xgp)
read(*,*)
!
end program eigen_junc
