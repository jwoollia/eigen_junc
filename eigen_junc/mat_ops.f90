!     Last change:  JW    7 Mar 2003    4:16 pm
!<<<<<<<<<<<<<
MODULE mat_ops
!<<<<<<<<<<<<<
implicit none
save
INTEGER, PARAMETER, PRIVATE :: dp=KIND(1.0d0)
interface gen_dot_prod ! X' Y Z
  module PROCEDURE sc_vec_mat_vec
  module PROCEDURE sc_vec_mat_mat
  module PROCEDURE sc_mat_mat_vec
  module PROCEDURE sc_mat_mat_mat
end interface
!
INTERFACE quad_prod_sc
  module PROCEDURE sc_mat_vec
  module PROCEDURE sc_mat_mat
END interface
!
interface quad_prod_ar
  module PROCEDURE ar_mat_mat
END interface
!
interface quad_prod
  module PROCEDURE ar_mat_mat
  module PROCEDURE sc_mat_vec
end interface
!
contains
!
!========================
function mat_product(a,c)
!========================
    REAL(KIND=dp), DIMENSION(:) :: a,c
    REAL(KIND=dp), DIMENSION(SIZE(a),SIZE(c)) :: mat_product
    REAL(KIND=dp), DIMENSION(SIZE(a),1) :: duma
    REAL(KIND=dp), DIMENSION(1,SIZE(c)) :: dumc
    duma(:,1)=a
    dumc(1,:)=c
    mat_product=MATMUL(duma,dumc)
  END function mat_product
!===========================================
REAL(KIND=dp) function sc_vec_mat_vec(a,b,c)
!===========================================
    REAL(KIND=dp), DIMENSION(:,:) :: b
    REAL(KIND=dp), DIMENSION(:) :: a,c
    REAL(KIND=dp), DIMENSION(SIZE(b,1)) :: d
      d=MATMUL(b,c)
      sc_vec_mat_vec=SUM(a*d)
  end function sc_vec_mat_vec
!===========================================
REAL(KIND=dp) function sc_vec_mat_mat(a,b,c)
!===========================================
    REAL(KIND=dp), DIMENSION(:,:) :: b,c
    REAL(KIND=dp), DIMENSION(:) :: a
    REAL(KIND=dp), DIMENSION(SIZE(b),1) :: d
      d=MATMUL(b,c)
      sc_vec_mat_mat=SUM(a*d(:,1))
  end function sc_vec_mat_mat
!===========================================
REAL(KIND=dp) function sc_mat_mat_vec(a,b,c)
!===========================================
    REAL(KIND=dp), DIMENSION(:,:) :: b,a
    REAL(KIND=dp), DIMENSION(:) :: c
    REAL(KIND=dp), DIMENSION(1,SIZE(b,2)) :: d
      d=MATMUL(TRANSPOSE(a),b)
      sc_mat_mat_vec=SUM(d(1,:)*c)
  end function sc_mat_mat_vec
!===========================================
REAL(KIND=dp) function sc_mat_mat_mat(a,b,c)
!===========================================
    REAL(KIND=dp), DIMENSION(:,:) :: b,a,c
    REAL(KIND=dp), DIMENSION(SIZE(b,1),1) :: d
      IF(SIZE(a,2)>1) PRINT *, 'gen_dot_prod: ar_mat_mat error'
      d=MATMUL(b,c)
      sc_mat_mat_mat=SUM(a*d)
  end function sc_mat_mat_mat
!=======================
function ar_mat_mat(a,c)
!=======================
    REAL(KIND=dp), DIMENSION(:,:) :: a,c
    REAL(KIND=dp), DIMENSION(SIZE(c,2),SIZE(c,2)) :: ar_mat_mat
    ar_mat_mat=MATMUL(TRANSPOSE(c),MATMUL(a,c))
  END function ar_mat_mat
!=====================================
REAL(KIND=dp) function sc_mat_vec(a,c)
!=====================================
    REAL(KIND=dp), DIMENSION(:,:) :: a
    REAL(KIND=dp), DIMENSION(:) :: c
    REAL(KIND=dp), DIMENSION(SIZE(a,1)) :: d
    IF(SIZE(a,2)/=SIZE(c)) STOP 'sc_mat_vec error'
    d=MATMUL(a,c)
    sc_mat_vec=SUM(c*d)
  END function sc_mat_vec
!=====================================
REAL(KIND=dp) function sc_mat_mat(a,c)
!=====================================
    REAL(KIND=dp), DIMENSION(:,:) :: a,c
    REAL(KIND=dp), DIMENSION(SIZE(a,1),1) :: d
    d=MATMUL(a,c)
    sc_mat_mat=SUM(c*d)
  END function sc_mat_mat
!================
function invrt(b)
! purpose    : invert a non-symmetric matrix a ...  matrix must be non-singular but can be non-positive definite
!              (program stops if singularity is encountered)
! strategy   : uses gauss-jordan algorithm (with partial pivoting)
!              time required is proportional to the cubic power of the order of the matrix ...
!              programmed after stoer,j. and bulirsch,r. 'introduction to numerical analysis',
!              springer verlag 1980, pp. 169-172 ...
!              11/2000 ... jaw f95 conversion maintaining variable names and labels
!              10/2001 ... jaw function
!================
    implicit none
    REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: b ! altered to allow matrix to be unchanged ... jaw 10/2001
    REAL(KIND=dp), DIMENSION(SIZE(b,1),SIZE(b,2)) :: invrt ! added ... jaw 10/2001
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: a
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: vec
    INTEGER, DIMENSION(:), ALLOCATABLE :: iflag
    REAL(KIND=dp) :: zero,diag,off,xx,zz,sve
    INTEGER :: i,j,imax,k,isve,ia,n
    ia=SIZE(b,1); n=SIZE(b,2) ! changed a to b ... jaw 10/2001
    IF(ia/=n) STOP 'invrt: matrix not square'
    ALLOCATE(a(n,n),vec(n),iflag(n)) ! changed b() to a() ... jaw 10/2001
    zero=1.d-12 ! minimum value to be distinguished from 0.d0
    diag=0.d0
    off=0.d0
    !label_1: DO i=1,n
    !  iflag(i)=i
    !END do label_1
    FORALL(i=1:SIZE(iflag)) iflag(i)=i
    a(:,:)=b(:,:) ! carry out operations on copy of b ... jaw 10/2001
    label_2: do i=1,n ! find maximum element in the each column
      xx=abs(a(i,i))
      imax=i
      label_3: do j=i+1,n ! start at i+1 element
        zz=abs(a(j,i))
        if(zz>xx) then
          xx=zz
          imax=j
        end if
      END do label_3
      if(xx<zero) stop 'subroutine invert : matrix is singular'
      if(imax>i) THEN ! interchange row i and row with max element in the column
        !label_4: do k=1,n
        !  sve=a(i,k)
        !  a(i,k)=a(imax,k)
        !  a(imax,k)=sve
        !END do label_4
        vec=a(i,:)
        a(i,:)=a(imax,:)
        a(imax,:)=vec
        isve=iflag(i)
        iflag(i)=iflag(imax)
        iflag(imax)=isve
      end if
      ! transform the matrix
      sve=1.d0/a(i,i)
      a(:,i)=a(:,i)*sve
      a(i,i)=sve
      label_6: do k=1,n
        if(k==i) CYCLE label_6
        label_7: do j=1,n
          if(j/=i) a(j,k)=a(j,k)-a(j,i)*a(i,k)
        end do label_7
        a(i,k)=-a(i,k)*sve
      END DO label_6
    end do label_2
    ! interchange columns (analogous to previous row changes )
    label_8: do i=1,n
      !label_9: do k=1,n
      vec(iflag(:))=a(i,:)
      !END do label_9
      a(i,:)=vec(:)
    end do label_8
    ! checking section has been moved to check_inverse
    invrt=a ! make function value the inverse ... jaw 10/2001
    DEALLOCATE(a,vec,iflag) ! changed b() to a() ... jaw 10/2001
    end function invrt
!======================================
REAL(KIND=dp) FUNCTION check_invrt(a,b)
!======================================
implicit none
REAL(KIND=dp), DIMENSION(:,:) :: a,b
REAL(KIND=dp) :: xx,diag,zero,off
INTEGER :: i,j,n
IF(SIZE(a,2)/=SIZE(b,1).or.SIZE(a,1)/=SIZE(b,2)) PRINT *, 'checking inverses size error'
zero=1.d-12
n=SIZE(b,1)
diag=0._dp
off=0._dp
!multiply matrix with its inverse, check elements
    label_33a: do i=1,n
      label_33b: do j=1,n
        xx=dot_product(a(:,j),b(i,:))
        if(i==j)then
          !if(dabs(xx-1.d0)>zero) print *,i,xx
          diag=diag+dabs(xx-1._dp)
        else
          !if(dabs(xx)>zero) print *,i,j,xx
          off=off+dabs(xx)
        end if
      END do label_33b
    end do label_33a
    xx=diag/REAL(n,KIND=dp)
    !print *,'diagonal : sum =',diag,'     average =',xx
    check_invrt=xx+off/real(n*(n-1),KIND=dp)
    !print *,'off-diag : sum =',off ,'     average =',xx
END function check_invrt
!===========================
SUBROUTINE symmetry(typ,arr)
!===========================
implicit none
REAL(KIND=dp), DIMENSION(:,:), INTENT(INOUT) :: arr
CHARACTER(LEN=2), INTENT(IN) :: typ
INTEGER :: i,j,sz
IF(SIZE(arr,2)/=SIZE(arr,1)) PRINT *, 'symmetry array error ... ',SIZE(arr,1),SIZE(arr,2)
sz=SIZE(arr,1)
select case (typ)
case ('ut','UT','uT','Ut') ! upper triangle given
  FORALL(i=1:sz,j=1:sz,i>j) arr(i,j)=arr(j,i)
case ('lt','LT','lT','Lt') ! lower triangle given
  FORALL(i=1:sz,j=1:sz,i<j) arr(i,j)=arr(j,i)
case default
  PRINT *, 'symmetry option not recognized'
end select
END SUBROUTINE symmetry
!>>>>>>>>>>>>>>>>>
END module mat_ops
!>>>>>>>>>>>>>>>>>
