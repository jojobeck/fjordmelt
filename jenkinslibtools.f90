module jenkinslibtools
  implicit none
contains
  subroutine change_array_size_2d(xin,N)
    double precision, allocatable, dimension(:,:), intent(inout) ::xin
    double precision, allocatable, dimension(:,:) :: dummy
    integer, intent(in):: N
    if(size(xin(:,4))<N)then
      allocate(dummy(size(xin(:,4)),4))
      ! allocate(dummy(shape(xin))
      dummy=xin
      deallocate(xin)
      allocate(xin(N,4))
      xin=0d0
      xin(1:size(dummy),:)=dummy
      deallocate(dummy)
    end if
  end subroutine change_array_size_2d
  subroutine change_array_size(xin,N)
    double precision, allocatable, dimension(:), intent(inout) ::xin
    double precision, allocatable, dimension(:) :: dummy
    integer, intent(in):: N
    if(size(xin)<N)then
      allocate(dummy(size(xin)))
      dummy=xin
      deallocate(xin)
      allocate(xin(N))
      xin=0d0
      xin(1:size(dummy))=dummy
      deallocate(dummy)
    end if
  end subroutine change_array_size

  ! function cumsum(xin)
  !   double precision,intent(in) :: xin(:)
  !   double precision,dimension((size(xin)+1)) ::cumsum
  !   integer ::i
  !   cumsum=0d0      !starting value is 0
  !   cumsum(2:)=xin
  !   do i=1,size(xin)
  !     cumsum(i+1)=sum(xin(1:i))
  !   end do
  ! end function
  function cumsum(xin)
    double precision,intent(in) :: xin(:)
    double precision,dimension((size(xin)+1)) ::cumsum
    double precision :: c
    integer ::i
    cumsum=0d0      !starting value is 0
    cumsum(2:)=xin
    c=0
    do i=1,size(xin)
      ! cumsum(i+1)=sum(xin(1:i))
      c=c+xin(i)
      cumsum(i+1)=c
    end do
  end function
  ! subroutine cumsum(xin):
  !   double precision,intent(inout),dimension(:):: xin
  !   double precision, allocatable,intent(out) :: x_cum(size(xin)+1)
  !   integer ::i
  !   x_cum(1)=0
  !   x_cum(2:)=x_in
  !   do i=1,size(xin)
  !     x_cum(i+1)=sum(x_in(1:i))
  !   end do
  ! end subroutine

    

  function polyfit(vx, vy, d)
    implicit none
    integer, intent(in)                   :: d
    integer, parameter                    :: dp = selected_real_kind(15, 307)
    real(dp), dimension(d+1)              :: polyfit
    real(dp), dimension(:), intent(in)    :: vx, vy

    real(dp), dimension(:,:), allocatable :: X
    real(dp), dimension(:,:), allocatable :: XT
    real(dp), dimension(:,:), allocatable :: XTX

    integer :: i, j

    integer     :: n, lda, lwork
    integer :: info
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work

    n = d+1
    lda = n
    lwork = n

    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n, size(vx)))
    allocate(X(size(vx), n))
    allocate(XTX(n, n))

    ! prepare the matrix
    do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
    end do

    XT  = transpose(X)
    XTX = matmul(XT, X)

    ! calls to LAPACK subs DGETRF and DGETRI
    call DGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
    call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if

    polyfit = matmul( matmul(XTX, XT), vy)

    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)

  end function

  subroutine polys_fine(x_new,degree,y_new,a)
    integer, intent(in)     :: degree
    integer, parameter      :: dp = selected_real_kind(15, 307)
    double precision,intent(in)     :: x_new(:)
    double precision,intent(out)    :: y_new(size(x_new))
    double precision,intent(in)     :: a(degree+1)
    double precision                :: b
    integer                 :: i,j
    y_new=0d0
    i=1
    do while(i<size(y_new))

      j=0
      b=0d0
      do while(j<degree+1)
       b=b+a(j+1)*(x_new(i)**j)
       j=j+1
      end do
      y_new(i)=a(1)+b 
      i=i+1
    end do

  end subroutine

  subroutine polys(x,y,degree)
    integer, intent(in)     :: degree
    integer, parameter      :: dp = selected_real_kind(15, 307)
    real(dp),intent(in)     :: x(:)
    real(dp),intent(inout)     :: y(:)
    real(dp)                :: y_new(size(y))
    real(dp)                :: a(degree+1),b
    integer                 :: i,j
    a=polyfit(x,y,degree)
    y_new=0d0
    i=1
    do while(i<size(y))

      ! y_new(i)=a(1)+a(2)*x(i)+a(3)*(x(i)**2)+a(4)*(x(i)**3)+a(5)*x(i)**4
      j=0
      b=0d0
      do while(j<degree+1)
       b=b+a(j+1)*(x(i)**j)
       j=j+1
      end do
      y_new(i)=a(1)+b 
      i=i+1
    end do
    y=y_new

  end subroutine

  subroutine make_positive_upwards(x,yy)
    double precision, intent(in)::x(:)
    double precision, intent(inout)         ::yy(:)
    double precision            ::y(size(yy))
    double precision            :: m
    integer::i,j,k
    y=0d0
    y=yy
    do i=1,size(y)-1
      if(y(i)>=y(i+1))then
        k=where_next_bigger(y(i),y(i:))
        if(k/=1)then
          m=(y(k-1+i)-y(i))/(x(k+i-1)-x(i))
          do j=1,k-1
            y(j-1+i)=m*(x(j-1+i)-x(i))+y(i)
          end do
        else
          do j=1,size(y(i:))
            y(i+j)=tan(asin(1e-6))*(x(i+j)+x(i))+y(i)
          end do
        end if
      end if 
    end do
    yy=y
  end subroutine

  Function where_next_bigger(y0,y)
    double precision, intent(in):: y0
    double precision, intent(in):: y(:)
    integer :: where_next_bigger
    integer ::i
    where_next_bigger=1
    do i=1,size(y)
      if(y(i)>y0)then
        where_next_bigger=i
        exit
      end if
    end do
  end Function
end module
