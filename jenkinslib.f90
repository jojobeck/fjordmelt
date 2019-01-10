
module jenkinslib

  use lib_array, only: interp1d, is_nan,linspace
  use jenkinslibtools, only: change_array_size, change_array_size_2d, &
    make_positive_upwards,where_next_bigger,polys,polys_fine,polyfit

  implicit None
  private
  public :: melting_process, ode, cal_Sb_melt,temperature_profile,interpolate,Ujenkins,&
    where_next_bigger,submarine_melt,submarine_melt_test,&
    floating_tongue_profile,vertical_tongue_profile,compute_oceanographich_profile,init_y0, &
    aggregate_submarine_melt_floating,aggregate_submarine_melt_angle,Ujenkins_cone,&
    floating_tongue_profile_half,aggregate_submarine_melt_angle_half&
    ,floating_tongue_profile_half_test,cumsum

  ! integer, parameter :: dp = kind(0.d0)

contains
  subroutine submarine_melt(glacier_x,glacier_dx,glacier_hb,glacier_zb,method,discharge, fjordmelt)
    double precision,dimension(:),intent(in)      :: glacier_x
    double precision,dimension(:),intent(in)      :: glacier_dx
    double precision,dimension(:),intent(in)      :: glacier_zb
    double precision,dimension(:),intent(in)      :: glacier_hb
    double precision,intent(out)                  :: fjordmelt(size(glacier_hb))
    double precision                              :: glacier_x_right(size(glacier_x))
    double precision                              :: x_half(2*size(glacier_x)-1),dx_half(2*size(glacier_x)-1)
    character(20),intent(in) :: method
    real,optional,intent(in) :: discharge
    double precision                 :: DU0,hbgl,ds,E0
    character(20)            :: glacier_type        ! tidewater or floating tongue
    double precision,dimension(:),allocatable                 :: s,x,hb,sina,melt,Ta,Sa !profile for jenkins calculation on fine grid ds
    double precision                              :: y0(4),yout(4)      !staring vlaues of plumer for velocity flux andT,S
    integer :: maxlen,n_gl,real_len

    n_gl=size(glacier_x)
    ds= 50d0                     !fine grid for jenkins mode: resolution  1 m


    DU0=1d0              !minimal discharge value
    E0=0.036
    hbgl=0.5*(glacier_zb(1)+glacier_zb(2))  !true depth of shelf
    glacier_x_right=glacier_x+0.5*glacier_dx
   
 
    
    !TODO: optional doesn't work when nothing is asigne
    !parameters
    glacier_type='tidewater'    ! default type
    if(size(glacier_hb)>1)glacier_type='floating'
    if(present(discharge))DU0=discharge
    

    !compute fine jenkins profile including interpolation
    maxlen=int((glacier_x(n_gl)+glacier_dx(n_gl)-glacier_x(1)+abs(glacier_hb(1)))/ds+1)
    allocate(s(maxlen))
    allocate(sina(maxlen))
    allocate(hb(maxlen))
    allocate(x(maxlen))
    s=-9999d0  !missing values
    sina=-9999d0
    x=-9999d0
    hb=-9999d0
    if(glacier_type=='floating')then
      !calculate floating tongue profile
      ! call floating_tongue_profile(ds,glacier_x,glacier_dx,glacier_hb,hbgl,size(glacier_x),s,x,hb,sina,real_len)
      call floating_tongue_profile_half(ds,glacier_x,glacier_dx,glacier_hb,hbgl,size(glacier_x),s,x,hb,sina,real_len,x_half,dx_half)
    else
      call vertical_tongue_profile(ds,hbgl,s,hb,sina,real_len)
    endif
    !calculate T/S profile
    allocate(Ta(real_len))
    allocate(Sa(real_len))
    call compute_oceanographich_profile(hb(1:real_len),Ta,Sa)  !ambient Temperature and salinity

    !melting:
    allocate(melt(real_len))
    call init_y0(Ta(1),Sa(1),sina(1),E0,DU0,y0,method)         !starting value for plue
    melt=0d0
    call melting_process(s(1:real_len),hb(1:real_len),Ta,Sa,E0,sina,melt,y0,real_len,method,yout)

    ! transform back of melting
    fjordmelt=0d0
    if(glacier_type=='floating')then  ! need to transform in horizontal and vertial melting

      ! call aggregate_submarine_melt_floating(ds,x(1:real_len),hb(1:real_len),glacier_x_right,glacier_dx,melt,fjordmelt)   ! important: use x of cell edge
      call aggregate_submarine_melt_angle(ds,x(1:real_len),hb(1:real_len),glacier_x_right,glacier_dx,melt,fjordmelt)   ! important: use x of cell edge
      write(*,*)'fjordmelt',fjordmelt
      fjordmelt=0d0
      call aggregate_submarine_melt_angle_half(ds,x(1:real_len),hb(1:real_len),glacier_x_right, &
                                                  glacier_dx,x_half,dx_half,melt,fjordmelt)   ! halfstep of glacier_x
      write(*,*)'fjordmelt half',fjordmelt

    else      !all vertical melting
      fjordmelt=sum(melt)*ds/glacier_dx(1)          !todo think about horizontal and vertical

    endif
    ! call aggregate_submarine_melt(s,x,hb,glacier_x_right,melt,glacier_type,fjordmelt)   ! important: use x of cell edge
    ! sum s til x > x_right


  end subroutine

  subroutine  aggregate_submarine_melt_floating(ds,x,hb,glacier_x,glacier_dx,melt,glacier_melt)
    !melting from gl+1:c
    double precision, intent(in) :: ds,x(:),hb(:),melt(:)
    double precision, intent(in) :: glacier_x(:),glacier_dx(:)
    double precision,intent(out) :: glacier_melt(size(glacier_x))
    double precision             :: melt_horizontal(size(glacier_x))   !horizotal with x
    double precision             :: m,hb0
    integer                      :: i,j
    melt_horizontal=0d0
    i=1
    j=2
    m=0d0
    hb0=hb(1)
    do while(i<=size(x) .and. j<=size(glacier_x))
      if(x(i)>glacier_x(j)) then             !aggregate in one glacier_melt cell
        melt_horizontal(j)=m*ds/glacier_dx(j)  !projection to x
        m=0d0
        hb0=hb(i)
        j=j+1
        if(melt(i+1)==0d0)then              ! skip to aggregate 0 melt
          j=size(glacier_x)+1
        end if
      end if
      m=m+melt(i)
      if(x(i)> glacier_x(size(glacier_x)-1))then         !skip to run loop through last glacier cell
        melt_horizontal(size(glacier_x))=sum(melt(i:))*ds/glacier_dx(size(glacier_x))
        j=size(glacier_x)+1
      end if
      i=i+1
    end do
    glacier_melt=melt_horizontal      
  end subroutine

  subroutine  aggregate_submarine_melt_angle_half(ds,x,hb,glacier_x,glacier_dx,x_coarse,dx_coarse,melt,glacier_melt)
    !going in half step size of glacier grid- massconservation meling to left cell for slope!
    !melt from gl:c
    double precision, intent(in) :: ds,x(:),hb(:),melt(:)
    double precision, intent(in) :: glacier_x(:),glacier_dx(:)
    double precision, intent(in) ::x_coarse(:),dx_coarse(:)
    double precision,intent(out) :: glacier_melt(size(glacier_x))
    double precision             :: melt_horizontal(size(x_coarse)), melt_vertical(size(x_coarse))
    double precision             :: m,hb0,m_s,dhb,ds_glacier
    integer                      :: i,j,gll
    gll=size(x_coarse)
    melt_vertical=0d0
    melt_horizontal=0d0
    i=1
    j=2
    m=0d0
    hb0=hb(1)
    do while(i<=size(x) .and. j<=gll)
      if(x(i)>x_coarse(j)) then             !aggregate in one glacier_melt cell
        dhb=hb(i-1)-hb0                         !glacier dhb
        ds_glacier=( dhb*dhb +dx_coarse(j)*dx_coarse(j) )**0.5  !on glacier_grid
        m_s=m*ds/ds_glacier                  !average melting on ds_glacier
        if(m_s*dhb/ds_glacier>x_coarse(j))then
          write(*,*) 'Warning to big timestep'
        end if
        melt_vertical(j-1)= 0.5*m_s*m_s* dhb/dx_coarse(j)  !melt area of triangel to neighbour left cell
        melt_horizontal(j)=m*ds*(1-melt_vertical(j-1))       !cumulative melt - left cell melting
        !intial values for next j+1 cell:
        m=0d0
        hb0=hb(i-1)
        j=j+1
        if(melt(i+1)==0d0)then              ! skip to aggregate 0 melt
          j=size(x_coarse)+1
        end if
      end if
      m=m+melt(i)
      if(x(i)>x_coarse(gll-1))then         !skip to run loop through last glacier cell
        dhb=abs(hb(i) )                        !glacier ghb to zero
        ds_glacier=(dhb*dhb +dx_coarse(gll)*dx_coarse(gll) )**0.5  !on glacier_grid
        m_s=sum(melt(i:))*ds/ds_glacier                  !average melting on ds_glacier
        melt_vertical(gll-1)=0.5*m_s*m_s*dhb/dx_coarse(gll)   !melt area of triangel to left call
        melt_horizontal(gll)=m_s*ds_glacier-melt_vertical(gll-1)
        ! melt_vertical(size(glacier_x)-1)=sum(melt(i:))*ds/abs( hb(size(hb)) - hb(1))
        ! j=j+1
        j=size(x_coarse)+1
      end if
      i=i+1
    end do
    i=1
    melt_horizontal=melt_horizontal+melt_vertical    !add both melts toegather
    glacier_melt(1)=melt_horizontal(1)
    do while (i<=size(glacier_x))                    !bring melt of coarse_x back to glacier grind
      glacier_melt(i+1)=melt_horizontal(2*i) +melt_horizontal(2*i+1)
      i=i+1
    end do
    glacier_melt=(melt_vertical+melt_horizontal)/glacier_dx      !average vertical/horitonal melt @ gl/gl+1,tongue has melting @gl!!!!
  end subroutine
  subroutine  aggregate_submarine_melt_angle(ds,x,hb,glacier_x,glacier_dx,melt,glacier_melt)
    double precision, intent(in) :: ds,x(:),hb(:),melt(:)
    double precision, intent(in) :: glacier_x(:),glacier_dx(:)
    double precision,intent(out) :: glacier_melt(size(glacier_x))
    ! double precision             :: melt_horizontal(size(glacier_x))   !horizotal with x
    double precision             :: melt_horizontal(size(glacier_x)), melt_vertical(size(glacier_x))
    double precision             :: m,hb0,m_s,dhb,ds_glacier
    integer                      :: i,j,gll
    gll=size(glacier_x)
    melt_vertical=0d0
    melt_horizontal=0d0
    i=1
    j=2
    m=0d0
    hb0=hb(1)
    do while(i<=size(x) .and. j<=gll)
      if(x(i)>glacier_x(j)) then             !aggregate in one glacier_melt cell
        dhb=hb(i-1)-hb0                         !glacier dhb
        ds_glacier=( dhb*dhb + glacier_dx(j)*glacier_dx(j) )**0.5  !on glacier_grid
        m_s=m*ds/ds_glacier                  !average melting on ds_glacier
        if(m_s*dhb/ds_glacier>glacier_x(j))then
          write(*,*) 'Warning to bis timestep'
        end if
        melt_vertical(j-1)= 0.5*m_s*m_s* dhb/glacier_dx(j)  !melt area of triangel to neighbour left cell
        melt_horizontal(j)=m*ds*(1-melt_vertical(j-1))       !cumulative melt - left cell melting
        !intial values for next j+1 cell:
        m=0d0
        hb0=hb(i-1)
        j=j+1
        if(melt(i+1)==0d0)then              ! skip to aggregate 0 melt
          j=size(glacier_x)+1
        end if
      end if

      m=m+melt(i)                            !cumualtive melt 

      if(x(i)> glacier_x(gll-1))then         !skip to run loop through last glacier cell
        dhb=abs(hb(i) )                        !glacier ghb to zero
        ds_glacier=(dhb*dhb + glacier_dx(gll)*glacier_dx(gll) )**0.5  !on glacier_grid
        m_s=sum(melt(i:))*ds/ds_glacier                  !average melting on ds_glacier
        melt_vertical(gll-1)=0.5*m_s*m_s*dhb/glacier_dx(gll)   !melt area of triangel to left call
        melt_horizontal(gll)=m_s*ds_glacier-melt_vertical(gll-1)
        ! melt_vertical(size(glacier_x)-1)=sum(melt(i:))*ds/abs( hb(size(hb)) - hb(1))
        ! j=j+1
        j=size(glacier_x)+1
      end if
      i=i+1
    end do
    glacier_melt=(melt_vertical+melt_horizontal)/glacier_dx      !average vertical/horitonal melt @ gl/gl+1,tongue has melting @gl!!!!
  end subroutine
  
  subroutine compute_oceanographich_profile(hb,Ta_in,Sa_in)
    double precision,intent(in)     ::hb(:)
    double precision,intent(inout)     ::Ta_in(:)
    double precision,intent(inout)     ::Sa_in(:)
    integer     ::steps
    double precision :: Ta_bottom,Ta_surface,d1,d2,Sa_bottom,Sa_surface
    double precision :: T(size(hb)),Z(size(hb)),S(size(hb))
    steps=size(hb)
    !++++dummy beracuse it will be a parameter TODO:
    Ta_bottom=4d0
    Ta_surface=4d0
    d1=-300d0
    d2=-100d0
    Sa_bottom=34.75d0
    Sa_surface=34.d0
    call temperature_profile(Ta_bottom,Ta_surface,d1,d2,hb(1),steps,T,Z)
    call temperature_profile(Sa_bottom,Sa_surface,d1,d2,hb(1),steps,S,Z)
    Sa_in=interpolate(Z,S,hb)
    Ta_in=interpolate(Z,T,hb)

  end subroutine
  subroutine floating_tongue_profile_half_test(ds,glacier_x,glacier_dx,glacier_hb,hbgl,n_gl,s_in,x_in,hb_in,sina_in,&
                                                              real_len,x_coarse,dx_coarse_back,x_learn,hb_learn)
    double precision,intent(in)    :: glacier_x(:)
    double precision,intent(in)    :: glacier_dx(:)
    double precision,intent(in)    :: glacier_hb(:)
    double precision,intent(inout)    :: s_in(:),x_in(:),hb_in(:),sina_in(:)
    integer,intent(inout)    :: real_len
    double precision, intent(in)   :: hbgl
    double precision,intent(in)    :: ds
    integer,intent(in)    :: n_gl    !=size(glacier_x))
    double precision               :: s_coarse(2*n_gl-1)
    double precision, intent(inout)               :: x_coarse(2*n_gl-1),dx_coarse_back(2*n_gl-1),hb_learn(n_gl+1),x_learn(n_gl+1)
    double precision               :: hb_coarse(2*n_gl-1),ds_coarse(2*n_gl-2),dhb_coarse(2*n_gl-2),dx_coarse(2*n_gl-2)
    double precision,dimension(:) ,allocatable                :: s,x,hb,sina,a,x_reg_dx,hb_reg_dx,s_reg_dx,ds_reg_dx,dhb_reg_dx!profile for jenkins calculation on fine grid ds
  
    double precision                :: dx_reg
    integer                                     :: end_coarse,steps_ds,i

    ! true floating part on coarse glacier grid:
  
    end_coarse=2*n_gl-1
    !for polynomial fit
    hb_learn(1:n_gl)=glacier_hb
    x_learn(1:n_gl)=glacier_x
    x_learn(1)=glacier_x(1)+0.5*glacier_dx(1)
    x_learn(n_gl+1)=glacier_x(n_gl)+0.5*glacier_dx(n_gl)
    x_coarse(1)=x_learn(1)
    hb_learn(1)=hbgl
    hb_learn(n_gl+1)=0
    ! dx_coarse=0.5d0*glacier_dx
    i=1
    do while (i<=end_coarse/2)
      x_coarse(2*i)=glacier_x(i+1)
      x_coarse(2*i+1)=glacier_x(i+1)+0.5*glacier_dx(i+1)
      i=i+1
    end do
    dx_coarse=x_coarse(2:)-x_coarse(1:end_coarse-1)


    !prepare for jenkins mode
    ! call make_positive_upwards(x_learn,hb_learn)      !only allow for postive angle
    allocate(a(4))
    a=polyfit(x_learn(1:n_gl),hb_learn(1:n_gl),4)
    ! a=polyfit(x_learn,hb_learn,4)
    call polys_fine(x_coarse,4,hb_coarse,a)                                    !polynomial fit for smooth shelf
    write(*,*)'hbcora', hb_coarse(1) 
    hb_coarse=hb_coarse-hb_coarse(1)+hbgl
    write(*,*)'hbcora after', hb_coarse(1) 
    write(*,*)'has tobe hbgl',hbgl
    write(*,*)'fine grid'
    where(hb_coarse>0)hb_coarse=0d0                                     !keep hb always under sealevel
    dhb_coarse=hb_coarse(2:)-hb_coarse(1:end_coarse-1)


    ds_coarse=(dx_coarse**2+dhb_coarse**2)**0.5
    s_coarse=cumsum(ds_coarse)

    !************determine all on Fine grind

    !*****determins regualar s with ds and last cell smaller step
    steps_ds=int(s_coarse(end_coarse)/ds)+1               
    dx_reg=((x_coarse(end_coarse)-x_coarse(1))/(steps_ds))   !step size for regular dx
    allocate(s(steps_ds+1))  
    allocate(hb(steps_ds+1))   
    allocate(x(steps_ds+1))   
    allocate(s_reg_dx(steps_ds+1))   !regular in dx (reg_x)
    allocate(x_reg_dx(steps_ds+1))     
    allocate(hb_reg_dx(steps_ds+1))
    allocate(dhb_reg_dx(steps_ds))
    allocate(ds_reg_dx(steps_ds))
    !create regualar x with belonging hb to detemin ds 
    x_reg_dx(1:steps_ds+1)=[(x_coarse(1)+i*dx_reg,i=0,steps_ds)]   ! stepzie=ds for regular ds
    ! write(*,*)' x_reg_dx(1)',x_reg_dx(1)
    x_reg_dx(steps_ds+1)=x_coarse(end_coarse)          ! last steo with stepsize < ds
    call polys_fine(x_reg_dx,4,hb_reg_dx,a)                                    !polynomial fit for smooth shelf
    hb_reg_dx=hb_reg_dx-hb_reg_dx(1)+hbgl
    write(*,*)'hb_reg(1)',hb_reg_dx(1)
    where(hb_reg_dx>0)hb_reg_dx=0d0                                     !kepp hb always under sealevel
    ! call make_positive_upwards(x_reg_dx,hb_reg_dx)
    dhb_reg_dx=hb_reg_dx(2:)-hb_reg_dx(1:steps_ds-1)
    ds_reg_dx=(dx_reg**2+dhb_reg_dx**2)**0.5
    s_reg_dx=cumsum(ds_reg_dx)
    !now interpolate to reg in s grid for jenkinsmode
    s(1:steps_ds)=[(0d0+i*ds,i=0,steps_ds-1)]   ! stepzie=ds for regular ds
    s(steps_ds+1)=s_coarse(end_coarse)          ! last steo with stepsize < ds
    hb=interpolate(s_reg_dx,hb_reg_dx,s)
    x=interpolate(s_reg_dx,x_reg_dx,s)
    ! write(*,*)' x__ds(1)',x(1)
    sina=( hb(2:) - hb(1:steps_ds) ) / ( x(2:) - x(1:steps_ds) )
    real_len=size(s)
    s_in(1:real_len)=s
    x_in(1:real_len)=x
    sina_in(1:real_len-1)=sina
    hb_in(1:real_len)=hb
    dx_coarse_back(2:)=dx_coarse
    dx_coarse_back(1)=dx_coarse(1)
   
  end subroutine
  subroutine floating_tongue_profile_half(ds,glacier_x,glacier_dx,glacier_hb,hbgl,n_gl,s_in,x_in,hb_in,sina_in,&
                                                                                    real_len,x_coarse,dx_coarse_back)
    double precision,intent(in)    :: glacier_x(:)
    double precision,intent(in)    :: glacier_dx(:)
    double precision,intent(in)    :: glacier_hb(:)
    double precision,intent(inout)    :: s_in(:),x_in(:),hb_in(:),sina_in(:)
    integer,intent(inout)    :: real_len
    double precision, intent(in)   :: hbgl
    double precision,intent(in)    :: ds
    integer,intent(in)    :: n_gl    !=size(glacier_x))
    double precision               :: s_coarse(2*n_gl-1),hb_learn(n_gl+1),x_learn(n_gl+1)
    double precision, intent(inout)               :: x_coarse(2*n_gl-1),dx_coarse_back(2*n_gl-1)
    double precision               :: hb_coarse(2*n_gl-1),ds_coarse(2*n_gl-2),dhb_coarse(2*n_gl-2),dx_coarse(2*n_gl-2)
    double precision,dimension(:) ,allocatable                :: s,x,hb,sina,a,x_reg_dx,hb_reg_dx,s_reg_dx,ds_reg_dx,dhb_reg_dx!profile for jenkins calculation on fine grid ds
  
    double precision                :: dx_reg
    integer                                     :: end_coarse,steps_ds,i

    ! true floating part on coarse glacier grid:
  
    end_coarse=2*n_gl-1
    !for polynomial fit
    hb_learn(1:n_gl)=glacier_hb
    x_learn(1:n_gl)=glacier_x
    x_learn(1)=glacier_x(1)+0.5*glacier_dx(1)
    x_learn(n_gl+1)=glacier_x(n_gl)+0.5*glacier_dx(n_gl)
    x_coarse(1)=x_learn(1)
    hb_learn(1)=hbgl
    hb_learn(n_gl+1)=0d0
    ! dx_coarse=0.5d0*glacier_dx
    i=1
    do while (i<=end_coarse/2)
      x_coarse(2*i)=glacier_x(i+1)
      x_coarse(2*i+1)=glacier_x(i+1)+0.5*glacier_dx(i+1)
      i=i+1
    end do
    dx_coarse=x_coarse(2:)-x_coarse(1:end_coarse-1)


    !prepare for jenkins mode
    call make_positive_upwards(x_learn,hb_learn)      !only allow for postive angle
    allocate(a(4))
    a=polyfit(x_learn,hb_learn,4)
    call polys_fine(x_coarse,4,hb_coarse,a)                                    !polynomial fit for smooth shelf
    hb_coarse=hb_coarse-hb_coarse(1)+hbgl
    where(hb_coarse>0)hb_coarse=0d0                                     !kepp hb always under sealevel
    dhb_coarse=hb_coarse(2:)-hb_coarse(1:end_coarse-1)


    ds_coarse=(dx_coarse**2+dhb_coarse**2)**0.5
    s_coarse=cumsum(ds_coarse)

    !************determine all on Fine grind

    !*****determins regualar s with ds and last cell smaller step
    steps_ds=int(s_coarse(end_coarse)/ds)+1               
    dx_reg=((x_coarse(end_coarse)-x_coarse(1))/(steps_ds))   !step size for regular dx
    allocate(s(steps_ds+1))  
    allocate(hb(steps_ds+1))   
    allocate(x(steps_ds+1))   
    allocate(s_reg_dx(steps_ds+1))   !regular in dx (reg_x)
    allocate(x_reg_dx(steps_ds+1))     
    allocate(hb_reg_dx(steps_ds+1))
    allocate(dhb_reg_dx(steps_ds))
    allocate(ds_reg_dx(steps_ds))
    !create regualar x with belonging hb to detemin ds 
    x_reg_dx(1:steps_ds+1)=[(x_coarse(1)+i*dx_reg,i=0,steps_ds)]   ! stepzie=ds for regular ds
    x_reg_dx(steps_ds+1)=x_coarse(end_coarse)          ! last steo with stepsize < ds
    call polys_fine(x_reg_dx,4,hb_reg_dx,a)                                    !polynomial fit for smooth shelf
    hb_reg_dx=hb_reg_dx-hb_reg_dx(1)+hbgl
    where(hb_reg_dx>0)hb_reg_dx=0d0                                     !kepp hb always under sealevel
    ! call make_positive_upwards(x_reg_dx,hb_reg_dx)
    dhb_reg_dx=hb_reg_dx(2:)-hb_reg_dx(1:steps_ds-1)
    ds_reg_dx=(dx_reg**2+dhb_reg_dx**2)**0.5
    s_reg_dx=cumsum(ds_reg_dx)
    !now interpolate to reg in s grid for jenkinsmode
    s(1:steps_ds)=[(0d0+i*ds,i=0,steps_ds-1)]   ! stepzie=ds for regular ds
    s(steps_ds+1)=s_coarse(end_coarse)          ! last steo with stepsize < ds
    hb=interpolate(s_reg_dx,hb_reg_dx,s)
    x=interpolate(s_reg_dx,x_reg_dx,s)
    sina=( hb(2:) - hb(1:steps_ds) ) / ( x(2:) - x(1:steps_ds) )
    real_len=size(s)
    s_in(1:real_len)=s
    x_in(1:real_len)=x
    sina_in(1:real_len-1)=sina
    hb_in(1:real_len)=hb
    dx_coarse_back(2:)=dx_coarse
    dx_coarse_back(1)=dx_coarse(1)
   
  end subroutine

  subroutine floating_tongue_profile(ds,glacier_x,glacier_dx,glacier_hb,hbgl,n_gl,s_in,x_in,hb_in,sina_in,real_len)
    double precision,intent(in)    :: glacier_x(:)
    double precision,intent(in)    :: glacier_dx(:)
    double precision,intent(in)    :: glacier_hb(:)
    double precision,intent(inout)    :: s_in(:),x_in(:),hb_in(:),sina_in(:)
    integer,intent(inout)    :: real_len
    double precision, intent(in)   :: hbgl
    double precision,intent(in)    :: ds
    integer,intent(in)    :: n_gl    !=size(glacier_x))
    double precision               :: hb_coarse(n_gl+1)
    double precision               :: x_coarse(n_gl+1),s_coarse(n_gl+1)
    double precision               :: dx_coarse(n_gl),ds_coarse(n_gl),dhb_coarse(n_gl)
    double precision,dimension(:) ,allocatable                :: s,x,hb,sina !profile for jenkins calculation on fine grid ds
    integer                                     :: end_coarse,steps_ds,i

    ! true floating part on coarse glacier grid:
    end_coarse=n_gl+1
    hb_coarse(1:n_gl)=glacier_hb
    x_coarse(1:n_gl)=glacier_x
    hb_coarse(1)=hbgl
    hb_coarse(end_coarse)=0d0
    dx_coarse=glacier_dx
    dx_coarse(1)=0.5d0*glacier_dx(1)

    dx_coarse(n_gl)=0.5d0*glacier_dx(size(glacier_dx))
    x_coarse(1)=glacier_x(1)+dx_coarse(1)
    x_coarse(end_coarse)=glacier_x(n_gl)+dx_coarse(n_gl)

    !prepare for jenkins mode
    call make_positive_upwards(x_coarse,hb_coarse)      !only allow for postive angle
    call polys(x_coarse,hb_coarse,4)                                    !polynomial fit for smooth shelf
    where(hb_coarse>0)hb_coarse=0d0                                     !kepp hb always under sealevel
    dhb_coarse=hb_coarse(2:)-hb_coarse(1:n_gl)
    
    !determine s on coarse grid
    
    ds_coarse=(dx_coarse**2+dhb_coarse**2)**0.5
    s_coarse=cumsum(ds_coarse)

    !************determine all on Fine grind

    !*****determins regualar s with ds and last cell smaller step
    steps_ds=int(s_coarse(end_coarse)/ds)+1
    allocate(s(steps_ds+1))
    allocate(x(steps_ds+1))
    allocate(hb(steps_ds+1))
    allocate(sina(steps_ds))
    s(1:steps_ds)=[(0d0+i*ds,i=0,steps_ds-1)]   ! stepzie=ds
    s(steps_ds+1)=s_coarse(end_coarse)          ! last steo with stepsize < ds
    
    hb=interpolate(s_coarse,hb_coarse,s)
    x=interpolate(s_coarse,x_coarse,s)
    sina=( hb(2:) - hb(1:steps_ds) ) / ( x(2:) - x(1:steps_ds) )
    real_len=size(s)
    s_in(1:real_len)=s
    x_in(1:real_len)=x
    sina_in(1:real_len-1)=sina
    hb_in(1:real_len)=hb
  end subroutine

  subroutine vertical_tongue_profile(ds,hbgl,s_in,hb_in,sina_in,real_len)
    double precision, intent(in)   :: hbgl
    double precision,intent(in)    :: ds
    double precision,intent(inout)    :: s_in(:),hb_in(:),sina_in(:)
    integer,intent(inout)    :: real_len
    double precision,dimension(:) ,allocatable                :: s,hb,sina !profile for jenkins calculation on fine grid ds
    integer                                     :: steps_ds,i

    ! true floating part on coarse glacier grid:
    !*****determins regualar s with ds and last cell smaller step
    steps_ds=int(abs(hbgl)/ds)+1
    allocate(s(steps_ds+1))
    allocate(hb(steps_ds+1))
    allocate(sina(steps_ds))
    s(1:steps_ds)=[(0d0+i*ds,i=0,steps_ds-1)]   ! stepzie=ds
    s(steps_ds+1)= abs(hbgl)         ! last steo with stepsize < ds
    hb=hbgl+s
    sina=1d0
    real_len=size(s)
    ! s_in(1:steps_ds+1)=s
    s_in(1:real_len)=s
    sina_in(1:real_len-1)=sina
    hb_in(1:real_len)=hb
  end subroutine

  subroutine temperature_profile(Ta_bottom,Ta_surface,d1,d2,Z0,steps,T,Z)

    real, parameter :: pi = 3.1415927
    double precision, intent(in):: Ta_bottom
    double precision, intent(in):: Ta_surface
    double precision, intent(in):: d1  !lower stratification
    double precision, intent(in):: d2  !upper stratifictaion
    double precision, intent(in):: Z0
    integer, intent(in)::steps
    double precision, intent(out):: T(steps)
    double precision, intent(out):: Z(steps)
    integer :: i
    integer :: c1,c2
    double precision :: step_size,Ta0
    if(d1>d2)then
      print*,'d1 has to be lower then d2'
      stop
    end if
    Z=Z0
    step_size=-1*Z0/(steps-1)
    Z(1:steps)=[(Z0+((i-1)*step_size),i=1,steps)]
    Z(steps)=0d0
    c1=1
    c2=steps
    T=Ta_bottom
    Ta0=Ta_bottom
    if(Z0<d2)then
      i=steps
      do i=steps,1,-1
        if(Z(i)<=d2)then
          c2=i
          exit
        end if
      end do
      if(Z0<=d1)then
        i=1
        do i=1,steps
          if(Z(i)>=d1)then
            c1=i
            Ta0=Ta_bottom
            exit
          end if
        end do
      end if
      if(Z0>d1)then
        c1=1
        Ta0=(Ta_surface-Ta_bottom)/(d2-d1)*(Z0-d1)+Ta_bottom
        T(1)=Ta0
      end if
      !!!!! linear
      ! amplitude=(Ta_surface-Ta0)/(c2-c1)
      ! c_xa=Ta_bottom-amplitude*c1
      i=1
      do i=1,steps
        if(i>c1 .and. i<c2)then
          T(i)=(Ta_surface-Ta0)/(c2-c1)*(i-c1)+Ta0
        end if  
        if(i>=c2)then
          T(i)=Ta_surface
        end if
      end do
    else
      T=Ta_surface
    end if
  end subroutine temperature_profile



  subroutine melting_process(S,Z,Ta,Sa,E0,sina,fjordmelt,y0,values,which_jenkins,yout)
    !all on fine grid

    double precision,intent(in)::  S(:)
    double precision, intent(in)::Sa(:)
    double precision, intent(in)::Ta(:)
    double precision, intent(in)::Z(:)
    double precision, intent(in)::y0(:)
    double precision, intent(in)::E0
    double precision, intent(in)::sina(:)
    double precision, intent(out):: fjordmelt(size(S))
    double precision, intent(out)::yout(size(y0)) 
    character(20), intent(in)::which_jenkins
    integer, intent(in)::values
    double precision::y(values,4)
    double precision:: Ufac(values)
    double precision:: U(values)
    double precision, allocatable ::u_cone(:),d_cone(:)
    double precision :: melt,Sb
    double precision:: S_plume,T_plume,k1(4),k2(4),k3(4),k4(4),h
    double precision,parameter:: beta_S      =  7.86d-04 ! haline density coefficient 
    double precision,parameter::beta_T      =  3.87d-05 ! [1/degC]thermal density coefficient
    integer::i
    y(1,:)=y0
    i=1
    fjordmelt=0
    Ufac=0d0
    U=0d0
    do while(i<values)
      T_plume=y(i,3)/y(i,1)
      S_plume=y(i,4)/y(i,1)
      melt=0
      yout=y(i,:)
      call cal_Sb_melt(T_plume,S_plume,Z(i),y(i,:),melt,Sb) !first calculate Sb and meltrate
      ! if(melt<=0)exit  !no refreezing
      if(Sb<0)exit
      if(y(i,4)/y(i,1)<0)exit
      if(y(i,2)<0)exit
      if(y(i,1)<0)exit
      h=S(i+1)-S(i)
      U(i)=y(i,2)/y(i,1)
      if(U(i)<=0)then
        write(*,*)'negativ velocoty i',i
      end if
      if(i>3)then    !#HACK!
        Ufac(i-1)=U(i)/U(i-1)
        Ufac(i-2)=U(i-1)/U(i-2)
        if(Ufac(i-1)/Ufac(i-2)>=2)then
          melt=0d0
          exit 
        end if
     end if

      fjordmelt(i)=melt
      if(which_jenkins=='cone')then
        !cone cant resolve negative veloceties du to sqrt
        k1= ode(Sa(i),Ta(i),Z(i),Sb,melt,E0,sina(i),y(i,:),which_jenkins)*h

        if((y(i,2)+k1(2)/2)<0)exit

        k2= ode(Sa(i),Ta(i),Z(i),Sb,melt,E0,sina(i),y(i,:)+k1/2,which_jenkins)*h

        if((y(i,2)+k2(2)/2)<0)exit

        k3= ode(Sa(i),Ta(i),Z(i),Sb,melt,E0,sina(i),y(i,:)+k2/2,which_jenkins)*h

        if((y(i,2)+k3(2))<0)exit

        k4= ode(Sa(i),Ta(i),Z(i),Sb,melt,E0,sina(i),y(i,:)+k3,which_jenkins)*h
        y(i+1,:)=y(i,:)+1d0/6d0*(k1+2*k2+2*k3+k4)
        i=i+1
      else
        k1= ode(Sa(i),Ta(i),Z(i),Sb,melt,E0,sina(i),y(i,:),which_jenkins)*h
        k2= ode(Sa(i),Ta(i),Z(i),Sb,melt,E0,sina(i),y(i,:)+k1/2,which_jenkins)*h
        k3= ode(Sa(i),Ta(i),Z(i),Sb,melt,E0,sina(i),y(i,:)+k2/2,which_jenkins)*h
        k4= ode(Sa(i),Ta(i),Z(i),Sb,melt,E0,sina(i),y(i,:)+k3,which_jenkins)*h
        y(i+1,:)=y(i,:)+1d0/6d0*(k1+2*k2+2*k3+k4)
        i=i+1
      end if
    
      if(i==values)then
        call cal_Sb_melt(T_plume,S_plume,Z(i),y(i,:),melt,Sb) !first calculate Sb and meltrate
        if(Sb<0)exit
        if(y(i,4)/y(i,1)<0)exit
        if(y(i,2)<0)exit
        if(y(i,1)<0)exit
        fjordmelt(i)=melt
      end if
    end do
    !fjordmelt is along s- coordinate!
    ! if(which_jenkins=='cone')then
    !   allocate(u_cone(i-1))
    !   allocate(d_cone(i-1))
    !   u_cone=y(1:i-1,2)/y(1:i-1,1)
    !   d_cone=(y(1:i-1,1)/u_cone)**0.5
    !   fjordmelt(1:i-1)=(fjordmelt(1:i-1)*d_cone)/d_cone(i-1)
    ! end if


  end subroutine melting_process
  subroutine init_y0(Ta,Sa,sina,E0,DU0,y0,method)
    double precision, intent(in):: Ta
    double precision, intent(in):: Sa
    double precision, intent(in):: E0
    double precision, intent(in):: sina
    double precision, intent(in):: DU0       !in m/s
    character(20),intent(in)    :: method
    double precision, intent(out)::y0(4)
    double precision            :: Uj
    double precision, parameter :: Pi=3.141592653589793238462643d0
    if(method== 'jenkins_1d')then
      Uj= Ujenkins(Ta,Sa,0d0,1d-6,sina,E0,DU0)
      ! setting plume temperture to 0
      !settinig plume salininty to almost 0 (1d-6) due to numerical resaons
      y0(1:4)=(/DU0,DU0*Uj,0d0,DU0*1d-6/)
    else
      Uj= Ujenkins_cone(Ta,Sa,0d0,1d-6,sina,E0,DU0)
      y0(1:4)=(/DU0*2d0/Pi,DU0*2d0/Pi*Uj,0d0,DU0*2d0/Pi*1d-6/)
    end if
  end subroutine

  Function Ujenkins_cone(Ta,Sa,T,S,sina,E0,DU0)
    double precision, intent(in)::Ta
    double precision, intent(in)::Sa
    double precision, intent(in)::E0
    double precision, intent(in)::sina
    double precision, intent(in)::T
    double precision, intent(in)::S
    double precision, intent(in)::DU0
    double precision,parameter:: beta_S      =  7.86d-04 ! haline density coefficient 
    double precision,parameter::beta_T      =  3.87d-05 ! [1/degC]thermal density coefficient
    double precision, parameter :: Cd=2.5d-3
    double precision, parameter :: Pi=3.141592653589793238462643d0
    double precision::rho 
    double precision::Ujenkins_cone
    rho=beta_S*(Sa-S)-beta_T*(Ta-T)
    Ujenkins_cone=(sina*9.81*rho*DU0/((DU0*2*Pi)**(1/2.)*E0*sina+2*Cd/Pi))**(2d0/5d0)

  end Function
  Function Ujenkins(Ta,Sa,T,S,sina,E0,DU0)
    double precision, intent(in)::Ta
    double precision, intent(in)::Sa
    double precision, intent(in)::E0
    double precision, intent(in)::sina
    double precision, intent(in)::T
    double precision, intent(in)::S
    double precision, intent(in)::DU0
    double precision,parameter:: beta_S      =  7.86d-04 ! haline density coefficient 
    double precision,parameter::beta_T      =  3.87d-05 ! [1/degC]thermal density coefficient
    double precision, parameter :: Cd=2.5d-3
    double precision::rho 
    double precision::Ujenkins
    rho=beta_S*(Sa-S)-beta_T*(Ta-T)
    Ujenkins=(sina*9.81*rho*DU0/(E0*sina+Cd))**(1d0/3d0)

  end Function
  Function interpolate(x_in,y_in,x)
    double precision, intent(in)::x_in(:)
    double precision, intent(in)::y_in(:)
    double precision, intent(in)::x(:)
    double precision,dimension(size(x))::interpolate


    interpolate= interp1d(x_in, y_in, x, &
                bounds_error = .false., & ! do not raise out-of-bounds error
                        fill_value = y_in(size(x_in)))

  end Function
  Function ode(Sa,Ta,Z,Sb,melt,E0,sina,y,name_dim)
    implicit none

    double precision, intent(in)::Sa
    double precision, intent(in)::Ta
    double precision, intent(in)::Z
    double precision, intent(in)::E0
    double precision, intent(in)::sina
    double precision, intent(in)::Sb
    double precision, intent(in)::melt
    double precision, intent(in)::y(:)
    character(20),intent(in)::name_dim
    double precision,parameter:: beta_S      =  7.86d-04 ! haline density coefficient 
    double precision,parameter::beta_T      =  3.87d-05 ! [1/degC]thermal density coefficient
    double precision, parameter :: lambda1=-5.73d-02
    double precision, parameter :: lambda2=8.32d-02
    double precision, parameter :: lambda3=7.61d-04
         
    double precision, parameter :: Cd=2.5d-3
    double precision, parameter :: c_a=3.974d03 !specifice heatspecifice heat capacity for seawater
    double precision, parameter :: ci=2.009d03 

    double precision, parameter :: Cd12GammaTS=5.9d-4
    double precision,parameter:: Cd12GammaT  =  1.1d-3
    double precision,parameter:: Cd12GammaS  =  3.1d-5
    double precision, parameter :: L=3.35d+05
    double precision, parameter :: Pi=3.141592653589793238462643d0
    double precision, parameter :: g=9.81d0
    double precision, parameter :: Ti=-10d0
    double precision, parameter :: Si=0d0
    double precision, dimension(4) :: ode
    double precision::y1,y2,y3,y4
    if(name_dim=='jenkins_1d') then

      y1=((E0*sina)*y(2)/y(1))+melt
      y2=((g*sina)*y(1)**2/y(2)*(beta_S*(Sa-y(4)/y(1))-beta_T*(Ta-y(3)/y(1)))-Cd*y(2)*y(2)/(y(1)*y(1)))
      y3=(E0*sina*y(2)/y(1)*Ta+melt*(lambda1*Sb+lambda2+lambda3*Z)-Cd12GammaT*y(2)/y(1)*(y(3)/y(1)-(lambda1*Sb+lambda2+lambda3*Z)))
      y4=(E0*sina*y(2)/y(1)*Sa+melt*Sb-Cd12GammaS*y(2)/y(1)*(y(4)/y(1)-Sb) )

    else if(name_dim=='cone') then
      ! if(y(1)<0 .or.y(2)<0)then
      !   y1=-888
      !   y2=-888
      !   y3=-888
      !   y4=-888
      ! else
      y1= (2*E0*sina)*abs(y(2)**(1d0/2d0))+melt*4/Pi*y(1)/abs(y(2)**(1d0/2d0))
      y2=((g*sina)*(y(1)**2d0)/y(2)*(beta_S*(Sa-y(4)/y(1))-beta_T*(Ta-y(3)/y(1)))-4d0*Cd/Pi*abs(y(2)**(3d0/2d0))/y(1))
      y3=(2*E0*sina*abs(y(2)**(1d0/2d0))*Ta+4d0/Pi*melt*(lambda1*Sb+lambda2+lambda3*Z)*y(1)/abs(y(2)**(1d0/2d0))-&
        4/Pi*Cd12GammaT*abs(y(2)**(1d0/2d0))*(y(3)/y(1)-(lambda1*Sb+lambda2+lambda3*Z)))
      y4=(2*E0*sina*abs(y(2)**(1d0/2d0))*Sa+4d0/Pi*melt*Sb*y(1)/abs(y(2)**(1d0/2d0))-&
        4/Pi*Cd12GammaS*abs(y(2)**(1d0/2d0))*(y(4)/y(1)-Sb) )
      ! end if
    else
      stop("no jenkins mode is chosen!")
    end if
    

    ode=(/y1,y2,y3,y4/)

  end Function ode

  subroutine cal_Sb_melt(T,S,Z,y0,melt,Sb)
    double precision, intent(in) ::T,S,Z
    double precision,intent(in)::y0(:)
    double precision,intent(out)::Sb
    double precision,intent(inout)::melt
    double precision, parameter :: lambda1=-5.73d-02
    double precision, parameter :: lambda2=8.32d-02
    double precision, parameter :: lambda3=7.61d-04
         
    double precision, parameter :: Cd=2.5d-3
    double precision, parameter :: c_a=3.974d03 !specifice heatspecifice heat capacity for seawater
    double precision, parameter :: ci=2.009d03 

    double precision, parameter :: Cd12GammaTS=5.9d-4
    double precision,parameter:: Cd12GammaT  =  1.1d-3
    double precision,parameter:: Cd12GammaS  =  3.1d-5
    double precision, parameter :: L=3.35d+05
    double precision, parameter :: Ti=-10d0
    double precision, parameter :: Si=0d0
    double precision :: A,B,C,p,q,Tb
    A=(c_a*Cd12GammaT-ci*Cd12GammaS)*lambda1
    B=Cd12GammaS*ci*(Ti+lambda1*S-lambda2-lambda3*Z)-Cd12GammaT*c_a*(T+lambda1*Si-lambda2-lambda3*Z)-Cd12GammaS*L
    C=(Cd12GammaS*ci*(lambda2+lambda3*Z-Ti)+Cd12GammaS*L)*S-Cd12GammaT*(lambda2+lambda3*Z-T)*Si
    p=B/A
    q=C/A

    if ((p/2)**2-q<0) then
      Sb=-100
    else
      Sb=-p/2+abs(sqrt((p/2)**2-q))
    end if
    Tb=lambda1*Sb+lambda2+lambda3*Z
    melt=Cd12GammaT*c_a*y0(2)/y0(1)*(y0(3)/y0(1)-Tb)/(L+ci*(Tb-Ti))
    ! /call the_ode_system()
  end subroutine cal_Sb_melt
  subroutine submarine_melt_test(E0,ds,glacier_x,glacier_dx,glacier_hb,glacier_zb,method,discharge, fjordtestmelt,length &
      ,s_t,x_t,hb_t,sina_t,sa_t,ta_t,y_t,the_length,x_learn,hb_learn)
    double precision,intent(in)      :: ds
    double precision,dimension(:),intent(in)      :: glacier_x
    double precision,dimension(:),intent(in)      :: glacier_dx
    double precision,dimension(:),intent(in)      :: glacier_zb
    double precision,dimension(:),intent(in)      :: glacier_hb
    ! double precision,intent(out)                  :: fjordmelt(size(glacier_hb))
    double precision                  :: fjordmelt(size(glacier_hb))
    double precision,dimension(length),intent(out)                  ::fjordtestmelt,s_t,x_t,hb_t,sina_t,sa_t,ta_t
    double precision                              :: x_half(2*size(glacier_x)-1),dx_half(2*size(glacier_x)-1)
    double precision,intent(out)                              :: x_learn(size(glacier_x)+1),hb_learn(size(glacier_x)+1)
    double precision                              :: glacier_x_right(size(glacier_x))
    character(20),intent(in) :: method
    integer,intent(in) ::length
    real,optional,intent(in) :: discharge
    double precision                 :: DU0,hbgl
    double precision,intent(in)                 ::E0 
    character(20)            :: glacier_type        ! tidewater or floating tongue
    double precision,dimension(:),allocatable                 :: s,x,hb,sina,melt,Ta,Sa !profile for jenkins calculation on fine grid ds
    double precision                              :: y0(4),yout(4)      !staring vlaues of plumer for vVypelocity flux andT,S
    double precision,intent(out)                              :: y_t(4)      !staring vlaues of plumer for vVypelocity flux andT,S
    integer,intent(out)                         ::the_length
    integer :: maxlen,n_gl,real_len

    n_gl=size(glacier_x)
    ! ds= 2d0                     !fine grid for jenkins mode: resolution  1 m


    DU0=1d0              !minimal discharge value
    ! E0=0.036
    hbgl=0.5*(glacier_zb(1)+glacier_zb(2))  !true depth of shelf
    glacier_x_right=glacier_x+0.5*glacier_dx

    !parameters
    glacier_type='tidewater'    ! default type
    if(size(glacier_hb)>1)glacier_type='floating'
    if(present(discharge))DU0=discharge


    !compute fine jenkins profile including interpolation
    maxlen=int((glacier_x(n_gl)+glacier_dx(n_gl)-glacier_x(1)+abs(glacier_hb(1)))/ds+1)
    allocate(s(maxlen))
    allocate(sina(maxlen))
    allocate(hb(maxlen))
    allocate(x(maxlen))
    s=-9999d0  !missing values
    sina=-9999d0
    x=-9999d0
    hb=-9999d0
    if(glacier_type=='floating')then
      !calculate floating tongue profile

      ! call floating_tongue_profile(ds,glacier_x,glacier_dx,glacier_hb,hbgl,size(glacier_x),s,x,hb,sina,real_len)
      ! call floating_tongue_profile_half(ds,glacier_x,glacier_dx,glacier_hb,hbgl,size(glacier_x),s,x,hb,sina,real_len,x_half,dx_half)
      call floating_tongue_profile_half_test(ds,glacier_x,glacier_dx,glacier_hb,hbgl,size(glacier_x),&
        s,x,hb,sina,real_len,x_half,dx_half, x_learn,hb_learn)
    else
      call vertical_tongue_profile(ds,hbgl,s,hb,sina,real_len)
    endif
    write(*,*)'real_len',real_len
    write(*,*)'x[1]',x(1)
    write(*,*)'glacier_x',glacier_x(1)
    hb_t(1:real_len)=hb
    x_t(1:real_len)=x
    s_t(1:real_len)=s
    sina_t(1:real_len)=sina
    !calculate T/S profile
    allocate(Ta(real_len))
    allocate(Sa(real_len))
    call compute_oceanographich_profile(hb(1:real_len),Ta,Sa)  !ambient Tempertatur and salinity
    sa_t(1:real_len)=Sa
    ta_t(1:real_len)=Ta

    !melting:
    allocate(melt(real_len))
    call init_y0(Ta(1),Sa(1),sina(1),E0,DU0,y0,method)
    y_t=y0
    melt=0d0
    call melting_process(s(1:real_len),hb(1:real_len),Ta(1:real_len),Sa(1:real_len),E0,sina,melt,y0,real_len,method,yout)
   
    fjordtestmelt=0d0
    fjordtestmelt(1:real_len)=melt

    the_length=real_len

    ! transform back of melting
    fjordmelt=0d0
    if(glacier_type=='floating')then  ! need to transform in horizontal and vertial melting

      call aggregate_submarine_melt_floating(ds,x(1:real_len),hb(1:real_len),glacier_x_right,glacier_dx,melt,fjordmelt)   ! important: use x of cell edge
      ! call aggregate_submarine_melt_angle(ds,x(1:real_len),hb(1:real_len),glacier_x_right,glacier_dx,melt,fjordmelt)

    else      !all vertical melting
      fjordmelt=sum(melt)*ds/glacier_dx(1)          !todo think about horizontal and vertical

    endif
    ! call aggregate_submarine_melt(s,x,hb,glacier_x_right,melt,glacier_type,fjordmelt)   ! important: use x of cell edge
    ! sum s til x > x_right


  end subroutine
!   SUBROUTINE polint(xa,ya,n,x,y,dy)
!     INTEGER n,NMAX
!     REAL dy,x,y,xa(n),ya(n)
!     PARAMETER (NMAX=10) 
!     INTEGER i,m,ns
!     REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
!     ns=1
!     dif=abs(x-xa(1))
!     do i=1,n 
!       dift=abs(x-xa(i))
!       if (dift.lt.dif) then
!         ns=i
!         dif=dift
!       endif
!       c(i)=ya(i) 
!       d(i)=ya(i)
!     end do
!     y=ya(ns) 
!     ns=ns-1
!     do  m=1,n-1
!       do i=1,n-m 
!         ho=xa(i)-x
!         hp=xa(i+m)-x
!         w=c(i+1)-d(i)
!         den=ho-hp
!         if(den.eq.0.)stop('failure in polint')
!         den=w/den
!         d(i)=hp*den
!         c(i)=ho*den
!       enddo 
!       if (2*ns.lt.n-m)then
!         dy=c(ns+1)
!       else
!         dy=d(ns)
!         ns=ns-1
!       endif
!       y=y+dy
!     enddo 
!     return
!     END Subroutine
!   subroutine temp_sinh(steps,T,Z)
!   implicit none
!   ! define constants
!
!   ! real, parameter :: k=0.73607166878564445  ! give j for (sinh(k*pi)) 
!   real, parameter :: k= 0.28054984976025232
!   ! real, parameter :: j = 5.
!   real, parameter :: j = 1d0
!   real, parameter :: pi = 3.1415927
!   real, parameter :: Ta_bottom=4.0
!   real, parameter :: Ta_suface=1.0
!   real, parameter :: d1=-200
!   real, parameter :: d2 = -100
!   real,parameter :: Z0=-500d0
!   integer :: i,c1,c2
!   integer,intent(in) ::steps
!   double precision, intent(out):: Z(steps),T(steps)
!   double precision :: step_size,step_pi,amplitude
!   T=Ta_bottom
!   Z=Z0
!   step_size=-1*Z0/(steps-1)
!   print*,'stepsize',step_size
!   Z(1:steps)=[(Z0+((i-1)*step_size),i=1,steps)]
!
!   c1=ceiling(d1-Z0)/step_size
!   print*,'c1',c1
!   print*,'Z(c1)',Z(c1)
!   c2=floor(d2-Z0)/step_size
!   print*,'c2',c2
!   print*,'Z(c2)',Z(c2)
!   step_pi=(pi*k*2)/(c2-c1)
!   print*,'pistet',step_pi
!   amplitude=(Ta_suface-Ta_bottom)/(2*j)
!   i=steps
!   do i=steps,1,-1
!     if(Z(i)<=d2)then
!       write(*,*)'the same c2 no',Z(i),'i',i
!       exit
!     end if
!   end do
!
!   i=1
!   do i=1,steps
!     if(Z(i)>=d1)then
!       write(*,*)'the same c1 no?',Z(i),'i',i
!       exit
!     end if
!   end do
!   write(*,*) 'Hello there'
!   write(*,*) 'zz',Z(0)
!   write(*,*) 'zz',Z(0:0)
!   write(*,*) 'len(z)',size(Z(0:0))
!   write(*,*) 'zz',Z(0:1)
!   write(*,*) 'len(z)',size(Z(0:1))
!   print*,'tan', tan(-1*pi/2)
!   print*,'tan', tan(-1*pi/4)
!   print*,'tan', tan(0*pi/4)
!   print*,'tan', tan(1*pi/4)
!   print*,'tan plus', tan(1*pi/2)
!   i=1
!
!   do i=1,steps
!     if(i>c1 .and. i<c2)then
!       T(i)=sinh(-1*k*pi+(i-c1)*step_pi)*amplitude&
!         +(Ta_bottom-sinh(-1*k*pi)*amplitude)
!
!     end if  
!     if(i>=c2)then
!       T(i)=Ta_suface
!     end if
!     ! print*,'z(i)',Z(i)
!   end do
!
! end subroutine
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
end module
