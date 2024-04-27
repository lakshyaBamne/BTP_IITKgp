    subroutine initialset_2DRT(rho,m1,m2,E)
    implicit none 
    integer,parameter :: Nmax=1400
    real*8,dimension (0:Nmax,0:Nmax) :: rho,m1,m2,E,p
    real*8,dimension (1:Nmax) :: x,y
    character(20) BoundC
    real*8 pi,u,v,c,grav,yint
    real*8 xW,xE,yN,yS
    integer i,j,Nx,Ny
    real*8 dx,dy,x0,y0,xlength,ylength,dt,Tfinal,t,theta,gamma
    common/parameter1/dx,dy,x0,y0,xlength,ylength,dt,Tfinal,t,theta,gamma
    common/parameter2/Nx,Ny
    common/parameter3/x,y
    common/parameter4/BoundC
    
    BoundC='solidwall'
    do i=1,Nx
        x(i)=x0+(i-0.5d0)*dx
    enddo
    do j=1,Ny
        y(j)=y0+(j-0.5d0)*dy
    enddo
    
    pi=4.d0*datan(1.d0)
    do j=1,Ny
        do i=1,Nx  
            if (y(j)>0.5d0) then
                rho(i,j)=1.d0
                p(i,j)=y(j)+1.5d0
                c=dsqrt(gamma*p(i,j)/rho(i,j))
                if (x(i)<0.125d0) then
                   v=-0.025d0*c*dcos(8.d0*pi*x(i))
                else
                    v=-0.025d0*c*dcos(8.d0*pi*(0.25d0-x(i)))
                endif
                m1(i,j)=0.d0
                m2(i,j)=rho(i,j)*v
                E(i,j)=p(i,j)/(gamma-1.d0)+0.5d0*rho(i,j)*v**2
            else
                rho(i,j)=2.d0
                p(i,j)=2.d0*y(j)+1.d0
                c=dsqrt(gamma*p(i,j)/rho(i,j))
                if (x(i)<0.125d0) then
                   v=-0.025d0*c*dcos(8.d0*pi*x(i))
                else
                    v=-0.025d0*c*dcos(8.d0*pi*(0.25d0-x(i)))
                endif
                m1(i,j)=0.d0
                m2(i,j)=rho(i,j)*v
                E(i,j)=p(i,j)/(gamma-1.d0)+0.5d0*rho(i,j)*v**2
            endif
        enddo
    enddo
    
    !!!=========Rayleigh-Taylor=========
    ! grav=0.1d0
    ! pi=4.d0*datan(1.d0)
    !do j=1,Ny
    !    do i=1,Nx  
    !        if (x(i)<0.d0) then
    !            yint=0.5d0+0.01d0*dcos(6.d0*pi*x(i))
    !        else
    !            yint=0.5d0+0.01d0*dcos(6.d0*pi*(2.d0/6.d0-x(i)))
    !        endif
    !            if (y(j)>yint) then
    !                rho(i,j)=2.d0
    !                m1(i,j)=0.d0
    !                m2(i,j)=0.d0
    !                E(i,j)=(0.1d0+2.d0*grav*(1.d0-y(j)))/(gamma-1.d0)
    !            else
    !                rho(i,j)=1.d0
    !                m1(i,j)=0.d0
    !                m2(i,j)=0.d0
    !                E(i,j)=(0.1d0+grav*(2.d0-yint-y(j)))/(gamma-1.d0)
    !            endif
    !            if (dabs(y(j)-yint)<=2.d0*dy) then
    !                rho(i,j)=1.5d0+0.25d0*(y(j)-yint)/dy
    !            endif
    !    enddo
    !enddo
    !

    call extendcell2D(rho,m1,m2,E)

    return
    end subroutine

