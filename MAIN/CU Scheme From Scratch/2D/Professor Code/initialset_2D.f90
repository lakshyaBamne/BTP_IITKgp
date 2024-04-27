    subroutine initialset_2D(rho,m1,m2,E)
    implicit none 
    integer,parameter :: Nmax=1400
    real*8,dimension (0:Nmax,0:Nmax) :: rho,m1,m2,E,p
    real*8,dimension (1:Nmax) :: x,y
    character(20) BoundC,Extype
    real*8 alpha1,alpha2,beta1,beta2,pi
    real*8 eps,J1,J2,I1,I2,Y1,Y2,u,v,ei
    real*8 grav,yint
    real*8 xW,xE,yN,yS
    real*8 ::Ysjn(3),Xwie(3)
    real*8,dimension (1:3,1:3,4) ::qk
    integer i,j,k1,k2,Nx,Ny
    real*8 dx,dy,x0,y0,xlength,ylength,dt,Tfinal,t,theta,gamma
    common/parameter1/dx,dy,x0,y0,xlength,ylength,dt,Tfinal,t,theta,gamma
    common/parameter2/Nx,Ny
    common/parameter3/x,y
    common/parameter4/BoundC,Extype

    do i=1,Nx
        x(i)=x0+(i-0.5d0)*dx
    enddo
    do j=1,Ny
        y(j)=y0+(j-0.5d0)*dy
    enddo

    select case(Extype)
    case ('MCW')
        !!=========Moving Contact Wave=========
        BoundC='free'

        do j=1,Ny
            do i=1,Nx   
                if ((x(i)<=0.1d0 .and. x(i)>=-0.1d0 .and. y(j)<=0.02d0) .or. &
                (x(i)<=0.02d0 .and. x(i)>=-0.02d0 .and. y(j)>=0.02d0 .and. y(j)<=0.1d0) .or. &
                (x(i)+0.02d0)**2+(y(j)-0.02d0)**2<=0.08d0**2 .or. &
                (x(i)-0.02d0)**2+(y(j)-0.02d0)**2<=0.08d0**2) then
                    rho(i,j)=1.4d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.28d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=1.0d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.2d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Exp') 
        !!=========Explosion=========
        BoundC='reflectfree'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)**2+y(j)**2<0.16d0-1.d-8) then
                    rho(i,j)=1.d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=1.d0/(gamma-1.d0)
                else
                    rho(i,j)=0.125d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=0.1d0/(gamma-1.d0)
                endif
            enddo
        enddo

    case ('Imp')
        !!=========implosion=========
        BoundC='solidwall'
        do j=1,Ny
            do i=1,Nx  
                if (dabs(x(i))+dabs(y(j))<0.15d0-1.0d-8) then
                    rho(i,j)=0.125d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=0.14d0/(gamma-1.d0)
                else
                    rho(i,j)=1.d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=1.d0/(gamma-1.d0)
                endif
            enddo
        enddo
    case ('KH')
        !!=========Kelvin�CHelmholtz=========
        BoundC='periodic'
        pi=4.d0*datan(1.d0)
        I1=0.00625d0
        do j=1,Ny
            do i=1,Nx  
                if (y(j)<-0.25d0) then
                    rho(i,j)=1.d0
                    u=-0.5d0+0.5d0*dexp((y(j)+0.25d0)/I1)
                    v=0.01d0*dsin(4.d0*pi*x(i))
                    ei=3.75d0
                    m1(i,j)=rho(i,j)*u
                    m2(i,j)=rho(i,j)*v
                    E(i,j)=rho(i,j)*(ei+0.5d0*(u**2+v**2))
                elseif(y(j)<0.d0) then
                    rho(i,j)=2.d0
                    u=0.5d0-0.5d0*dexp((-y(j)-0.25d0)/I1)
                    v=0.01d0*dsin(4.d0*pi*x(i))
                    ei=1.875d0
                    m1(i,j)=rho(i,j)*u
                    m2(i,j)=rho(i,j)*v
                    E(i,j)=rho(i,j)*(ei+0.5d0*(u**2+v**2))
                elseif (y(j)<0.25d0) then
                    rho(i,j)=2.d0
                    u=0.5d0-0.5d0*dexp((y(j)-0.25d0)/I1)
                    v=0.01d0*dsin(4.d0*pi*x(i))
                    ei=1.875d0
                    m1(i,j)=rho(i,j)*u
                    m2(i,j)=rho(i,j)*v
                    E(i,j)=rho(i,j)*(ei+0.5d0*(u**2+v**2))
                else
                    rho(i,j)=1.d0
                    u=-0.5d0+0.5d0*dexp((-y(j)+0.25d0)/I1)
                    v=0.01d0*dsin(4.d0*pi*x(i))
                    ei=3.75d0
                    m1(i,j)=rho(i,j)*u
                    m2(i,j)=rho(i,j)*v
                    E(i,j)=rho(i,j)*(ei+0.5d0*(u**2+v**2))
                endif
            enddo
        enddo
    end select 

    call extendcell2D(rho,m1,m2,E)

    return
    end subroutine

