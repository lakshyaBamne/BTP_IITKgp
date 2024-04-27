    subroutine ini_Riemann(rho,m1,m2,E)
    implicit none 
    integer,parameter :: Nmax=1400
    real*8,dimension (0:Nmax,0:Nmax) :: rho,m1,m2,E,p
    real*8,dimension (1:Nmax) :: x,y
    real*8 ::Ysjn(9),Xwie(9)
    real*8 yN,yS,xW,xE
    real*8,dimension (1:9,1:4) ::qk
    character(20) BoundC,Extype
    integer i,j,Nx,Ny,k
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
    case ('Cfg3')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>1.0d0 .and. y(j)>1.0d0) then
                    rho(i,j)=1.5d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=1.5d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>1.0d0 .and. y(j)<1.0d0) then
                    rho(i,j)=33.d0/62.d0
                    m1(i,j)=0.d0
                    m2(i,j)=rho(i,j)*4.d0/dsqrt(11.d0)
                    E(i,j)=0.3d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<1.0d0 .and. y(j)<1.0d0)then
                    rho(i,j)=77.d0/558.d0
                    m1(i,j)=rho(i,j)*4.d0/dsqrt(11.d0)
                    m2(i,j)=rho(i,j)*4.d0/dsqrt(11.d0)
                    E(i,j)=9.d0/310.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=33.d0/62.d0
                    m1(i,j)=rho(i,j)*4.d0/dsqrt(11.d0)
                    m2(i,j)=0.0d0
                    E(i,j)=0.3d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Cfg5')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx 
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=1.d0
                    m1(i,j)=-0.75d0
                    m2(i,j)=-0.5d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0 .and. y(j)<0.5d0) then
                    rho(i,j)=3.d0
                    m1(i,j)=2.25d0
                    m2(i,j)=-1.5d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=1.d0
                    m1(i,j)=0.75d0
                    m2(i,j)=0.5d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=2.d0
                    m1(i,j)=-1.5d0
                    m2(i,j)=1.0d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Cfg6')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=1.d0
                    m1(i,j)=0.75d0
                    m2(i,j)=-0.5d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=3.d0
                    m1(i,j)=-2.25d0
                    m2(i,j)=-1.5d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=1.d0
                    m1(i,j)=-0.75d0
                    m2(i,j)=0.5d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=2.d0
                    m1(i,j)=1.5d0
                    m2(i,j)=1.0d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
        case ('Cfg7')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=1.0d0
                    m1(i,j)=0.1d0
                    m2(i,j)=0.1d0
                    E(i,j)=1.0d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0 .and. y(j)<0.5d0) then
                    rho(i,j)=0.5197d0
                    m1(i,j)=0.05197d0
                    m2(i,j)=-rho(i,j)*0.6259d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=0.8d0
                    m1(i,j)=0.08d0
                    m2(i,j)=0.08d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=0.5197d0
                    m1(i,j)=-rho(i,j)*0.6259d0
                    m2(i,j)=rho(i,j)*0.1d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Cfg8')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=0.5197d0
                    m1(i,j)=0.05197d0
                    m2(i,j)=0.05197d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=1.d0
                    m1(i,j)=0.1d0
                    m2(i,j)=-0.6259d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=0.8d0
                    m1(i,j)=0.08d0
                    m2(i,j)=0.08d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=1.d0
                    m1(i,j)=-0.6259d0
                    m2(i,j)=0.1d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Cfg11')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=1.d0
                    m1(i,j)=0.1d0
                    m2(i,j)=0.d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0+1.0d-12 .and. y(j)<0.5d0+1.0d-12)then
                    rho(i,j)=0.5313d0
                    m1(i,j)=0.05313d0
                    m2(i,j)=0.3866d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0+1.0d-12 .and. y(j)<0.5d0+1.0d-12)then
                    rho(i,j)=0.8d0
                    m1(i,j)=0.08d0
                    m2(i,j)=0.d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=0.5313d0
                    m1(i,j)=0.4379d0
                    m2(i,j)=0.0d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Cfg12')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=0.5313d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0+1.0d-12 .and. y(j)<0.5d0+1.0d-12) then
                    rho(i,j)=1.0d0
                    m1(i,j)=0.d0
                    m2(i,j)=rho(i,j)*0.7276d0
                    E(i,j)=1.0d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0+1.0d-12 .and. y(j)<0.5d0+1.0d-12)then
                    rho(i,j)=0.8d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=1.0d0
                    m1(i,j)=rho(i,j)*0.7276d0
                    m2(i,j)=0.0d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Cfg13')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=1.d0
                    m1(i,j)=0.d0
                    m2(i,j)=-0.3d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0+1.0d-12 .and. y(j)<0.5d0+1.0d-12)then
                    rho(i,j)=0.5313d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.2272d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0+1.0d-12 .and. y(j)<0.5d0+1.0d-12)then
                    rho(i,j)=1.0625d0
                    m1(i,j)=0.0d0
                    m2(i,j)=0.8654d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=2.d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.6d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Cfg15')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=1.0d0
                    m1(i,j)=0.1d0
                    m2(i,j)=-0.3d0
                    E(i,j)=1.0d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0+1.0d-12 .and. y(j)<0.5d0+1.0d-12) then
                    rho(i,j)=0.5313d0
                    m1(i,j)=0.05313d0
                    m2(i,j)=rho(i,j)*0.4276d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0+1.0d-12 .and. y(j)<0.5d0+1.0d-12)then
                    rho(i,j)=0.8d0
                    m1(i,j)=0.08d0
                    m2(i,j)=-0.24d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=0.5197d0
                    m1(i,j)=-rho(i,j)*0.6259d0
                    m2(i,j)=-rho(i,j)*0.3d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Cfg17')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=1.d0
                    m1(i,j)=0.d0
                    m2(i,j)=-0.4d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=0.5197d0
                    m1(i,j)=0.d0
                    m2(i,j)=-0.5851d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=1.0625d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.2279d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=2.d0
                    m1(i,j)=0.d0
                    m2(i,j)=-0.6d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    case ('Cfg19')
        BoundC='free'
        do j=1,Ny
            do i=1,Nx  
                if (x(i)>0.5d0 .and. y(j)>0.5d0) then
                    rho(i,j)=1.d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.3d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)>0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=0.5197d0
                    m1(i,j)=0.d0
                    m2(i,j)=-0.2213d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                elseif (x(i)<0.5d0 .and. y(j)<0.5d0)then
                    rho(i,j)=1.0625d0
                    m1(i,j)=0.d0
                    m2(i,j)=0.2279d0
                    E(i,j)=0.4d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                else
                    rho(i,j)=2.d0
                    m1(i,j)=0.d0
                    m2(i,j)=-0.6d0
                    E(i,j)=1.d0/(gamma-1.d0)+0.5d0*(m1(i,j)**2+m2(i,j)**2)/rho(i,j)
                endif
            enddo
        enddo
    end select

    call extendcell2D(rho,m1,m2,E)

    return
    end subroutine

