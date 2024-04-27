    subroutine extendcell2D(rho,m1,m2,E)
    implicit none
    integer, parameter :: Nmax=1400
    real*8, dimension(0:Nmax,0:Nmax) ::rho,m1,m2,E
    integer i,j,Nx,Ny
    character(20) BoundC
    common/parameter2/Nx,Ny
    common/parameter4/BoundC

    select case (BoundC)
    case('periodic')
       ! !==periodic boundary====
        do i=1,Nx
            rho(i,0)=rho(i,Ny)
            rho(i,Ny+1)=rho(i,1)
            m1(i,0)=m1(i,Ny)
            m1(i,Ny+1)=m1(i,1)
            m2(i,0)=m2(i,Ny)
            m2(i,Ny+1)=m2(i,1)
            E(i,0)=E(i,Ny)
            E(i,Ny+1)=E(i,1)
        enddo
        do j=1,Ny
            rho(0,j)=rho(Nx,j)
            rho(Nx+1,j)=rho(1,j)
            m1(0,j)=m1(Nx,j)
            m1(Nx+1,j)=m1(1,j)
            m2(0,j)=m2(Nx,j)
            m2(Nx+1,j)=m2(1,j)
            E(0,j)=E(Nx,j)
            E(Nx+1,j)=E(1,j)
        enddo
    case('reflectfree')
        !! ==reflecting boundary(solid wall)+free====
        do i=1,Nx
            rho(i,0)=rho(i,1)
            rho(i,Ny+1)=rho(i,Ny)

            m1(i,0)=m1(i,1)
            m1(i,Ny+1)=m1(i,Ny)
            
            m2(i,0)=-m2(i,1)
            m2(i,Ny+1)=m2(i,Ny)
            
            E(i,0)=E(i,1)
            E(i,Ny+1)=E(i,Ny)
        enddo
        do j=1,Ny
            rho(0,j)=rho(1,j)
            rho(Nx+1,j)=rho(Nx,j)
            
            m1(0,j)=-m1(1,j)
            m1(Nx+1,j)=m1(Nx,j)
            
            m2(0,j)=m2(1,j)
            m2(Nx+1,j)=m2(Nx,j)
            
            E(0,j)=E(1,j)
            E(Nx+1,j)=E(Nx,j)
        enddo
        
    case('free')
        !!==free====
        do i=1,Nx
            rho(i,0)=rho(i,1)
            rho(i,Ny+1)=rho(i,Ny)

            m1(i,0)=m1(i,1)
            m1(i,Ny+1)=m1(i,Ny)
            
            m2(i,0)=m2(i,1)
            m2(i,Ny+1)=m2(i,Ny)
            
            E(i,0)=E(i,1)
            E(i,Ny+1)=E(i,Ny)
        enddo
        do j=1,Ny
            rho(0,j)=rho(1,j)
            rho(Nx+1,j)=rho(Nx,j)
            
            m1(0,j)=m1(1,j)
            m1(Nx+1,j)=m1(Nx,j)
            
            m2(0,j)=m2(1,j)
            m2(Nx+1,j)=m2(Nx,j)
            
            E(0,j)=E(1,j)
            E(Nx+1,j)=E(Nx,j)
        enddo
    case('solidwall')
       !! ==reflecting boundary(solid wall)====
        do i=1,Nx
            rho(i,0)=rho(i,1)
            rho(i,Ny+1)=rho(i,Ny)
            m1(i,0)=m1(i,1)
            m1(i,Ny+1)=m1(i,Ny)
            m2(i,0)=-m2(i,1)
            m2(i,Ny+1)=-m2(i,Ny)
            E(i,0)=E(i,1)
            E(i,Ny+1)=E(i,Ny)
        enddo
        do j=1,Ny
            rho(0,j)=rho(1,j)
            rho(Nx+1,j)=rho(Nx,j)
            m1(0,j)=-m1(1,j)
            m1(Nx+1,j)=-m1(Nx,j)
            m2(0,j)=m2(1,j)
            m2(Nx+1,j)=m2(Nx,j)
            E(0,j)=E(1,j)
            E(Nx+1,j)=E(Nx,j)
        enddo
    end select
    
    return
    end subroutine