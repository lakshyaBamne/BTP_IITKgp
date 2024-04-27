    subroutine extendcell1D(rho,m1,E)
    implicit none
    integer, parameter :: Nmax=10000
    real*8, dimension(0:Nmax,1) ::rho,m1,E
    integer Nx
    character(20) BoundC
    common/parameter2/Nx
    common/parameter4/BoundC

    select case (BoundC)
    case('periodic')
       ! !==periodic boundary====
            rho(0,1)=rho(Nx,1)
            rho(Nx+1,1)=rho(1,1)
            m1(0,1)=m1(Nx,1)
            m1(Nx+1,1)=m1(1,1)
            E(0,1)=E(Nx,1)
            E(Nx+1,1)=E(1,1)
    case('reflectfree')
       !! ==reflecting boundary(solid wall)+free====
            rho(0,1)=rho(1,1)
            rho(Nx+1,1)=rho(Nx,1)
            m1(0,1)=-m1(1,1)
            m1(Nx+1,1)=m1(Nx,1)
            E(0,1)=E(1,1)
            E(Nx+1,1)=E(Nx,1)
        case('free')
        !!==free====
            rho(0,1)=rho(1,1)
            rho(Nx+1,1)=rho(Nx,1)
            m1(0,1)=m1(1,1)
            m1(Nx+1,1)=m1(Nx,1)
            E(0,1)=E(1,1)
            E(Nx+1,1)=E(Nx,1)
        case('solidwall')
       !! ==reflecting boundary(solid wall)====
            rho(0,1)=rho(1,1)
            rho(Nx+1,1)=rho(Nx,1)
            m1(0,1)=-m1(1,1)
            m1(Nx+1,1)=-m1(Nx,1)
            E(0,1)=E(1,1)
            E(Nx+1,1)=E(Nx,1)
    end select
    
    return
    end subroutine