    program Euler_gas_main_1D
    implicit none
    integer,parameter :: Nmax=10000
    real*8,dimension (0:Nmax,1) ::rho0,mx0,E0,rho,mx,E
    real*8,dimension (1:Nmax,1) :: x
    real*8,dimension (0:Nmax,1):: F1,F2,F3,F4
    integer j,Nx,newdtneed,nstep,istep
    character(40) filename,rhon,un,En
    character(20) BoundC,Extype,shrinkSw
    real*8 gamma,dx,x0,xlength,dt,Tfinal,t,theta,dtstep
    common/parameter1/dx,x0,xlength,dt,Tfinal,t,theta,gamma
    common/parameter2/Nx
    common/parameter3/x
    common/parameter4/BoundC,Extype,shrinkSw
    !=======input the data========================

    Nx=1000
    ! Tfinal=0.035d0 !! test
    ! Tfinal=2.0d0 !! MCW
    Tfinal=0.038d0 !! Blast
    ! Tfinal=1.8d0 !! Shu Osher
    ! Tfinal=0.2d0 !! SOD
    ! Tfinal=0.2d0 !! SPP
    ! Tfinal=1.3d0 !! Lax

    x0=0.d0
    xlength=1.d0
    dx=xlength/Nx
    theta=1.3d0
    gamma=1.4d0
    nstep=1
    dtstep=Tfinal/nstep

    ! Extype='MCW'
    ! Extype='Lax'
    Extype='Blast'
    ! Extype='ShuOS'
    ! Extype='SOD'
    ! Extype='SPP'
    ! Extype='test'
    !Extype='SMS'

    shrinkSw='off' ! comment out for CU scheme


    !!===========================initial data=====================================================
    call initialset_1D(rho0,mx0,E0)
    
    t=0.0d0
    do istep=0,nstep
        Tfinal=istep*dtstep
        do while (t<Tfinal-1.0d-7)
            newdtneed=1
            call numerical_flux_Euler1D(rho0,mx0,E0,newdtneed,F1,F2,F3)
            call SSPRK_Eulergas1D(rho0,mx0,E0,F1,F2,F3,rho,mx,E)
            call extendcell1d(rho,mx,E)
            t=t+dt
            write(*,*) 't=', t
            rho0=rho
            mx0=mx
            E0=E
        enddo

        if(istep == 1)then
        write(rhon,101) istep
        open(11,file=rhon)
        !! .....Writing the solution.............................................
        do j=1,Nx
            write(11,*) x(j,1),rho0(j,1)
        enddo
        close(11)
        endif
    enddo


    ! 101 format("out1DMCWrho-CU-ref",i3.3)
    ! 101 format("out1DLaxrho-CU-ref",i3.3)
    101 format("Blast",i3.3)
    ! 101 format("out1DShuOSrho-CU-ref",i3.3)
    ! 101 format("out1DSODrho-CU-ref",i3.3)
    ! 101 format("out1DSPPrhoR",i3.3)
    ! 101 format("out1Dtest41rho2",i3.3)
    !101 format("out1DSMSrrho",i3.3)
    
100 format (' ', f24.16,' ', f24.16,' ', f24.16,' ', f24.16,' ', f24.16,' ', f24.16)
    end program

    real*8 function dminmod(v1,v2,v3)
    implicit none
    real*8 v1,v2,v3,u
      dminmod=0.5d0*(dsign(1.d0,v1)+dsign(1.d0,v3))*dmin1(dabs(v1),dabs(v2),dabs(v3))
    return
    end function


    real*8 function minmod(x,y)
    implicit none
    real*8 x,y
    minmod=0.5d0*(dsign(1.d0,x)+dsign(1.d0,y))*dmin1(dabs(x),dabs(y))
    return
    end function

