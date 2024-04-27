program Euler_gas_main_2D
    implicit none 
    integer,parameter :: Nmax=1400
    real*8,dimension (0:Nmax,0:Nmax) ::rho0,mx0,my0,E0,rho,mx,my,E 
    real*8,dimension (1:Nmax) :: x,y
    real*8,dimension (0:Nmax,1:Nmax):: F1,F2,F3,F4
    real*8,dimension (1:Nmax,0:Nmax):: G1,G2,G3,G4
    integer i,j,Nx,Ny,newdtneed,nstep,istep
    character(40) filename,rhon,un,vn,En
    character(20) BoundC,Extype,scaleSw,shrinkSw
    real*8 pstar,lambda_max,p1,p2
    real*8 gamma,dx,dy,x0,y0,xlength,ylength,dt,Tfinal,t,theta,dtstep
    common/parameter1/dx,dy,x0,y0,xlength,ylength,dt,Tfinal,t,theta,gamma
    common/parameter2/Nx,Ny
    common/parameter3/x,y
    common/parameter4/BoundC,Extype,scaleSw,shrinkSw
    !=======input the data========================
    
    Nx=100
    Ny=100
    tfinal=2.d0
    x0=-0.2d0
    y0=0.d0
    xlength=0.4d0
    ylength=0.8d0
    
    !!RT 
    ! Nx=128
    ! Ny=1024
    ! Tfinal=2.95d0
    ! x0=0.d0
    ! y0=0.d0
    ! xlength=0.25d0
    ! ylength=1.d0

    dx=xlength/Nx
    dy=ylength/Ny
    theta=1.3d0
    gamma=1.4d0
    ! gamma=5.d0/3.d0 !!RT example
    nstep=1
    dtstep=Tfinal/nstep

    Extype='MCW'
    ! Extype='Exp'
    ! Extype='Imp'
    ! Extype='KH'
    ! Extype='RT'
    ! Extype='Cfg3'

    !! CURH Optimization to be included or excluded
    ! shrinkSw='on'
    ! scaleSw='on'
    
    shrinkSw='off'
    scaleSw='off'
    
    !!===========================initial data=====================================================
    call initialset_2D(rho0,mx0,my0,E0)
    ! call initialset_2DRT(rho0,mx0,my0,E0)
    ! call ini_Riemann(rho0,mx0,my0,E0)
    
    t=0.0d0 
    do istep=0,nstep
        Tfinal=istep*dtstep
        do while (t<Tfinal-1.0d-7)
            newdtneed=1
            call  numerical_flux_Euler2D(rho0,mx0,my0,E0,newdtneed,F1,F2,F3,F4,G1,G2,G3,G4)
            call  SSPRK_Eulergas2D(rho0,mx0,my0,E0,F1,F2,F3,F4,G1,G2,G3,G4,rho,mx,my,E)
            ! call  numerical_flux_RT(rho0,mx0,my0,E0,newdtneed,F1,F2,F3,F4,G1,G2,G3,G4)
            ! call  SSPRK_EulergasRT(rho0,mx0,my0,E0,F1,F2,F3,F4,G1,G2,G3,G4,rho,mx,my,E)
            call extendcell2d(rho,mx,my,E)
            t=t+dt
            write(*,*) 't=', t, ' | dt=', dt
            rho0=rho
            mx0=mx
            my0=my
            E0=E
        enddo
        
        write(rhon,101) istep

        open(11,file=rhon)
        !! .....Writing the solution.............................................
        do j=1,Ny
            do i=1,Nx
                write(11,*) x(i),y(j),rho0(i,j)
            enddo
        enddo
        close(11)

    enddo


101 format("MCW-rho-CURH",i3.3)
! 101 format("EXP-rho-CURH",i3.3)
! 101 format("IMP1-rho-CURH",i3.3)
! 101 format("KHI3-rho-CURH",i3.3)
! 101 format("RTI-rho-CURH",i3.3)
! 101 format("CFG3-rho-CURH",i3.3)

100 format (' ', f24.16,' ', f24.16,' ', f24.16,' ', f24.16,' ', f24.16,' ', f24.16)
end program

real*8 function dminmod(v1,v2,v3)
    implicit none
    real*8 v1,v2,v3,u
    dminmod=0.5d0*(dsign(1.d0,v1)+dsign(1.d0,v3))*dmin1(dabs(v1),dabs(v2),dabs(v3))
    return
end function
    
real*8 function dminmod4(a,b,c,d)
    implicit none
    real*8 a,b,c,d,e,f
    e = dmin1(a,b,c,d)
    f = dmax1(a,b,c,d)
    dminmod4 = 0.d0
    if(e>0.d0) then
        dminmod4 = e
    elseif(f<0.d0) then
        dminmod4 = f
    endif
    return
end function
    
real*8 function minmod(x,y)
    implicit none
    real*8 x,y
    minmod=0.5d0*(dsign(1.d0,x)+dsign(1.d0,y))*dmin1(dabs(x),dabs(y))
    return
end function

