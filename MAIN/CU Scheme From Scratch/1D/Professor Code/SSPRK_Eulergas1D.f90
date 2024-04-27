    subroutine SSPRK_Eulergas1D(rho0,mx0,E0,F1,F2,F3,rhonew,mxnew,Enew)
    implicit none
    integer,parameter :: Nmax=10000
    real*8,dimension (0:Nmax,1) :: rho0,mx0,E0,rho1,mx1,E1,rho2,mx2,E2,rhonew,mxnew,Enew
    real*8,dimension (0:Nmax,1):: F1,F2,F3
    integer i,Nx,newdtneed
    real*8 dx,x0,xlength,dt,Tfinal,t,theta,gamma,lambda
    common/parameter1/dx,x0,xlength,dt,Tfinal,t,theta,gamma
    common/parameter2/Nx

    lambda=dt/dx

    do i=1,Nx
        rho1(i,1)=rho0(i,1)-lambda*(F1(i,1)-F1(i-1,1))
        mx1(i,1)=mx0(i,1)-lambda*(F2(i,1)-F2(i-1,1))
        E1(i,1)=E0(i,1)-lambda*(F3(i,1)-F3(i-1,1))
    enddo

    call extendcell1D(rho1,mx1,E1)
    newdtneed=0
    call numerical_flux_Euler1D(rho1,mx1,E1,newdtneed,F1,F2,F3)
    do i=1,Nx
        rho2(i,1)=(3.d0*rho0(i,1)+rho1(i,1)-lambda*(F1(i,1)-F1(i-1,1)))/4.d0
        mx2(i,1)=(3.d0*mx0(i,1)+mx1(i,1)-lambda*(F2(i,1)-F2(i-1,1)))/4.d0
        E2(i,1)=(3.d0*E0(i,1)+E1(i,1)-lambda*(F3(i,1)-F3(i-1,1)))/4.d0
    enddo

    call extendcell1D(rho2,mx2,E2)
    newdtneed=0
    call numerical_flux_Euler1D(rho2,mx2,E2,newdtneed,F1,F2,F3)
    do i=1,Nx
        rhonew(i,1)=(rho0(i,1)+2.d0*rho2(i,1)-2.d0*lambda*(F1(i,1)-F1(i-1,1)))/3.d0
        mxnew(i,1)=(mx0(i,1)+2.d0*mx2(i,1)-2.d0*lambda*(F2(i,1)-F2(i-1,1)))/3.d0
        Enew(i,1)=(E0(i,1)+2.d0*E2(i,1)-2.d0*lambda*(F3(i,1)-F3(i-1,1)))/3.d0
    enddo

    return
    end subroutine
