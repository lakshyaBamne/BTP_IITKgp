    subroutine numerical_flux_Euler1D(rho,mx,E,newdtneed,F1,F2,F3)
    implicit none
    integer,parameter :: Nmax=10000
    real*8, dimension (0:Nmax,1) :: rho,mx,E
    real*8,dimension (0:Nmax,1) :: ap,am,ap2,am2
    real*8,dimension (0:Nmax,1) :: rhoE,mxE,EE
    real*8,dimension (1:Nmax,1) :: rhoW,mxW,EW
    real*8,dimension (0:Nmax,1) :: F1,F2,F3,F4
    real*8, external :: minmod, dminmod
    real*8 pE,uE,pW,uW,f1E,f2E,f3E,f1W,f2W,f3W
    character(20) BoundC,Extype,shrinkSw
    real*8 prod_a,ratio1,ratio2,ratio3,ratiomax,ratiomin
    real*8 df1,df2,df3,du1,du2,du3,du1eps,du2eps,du3eps
    real*8 d1,d2,d3,slx,amax,dist_a,drho,dmx,dE,rho_star,mx_star,E_star
    real*8 eps_tol
    integer i,Nx,newdtneed
    real*8 dx,x0,xlength,dt,Tfinal,t,theta,gamma
    common/parameter1/dx,x0,xlength,dt,Tfinal,t,theta,gamma
    common/parameter2/Nx
    common/parameter4/BoundC,Extype,shrinkSw
    !===using the linear appriximation or 2nd order approximation=================

    !! Calculating the Piecewise Linear Reconstruction
    do i=1,Nx
        d1=theta*(rho(i+1,1)-rho(i,1))
        d2=0.5d0*(rho(i+1,1)-rho(i-1,1))
        d3=theta*(rho(i,1)-rho(i-1,1))
        slx=0.5d0*dminmod(d1,d2,d3)
        rhoE(i,1)=rho(i,1)+slx
        rhoW(i,1)=rho(i,1)-slx

        d1=theta*(mx(i+1,1)-mx(i,1))
        d2=0.5d0*(mx(i+1,1)-mx(i-1,1))
        d3=theta*(mx(i,1)-mx(i-1,1))
        slx=0.5d0*dminmod(d1,d2,d3)
        mxE(i,1)=mx(i,1)+slx
        mxW(i,1)=mx(i,1)-slx

        d1=theta*(E(i+1,1)-E(i,1))
        d2=0.5d0*(E(i+1,1)-E(i-1,1))
        d3=theta*(E(i,1)-E(i-1,1))
        slx=0.5d0*dminmod(d1,d2,d3)
        EE(i,1)=E(i,1)+slx
        EW(i,1)=E(i,1)-slx
    enddo

    !! Extend the cells for the Piecewise Linear Reconstructions calculated earlier
    select case (BoundC)
    case('periodic')
        !!===periodic boundary======
        rhoE(0,1)=rhoE(Nx,1)
        mxE(0,1)=mxE(Nx,1)
        EE(0,1)=EE(Nx,1)
        rhoW(Nx+1,1)=rhoW(1,1)
        mxW(Nx+1,1)=mxW(1,1)
        EW(Nx+1,1)=EW(1,1)
    case('reflectfree')
        !! ==reflecting boundary(solid wall)+free====
        rhoE(0,1)=rhoW(1,1)
        mxE(0,1)=-mxW(1,1)
        EE(0,1)=EW(1,1)
        rhoW(Nx+1,1)=rhoE(Nx,1)
        mxW(Nx+1,1)=mxE(Nx,1)
        EW(Nx+1,1)=EE(Nx,1)
    case('free')
        ! !==free====
        rhoE(0,1)=rhoW(1,1)
        rhoW(Nx+1,1)=rhoE(Nx,1)

        mxE(0,1)=mxW(1,1)
        mxW(Nx+1,1)=mxE(Nx,1)
        
        EE(0,1)=EW(1,1)
        EW(Nx+1,1)=EE(Nx,1)
    case('solidwall')
        !! ==reflecting boundary(solid wall)====
        rhoE(0,1)=rhoW(1,1)
        rhoW(Nx+1,1)=rhoE(Nx,1)
        mxE(0,1)=-mxW(1,1)
        mxW(Nx+1,1)=-mxE(Nx,1)
        EE(0,1)=EW(1,1)
        EW(Nx+1,1)=EE(Nx,1)
    end select

    eps_tol=1.0d-12
    do i=0,Nx

        !! calculating primitive variables which will be used later on
        uE=mxE(i,1)/rhoE(i,1)
        uW=mxW(i+1,1)/rhoW(i+1,1)
        pE=(gamma-1.d0)*(-0.5d0*rhoE(i,1)*uE**2+EE(i,1))
        pW=(gamma-1.d0)*(-0.5d0*rhoW(i+1,1)*uW**2+EW(i+1,1))
        
        ap(i,1)=dmax1(uW+dsqrt(gamma*pW/rhoW(i+1,1)),uE+dsqrt(gamma*pE/rhoE(i,1)),0.d0)
        am(i,1)=dmin1(uW-dsqrt(gamma*pW/rhoW(i+1,1)),uE-dsqrt(gamma*pE/rhoE(i,1)),0.d0)
        
        ap2(i,1)=ap(i,1)
        am2(i,1)=am(i,1)

        !! calculating flux variables which will be used later on
        f1E=mxE(i,1)
        f2E=mxE(i,1)*uE+pE
        f3E=uE*(EE(i,1)+pE)
        f1W=mxW(i+1,1)
        f2W=mxW(i+1,1)*uW+pW
        f3W=uW*(EW(i+1,1)+pW)

        !! CURH optimization block
        if (shrinkSw=='on') then
            df1=f1W-f1E
            du1=rhoW(i+1,1)-rhoE(i,1)
            df2=f2W-f2E
            du2=mxW(i+1,1)-mxE(i,1)
            df3=f3W-f3E
            du3=EW(i+1,1)-EE(i,1)
            if(du1>-1.0d-10) then
                du1eps=dmax1(du1,1.d-10)
            else
                du1eps=dmin1(du1,-1.d-10)
            endif
            if(du2>-1.0d-10) then
                du2eps=dmax1(du2,1.d-10)
            else
                du2eps=dmin1(du2,-1.d-10)
            endif
            if(du3>-1.0d-10) then
                du3eps=dmax1(du3,1.d-10)
            else
                du3eps=dmin1(du3,-1.d-10)
            endif
            ratio1=2.d0*df1/(du1+du1eps)
            ratio2=2.d0*df2/(du2+du2eps)
            ratio3=2.d0*df3/(du3+du3eps)
            ratiomin=dmin1(ratio1,ratio2,ratio3)
            ratiomax=dmax1(ratio1,ratio2,ratio3)
          
            if(ratiomax>-1.0d-10) then
                ap(i,1)=dmin1(ap2(i,1),ratiomax)
                am(i,1)=dmax1(am2(i,1),-ratiomax)
            endif
            if(ratiomin<1.0d-10) then
                ap(i,1)=dmin1(ap2(i,1),-ratiomin)
                am(i,1)=dmax1(am2(i,1),ratiomin)
            endif
        endif

        !! calculating the Antidiffusion terms and the final CU Numerical Flux
        dist_a=ap(i,1)-am(i,1)
        prod_a=ap(i,1)*am(i,1)
        if (dist_a>eps_tol) then
            ! rho_star=((ap(i,1)*rhoW(i+1,1)-am(i,1)*rhoE(i,1))-(f1W-f1E))/dist_a
            ! mx_star=((ap(i,1)*mxW(i+1,1)-am(i,1)*mxE(i,1))-(f2W-f2E))/dist_a
            ! E_star=((ap(i,1)*EW(i+1,1)-am(i,1)*EE(i,1))-(f3W-f3E))/dist_a
            ! drho=minmod(rhoW(i+1,1)-rho_star,rho_star-rhoE(i,1)) !anti-diffusion term
            ! dmx=minmod(mxW(i+1,1)-mx_star,mx_star-mxE(i,1)) !anti-diffusion term
            ! dE=minmod(EW(i+1,1)-E_star,E_star-EE(i,1)) !anti-diffusion term

            ! F1(i,1)=(ap(i,1)*f1E-am(i,1)*f1W+prod_a*(rhoW(i+1,1)-rhoE(i,1)-drho))/dist_a
            ! F2(i,1)=(ap(i,1)*f2E-am(i,1)*f2W+prod_a*(mxW(i+1,1)-mxE(i,1)-dmx))/dist_a
            ! F3(i,1)=(ap(i,1)*f3E-am(i,1)*f3W+prod_a*(EW(i+1,1)-EE(i,1)-dE))/dist_a

            ! F1(i,1)=(ap(i,1)*f1E-am(i,1)*f1W+prod_a*(rhoW(i+1,1)-rhoE(i,1)))/dist_a
            ! F2(i,1)=(ap(i,1)*f2E-am(i,1)*f2W+prod_a*(mxW(i+1,1)-mxE(i,1)))/dist_a
            ! F3(i,1)=(ap(i,1)*f3E-am(i,1)*f3W+prod_a*(EW(i+1,1)-EE(i,1)))/dist_a

            F1(i,1)=(ap(i,1)*f1E-am(i,1)*f1W+prod_a*(rhoW(i+1,1)-rhoE(i,1)))/dist_a
            F2(i,1)=(ap(i,1)*f2E-am(i,1)*f2W+prod_a*(mxW(i+1,1)-mxE(i,1)))/dist_a
            F3(i,1)=(ap(i,1)*f3E-am(i,1)*f3W+prod_a*(EW(i+1,1)-EE(i,1)))/dist_a
        else
            F1(i,1)=0.5d0*(f1E+f1W)
            F2(i,1)=0.5d0*(f2E+f2W)
            F3(i,1)=0.5d0*(f3E+f3W)
        endif
    enddo


    !CFL
    if (newdtneed==1) then
        amax=0.d0
        do i=0,Nx
            amax=dmax1(amax,ap2(i,1),-am2(i,1))
        enddo
        dt=0.475d0*dx/amax
        if (t+dt>Tfinal)  then
            dt=Tfinal-t
        endif
    endif

    return
    end subroutine


