    subroutine numerical_flux_Euler2D(rho,mx,my,E,newdtneed,F1,F2,F3,F4,G1,G2,G3,G4)
    implicit none 
    integer,parameter :: Nmax=1400
    real*8, dimension (0:Nmax,0:Nmax) :: rho,mx,my,E
    real*8,dimension (0:Nmax,1:Nmax) :: ap,am,ap2,am2
    real*8,dimension (1:Nmax,0:Nmax) :: bp,bm,bp2,bm2
    real*8,dimension (0:Nmax,1:Nmax) :: rhoE,mxE,myE,EE
    real*8,dimension (1:Nmax,0:Nmax) :: rhoN,mxN,myN,EN
    real*8,dimension (1:Nmax,1:Nmax) :: rhoW,mxW,myW,EW,rhoS,mxS,myS,ES
    real*8,dimension (0:Nmax,0:Nmax) :: rhoNE,rhoNW,rhoSE,rhoSW,mxNE,mxNW,mxSE,mxSW,myNE,myNW,mySE,mySW,ENE,ENW,ESE,ESW
    real*8,dimension (0:Nmax,1:Nmax):: F1,F2,F3,F4
    real*8,dimension (1:Nmax,0:Nmax):: G1,G2,G3,G4
    real*8, external :: minmod, dminmod,dminmod4
    real*8 pE,uE,vE,pW,uW,vW,pS,uS,vS,pN,uN,vN,f1E,f2E,f3E,f4E,f1W,f2W,f3W,f4W,g1S,g2S,g3S,g4S,g1N,g2N,g3N,g4N
    character(20) BoundC,Extype,scaleSw,shrinkSw
    real*8 differ1,differ2,differ,prod_a,prod_b,alpha,ratio1,ratio2,ratio3,ratio4,ratiomax,ratiomin
    real*8 df1,df2,df3,df4,dg1,dg2,dg3,dg4,du1,du2,du3,du4,du1eps,du2eps,du3eps,du4eps
    real*8 d1,d2,d3,slx,sly,amax,bmax,dist_a,dist_b,drho,dmx,dmy,dE,rho_star,mx_star,my_star,E_star
    real*8 eps_tol
    integer i,j,k,Nx,Ny,newdtneed
    real*8 dx,dy,x0,y0,xlength,ylength,dt,Tfinal,t,theta,gamma
    common/parameter1/dx,dy,x0,y0,xlength,ylength,dt,Tfinal,t,theta,gamma
    common/parameter2/Nx,Ny
    common/parameter4/BoundC,Extype,scaleSw,shrinkSw
    !===using the linear appriximation or 2nd order approximation=================
    do j=1,Ny
        do i=1,Nx
            !! Density PLR calculation
            ! SLX
            d1=theta*(rho(i+1,j)-rho(i,j))
            d2=0.5d0*(rho(i+1,j)-rho(i-1,j))
            d3=theta*(rho(i,j)-rho(i-1,j))  
            slx=0.5d0*dminmod(d1,d2,d3)
            ! SLY
            d1=theta*(rho(i,j+1)-rho(i,j))
            d2=0.5d0*(rho(i,j+1)-rho(i,j-1))
            d3=theta*(rho(i,j)-rho(i,j-1))
            sly=0.5d0*dminmod(d1,d2,d3)  
            ! PLR
            rhoE(i,j)=rho(i,j)+slx
            rhoW(i,j)=rho(i,j)-slx
            rhoN(i,j)=rho(i,j)+sly
            rhoS(i,j)=rho(i,j)-sly            
            ! rhoNE(i,j)=(slx+sly)+rho(i,j)
            ! rhoNW(i,j)=(-slx+sly)+rho(i,j)
            ! rhoSE(i,j)=(slx-sly)+rho(i,j)
            ! rhoSW(i,j)=(-slx-sly)+rho(i,j)
            
            !! MomentumX PLR calculation
            ! SLX
            d1=theta*(mx(i+1,j)-mx(i,j))
            d2=0.5d0*(mx(i+1,j)-mx(i-1,j))
            d3=theta*(mx(i,j)-mx(i-1,j))   
            slx=0.5d0*dminmod(d1,d2,d3)
            ! SLY            
            d1=theta*(mx(i,j+1)-mx(i,j))
            d2=0.5d0*(mx(i,j+1)-mx(i,j-1))
            d3=theta*(mx(i,j)-mx(i,j-1))
            sly=0.5d0*dminmod(d1,d2,d3)  
            ! PLR
            mxE(i,j)=mx(i,j)+slx
            mxW(i,j)=mx(i,j)-slx
            mxN(i,j)=mx(i,j)+sly
            mxS(i,j)=mx(i,j)-sly
            ! mxNE(i,j)=(slx+sly)+mx(i,j)
            ! mxNW(i,j)=(-slx+sly)+mx(i,j)
            ! mxSE(i,j)=(slx-sly)+mx(i,j)
            ! mxSW(i,j)=(-slx-sly)+mx(i,j)
            
            !! MomentumY PLR calculation
            ! SLX
            d1=theta*(my(i+1,j)-my(i,j))
            d2=0.5d0*(my(i+1,j)-my(i-1,j))
            d3=theta*(my(i,j)-my(i-1,j))   
            slx=0.5d0*dminmod(d1,d2,d3)
            ! SLY
            d1=theta*(my(i,j+1)-my(i,j))
            d2=0.5d0*(my(i,j+1)-my(i,j-1))
            d3=theta*(my(i,j)-my(i,j-1))
            sly=0.5d0*dminmod(d1,d2,d3)  
            ! PLR
            myE(i,j)=my(i,j)+slx
            myW(i,j)=my(i,j)-slx
            myN(i,j)=my(i,j)+sly
            myS(i,j)=my(i,j)-sly            
            ! myNE(i,j)=(slx+sly)+my(i,j)
            ! myNW(i,j)=(-slx+sly)+my(i,j)
            ! mySE(i,j)=(slx-sly)+my(i,j)
            ! mySW(i,j)=(-slx-sly)+my(i,j)
            
            !! Energy PLR calculation
            ! SLX
            d1=theta*(E(i+1,j)-E(i,j))
            d2=0.5d0*(E(i+1,j)-E(i-1,j))
            d3=theta*(E(i,j)-E(i-1,j))   
            slx=0.5d0*dminmod(d1,d2,d3)
            ! SLY
            d1=theta*(E(i,j+1)-E(i,j))
            d2=0.5d0*(E(i,j+1)-E(i,j-1))
            d3=theta*(E(i,j)-E(i,j-1))
            sly=0.5d0*dminmod(d1,d2,d3)  
            ! PLR
            EE(i,j)=E(i,j)+slx
            EW(i,j)=E(i,j)-slx
            EN(i,j)=E(i,j)+sly
            ES(i,j)=E(i,j)-sly
            ! ENE(i,j)=(slx+sly)+E(i,j)
            ! ENW(i,j)=(-slx+sly)+E(i,j)
            ! ESE(i,j)=(slx-sly)+E(i,j)
            ! ESW(i,j)=(-slx-sly)+E(i,j)
        enddo
    enddo

    ! Extend cells for the PLR calculated earlier
    select case (BoundC)
    case('periodic')
        !!===periodic boundary======
        do j=1,Ny
            rhoE(0,j)=rhoE(Nx,j)
            mxE(0,j)=mxE(Nx,j)
            myE(0,j)=myE(Nx,j)
            EE(0,j)=EE(Nx,j)

            rhoW(Nx+1,j)=rhoW(1,j)
            mxW(Nx+1,j)=mxW(1,j)
            myW(Nx+1,j)=myW(1,j)
            EW(Nx+1,j)=EW(1,j)
            
            rhoNE(0,j)=rhoNE(Nx,j)
            rhoNW(Nx+1,j)=rhoNW(1,j)
            rhoSE(0,j)=rhoSE(Nx,j)
            rhoSW(Nx+1,j)=rhoSW(1,j)
            
            mxNE(0,j)=mxNE(Nx,j)
            mxNW(Nx+1,j)=mxNW(1,j)
            mxSE(0,j)=mxSE(Nx,j)
            mxSW(Nx+1,j)=mxSW(1,j) 
            
            myNE(0,j)=myNE(Nx,j)
            myNW(Nx+1,j)=myNW(1,j)
            mySE(0,j)=mySE(Nx,j)
            mySW(Nx+1,j)=mySW(1,j)
            
            ENE(0,j)=ENE(Nx,j)
            ENW(Nx+1,j)=ENW(1,j)
            ESE(0,j)=ESE(Nx,j)
            ESW(Nx+1,j)=ESW(1,j) 
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoN(i,Ny)
            mxN(i,0)=mxN(i,Ny)
            myN(i,0)=myN(i,Ny)
            EN(i,0)=EN(i,Ny)
            rhoS(i,Ny+1)=rhoS(i,1)
            mxS(i,Ny+1)=mxS(i,1)
            myS(i,Ny+1)=myS(i,1)
            ES(i,Ny+1)=ES(i,1)
    
            rhoSW(i,Ny+1)=rhoSW(i,1) 
            rhoSE(i,Ny+1)=rhoSE(i,1)
            rhoNE(i,0)=rhoNE(i,Ny)
            rhoNW(i,0)=rhoNW(i,Ny)
            mxSW(i,Ny+1)=mxSW(i,1)
            rhoSE(i,Ny+1)=mxSE(i,1) 
            mxNE(i,0)=mxNE(i,Ny)
            mxNW(i,0)=mxNW(i,Ny)
            mySW(i,Ny+1)=mySW(i,1) 
            mySE(i,Ny+1)=mySE(i,1) 
            myNE(i,0)=myNE(i,Ny) 
            myNW(i,0)=myNW(i,Ny) 
            ESW(i,Ny+1)=ESW(i,1) 
            ESE(i,Ny+1)=ESE(i,1)
            ENE(i,0)=ENE(i,Ny)
            ENW(i,0)=ENW(i,Ny) 
        enddo
    case('reflectfree')
        !! ==reflecting boundary(solid wall)+free====
        do j=1,Ny
            rhoE(0,j)=rhoW(1,j)
            rhoW(Nx+1,j)=rhoE(Nx,j)

            ! rhoNE(0,j)=rhoNW(1,j)
            ! rhoNW(Nx+1,j)=rhoNE(Nx,j)
            ! rhoSE(0,j)=rhoSW(1,j)
            ! rhoSW(Nx+1,j)=rhoSE(Nx,j)

            mxE(0,j)=-mxW(1,j)
            mxW(Nx+1,j)=mxE(Nx,j)

            ! mxNE(0,j)=-mxNW(1,j)
            ! mxNW(Nx+1,j)=mxNE(Nx,j)
            ! mxSE(0,j)=-mxSW(1,j)
            ! mxSW(Nx+1,j)=mxSE(Nx,j)
            
            myE(0,j)=myW(1,j)
            myW(Nx+1,j)=myE(Nx,j)

            ! myNE(0,j)=myNW(1,j)
            ! myNW(Nx+1,j)=myNE(Nx,j)
            ! mySE(0,j)=mySW(1,j)
            ! mySW(Nx+1,j)=mySE(Nx,j)
            
            EE(0,j)=EW(1,j)
            EW(Nx+1,j)=EE(Nx,j)

            ! ENE(0,j)=ENW(1,j)
            ! ENW(Nx+1,j)=ENE(Nx,j)
            ! ESE(0,j)=ESW(1,j)
            ! ESW(Nx+1,j)=ESE(Nx,j)     
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoS(i,1)
            rhoS(i,Ny+1)=rhoN(i,Ny)
            ! rhoSW(i,Ny+1)=rhoNW(i,Ny) 
            ! rhoSE(i,Ny+1)=rhoNE(i,Ny)
            ! rhoNE(i,0)=rhoSE(i,1)
            ! rhoNW(i,0)=rhoSW(i,1)

            mxN(i,0)=mxS(i,1)
            mxS(i,Ny+1)=mxN(i,Ny)
            ! mxSW(i,Ny+1)=mxNW(i,Ny)
            ! mxSE(i,Ny+1)=mxNE(i,Ny) 
            ! mxNE(i,0)=mxSE(i,1)
            ! mxNW(i,0)=mxSW(i,1)
            
            myN(i,0)=-myS(i,1)
            myS(i,Ny+1)=myN(i,Ny)
            ! mySW(i,Ny+1)=myNW(i,Ny) 
            ! mySE(i,Ny+1)=myNE(i,Ny) 
            ! myNE(i,0)=-mySE(i,1) 
            ! myNW(i,0)=-mySW(i,1) 
            
            EN(i,0)=ES(i,1)
            ES(i,Ny+1)=EN(i,Ny)
            ! ESW(i,Ny+1)=ENW(i,Ny) 
            ! ESE(i,Ny+1)=ENE(i,Ny)
            ! ENE(i,0)=ESE(i,1)
            ! ENW(i,0)=ESW(i,1) 
        enddo
    case('free')
        ! !==free====
        do j=1,Ny
            ! extend Density PLR
            rhoE(0,j)=rhoW(1,j)
            rhoW(Nx+1,j)=rhoE(Nx,j)

            ! rhoNE(0,j)=rhoNW(1,j)
            ! rhoNW(Nx+1,j)=rhoNE(Nx,j)
            
            ! rhoSE(0,j)=rhoSW(1,j)
            ! rhoSW(Nx+1,j)=rhoSE(Nx,j)
            
            ! extend MomentumX PLR
            mxE(0,j)=mxW(1,j)
            mxW(Nx+1,j)=mxE(Nx,j)
            
            ! mxNE(0,j)=mxNW(1,j)
            ! mxNW(Nx+1,j)=mxNE(Nx,j)
            
            ! mxSE(0,j)=mxSW(1,j)
            ! mxSW(Nx+1,j)=mxSE(Nx,j)
            
            ! extend MomentumY PLR
            myE(0,j)=myW(1,j)
            myW(Nx+1,j)=myE(Nx,j)
            ! myNE(0,j)=myNW(1,j)
            ! myNW(Nx+1,j)=myNE(Nx,j)
            ! mySE(0,j)=mySW(1,j)
            ! mySW(Nx+1,j)=mySE(Nx,j)
            
            ! extend Energy PLR
            EE(0,j)=EW(1,j)
            EW(Nx+1,j)=EE(Nx,j)
            ! ENE(0,j)=ENW(1,j)
            ! ENW(Nx+1,j)=ENE(Nx,j)
            ! ESE(0,j)=ESW(1,j)
            ! ESW(Nx+1,j)=ESE(Nx,j)
        enddo
        do i=1,Nx
            ! Density PLR
            rhoN(i,0)=rhoS(i,1)
            rhoS(i,Ny+1)=rhoN(i,Ny)

            ! rhoSW(i,Ny+1)=rhoNW(i,Ny) 
            ! rhoSE(i,Ny+1)=rhoNE(i,Ny)
    
            ! rhoNE(i,0)=rhoSE(i,1)
            ! rhoNW(i,0)=rhoSW(i,1)
            
            ! MomentumX PLR
            mxN(i,0)=mxS(i,1)
            mxS(i,Ny+1)=mxN(i,Ny)
            ! mxSW(i,Ny+1)=mxNW(i,Ny)
            ! mxSE(i,Ny+1)=mxNE(i,Ny) 
            ! mxNE(i,0)=mxSE(i,1)
            ! mxNW(i,0)=mxSW(i,1)
            
            ! MomentumY PLR
            myN(i,0)=myS(i,1)
            myS(i,Ny+1)=myN(i,Ny)
            ! mySW(i,Ny+1)=myNW(i,Ny) 
            ! mySE(i,Ny+1)=myNE(i,Ny) 
            ! myNE(i,0)=mySE(i,1) 
            ! myNW(i,0)=mySW(i,1) 
            
            ! Energy PLR
            EN(i,0)=ES(i,1)            
            ES(i,Ny+1)=EN(i,Ny)        
            ! ESW(i,Ny+1)=ENW(i,Ny) 
            ! ESE(i,Ny+1)=ENE(i,Ny)
            ! ENE(i,0)=ESE(i,1)
            ! ENW(i,0)=ESW(i,1) 
        enddo  
    case('solidwall')
        !! ==reflecting boundary(solid wall)====
        do j=1,Ny
            rhoE(0,j)=rhoW(1,j)
            mxE(0,j)=-mxW(1,j)
            myE(0,j)=myW(1,j)
            EE(0,j)=EW(1,j)
            rhoW(Nx+1,j)=rhoE(Nx,j)
            mxW(Nx+1,j)=-mxE(Nx,j)
            myW(Nx+1,j)=myE(Nx,j)
            EW(Nx+1,j)=EE(Nx,j)
            rhoNE(0,j)=rhoNW(1,j)
            rhoNW(Nx+1,j)=rhoNE(Nx,j)
            rhoSE(0,j)=rhoSW(1,j)
            rhoSW(Nx+1,j)=rhoSE(Nx,j)
            mxNE(0,j)=-mxNW(1,j)
            mxNW(Nx+1,j)=-mxNE(Nx,j)
            mxSE(0,j)=-mxSW(1,j)
            mxSW(Nx+1,j)=-mxSE(Nx,j)
            myNE(0,j)=myNW(1,j)
            myNW(Nx+1,j)=myNE(Nx,j)
            mySE(0,j)=mySW(1,j)
            mySW(Nx+1,j)=mySE(Nx,j)
            ENE(0,j)=ENW(1,j)
            ENW(Nx+1,j)=ENE(Nx,j)
            ESE(0,j)=ESW(1,j)
            ESW(Nx+1,j)=ESE(Nx,j)     
        enddo
        do i=1,Nx
            rhoN(i,0)=rhoS(i,1)
            mxN(i,0)=mxS(i,1)
            myN(i,0)=-myS(i,1)
            EN(i,0)=ES(i,1)
            rhoS(i,Ny+1)=rhoN(i,Ny)
            mxS(i,Ny+1)=mxN(i,Ny)
            myS(i,Ny+1)=-myN(i,Ny)
            ES(i,Ny+1)=EN(i,Ny)
            rhoSW(i,Ny+1)=rhoNW(i,Ny) 
            rhoSE(i,Ny+1)=rhoNE(i,Ny)
            rhoNE(i,0)=rhoSE(i,1)
            rhoNW(i,0)=rhoSW(i,1)
            mxSW(i,Ny+1)=mxNW(i,Ny)
            mxSE(i,Ny+1)=mxNE(i,Ny) 
            mxNE(i,0)=mxSE(i,1)
            mxNW(i,0)=mxSW(i,1)
            mySW(i,Ny+1)=-myNW(i,Ny) 
            mySE(i,Ny+1)=-myNE(i,Ny) 
            myNE(i,0)=-mySE(i,1) 
            myNW(i,0)=-mySW(i,1) 
            ESW(i,Ny+1)=ENW(i,Ny) 
            ESE(i,Ny+1)=ENE(i,Ny)
            ENE(i,0)=ESE(i,1)
            ENW(i,0)=ESW(i,1) 
        enddo
    end select
    
    !! CURH Numerical Flux calculation starts here
    eps_tol=1.0d-12
    do j=1,Ny
        do i=0,Nx
            uE=mxE(i,j)/rhoE(i,j)
            uW=mxW(i+1,j)/rhoW(i+1,j)
            vE=myE(i,j)/rhoE(i,j)
            vW=myW(i+1,j)/rhoW(i+1,j)
            pE=(gamma-1.d0)*(-0.5d0*rhoE(i,j)*(uE**2+vE**2)+EE(i,j))
            pW=(gamma-1.d0)*(-0.5d0*rhoW(i+1,j)*(uW**2+vW**2)+EW(i+1,j))
            
            ap2(i,j)=dmax1(uW+dsqrt(gamma*pW/rhoW(i+1,j)),uE+dsqrt(gamma*pE/rhoE(i,j)),0.d0)
            am2(i,j)=dmin1(uW-dsqrt(gamma*pW/rhoW(i+1,j)),uE-dsqrt(gamma*pE/rhoE(i,j)),0.d0)

            ! CURH Optimization block
            select case (scaleSw)
            case('on')
                differ1=dabs(0.5d0*(rhoW(i+1,j)*uW**2-rhoE(i,j)*uE**2)+(pW-pE)/(gamma-1.d0))
                differ2=0.5d0*dabs(rhoW(i+1,j)*vW**2-rhoE(i,j)*vE**2)
                differ=dsqrt(differ1**2+differ2**2)
                if (differ>1.0d-12) then
                    alpha=differ1/differ
                else
                    alpha=0.0d0
                endif
                ap(i,j)=dmax1(uW+alpha*dsqrt(gamma*pW/rhoW(i+1,j)),uE+alpha*dsqrt(gamma*pE/rhoE(i,j)),0.d0)
                am(i,j)=dmin1(uW-alpha*dsqrt(gamma*pW/rhoW(i+1,j)),uE-alpha*dsqrt(gamma*pE/rhoE(i,j)),0.d0)
            case('off')
                ap(i,j)=ap2(i,j)
                am(i,j)=am2(i,j)
            end select
    
            f1E=mxE(i,j)
            f2E=mxE(i,j)*uE+pE
            f3E=mxE(i,j)*vE
            f4E=uE*(EE(i,j)+pE)

            f1W=mxW(i+1,j)
            f2W=mxW(i+1,j)*uW+pW
            f3W=vW*mxW(i+1,j)
            f4W=uW*(EW(i+1,j)+pW)
    
            ! CURH Optimization Block
            if (shrinkSw=='on') then
                df1=f1w-f1e
                df2=f2w-f2e
                df3=f3w-f3e
                df4=f4w-f4e

                du1=rhow(i+1,j)-rhoe(i,j)
                du2=mxw(i+1,j)-mxe(i,j)
                du3=myw(i+1,j)-mye(i,j)
                du4=ew(i+1,j)-ee(i,j)
                
                if(du1>0.d0) then
                    du1eps=dmax1(du1,1.d-10)
                else
                    du1eps=dmin1(du1,-1.d-10)
                endif
                if(du2>0.d0) then
                    du2eps=dmax1(du2,1.d-10)
                else
                    du2eps=dmin1(du2,-1.d-10)
                endif
                if(du3>0.d0) then
                    du3eps=dmax1(du3,1.d-10)
                else
                    du3eps=dmin1(du3,-1.d-10)
                endif
                if(du4>0.d0) then
                    du4eps=dmax1(du4,1.d-10)
                else
                    du4eps=dmin1(du4,-1.d-10)
                endif

                ratio1=2.d0*df1/(du1+du1eps)
                ratio2=2.d0*df2/(du2+du2eps)
                ratio3=2.d0*df3/(du3+du3eps)
                ratio4=2.d0*df4/(du4+du4eps)
                ratiomin=dmin1(ratio1,ratio2,ratio3,ratio4)
                ratiomax=dmax1(ratio1,ratio2,ratio3,ratio4)

                if(ratiomax>0.d0) then
                    ap(i,j)=dmin1(ap(i,j),ratiomax)
                    am(i,j)=dmax1(am(i,j),-ratiomax)
                endif
                
                if(ratiomin<0.d0) then
                    ap(i,j)=dmin1(ap(i,j),-ratiomin)
                    am(i,j)=dmax1(am(i,j),ratiomin)
                endif
            endif
    
            dist_a=ap(i,j)-am(i,j)
            prod_a=ap(i,j)*am(i,j)
            if (dist_a>eps_tol) then
                ! rho_star=((ap(i,j)*rhoW(i+1,j)-am(i,j)*rhoE(i,j))-(f1W-f1E))/dist_a
                ! mx_star=((ap(i,j)*mxW(i+1,j)-am(i,j)*mxE(i,j))-(f2W-f2E))/dist_a
                ! my_star=((ap(i,j)*myW(i+1,j)-am(i,j)*myE(i,j))-(f3W-f3E))/dist_a
                ! E_star=((ap(i,j)*EW(i+1,j)-am(i,j)*EE(i,j))-(f4W-f4E))/dist_a
                
                ! drho=dminmod4((rhoSW(i+1,j)-rho_star),(rho_star-rhoSE(i,j)),(rhoNW(i+1,j)-rho_star),(rho_star-rhoNE(i,j))) !anti-diffusion term  
                ! dmx=dminmod4((mxSW(i+1,j)-mx_star),(mx_star-mxSE(i,j)),(mxNW(i+1,j)-mx_star),(mx_star-mxNE(i,j))) !anti-diffusion term 
                ! dmy=dminmod4((mySW(i+1,j)-my_star),(my_star-mySE(i,j)),(myNW(i+1,j)-my_star),(my_star-myNE(i,j))) !anti-diffusion term 
                ! dE=dminmod4((ESW(i+1,j)-E_star),(E_star-ESE(i,j)),(ENW(i+1,j)-E_star),(E_star-ENE(i,j))) !anti-diffusion term  
    
                ! F1(i,j)=((ap(i,j)*f1E-am(i,j)*f1W)+prod_a*((rhoW(i+1,j)-rhoE(i,j))-drho))/dist_a
                ! F2(i,j)=((ap(i,j)*f2E-am(i,j)*f2W)+prod_a*((mxW(i+1,j)-mxE(i,j))-dmx))/dist_a
                ! F3(i,j)=((ap(i,j)*f3E-am(i,j)*f3W)+prod_a*((myW(i+1,j)-myE(i,j))-dmy))/dist_a
                ! F4(i,j)=((ap(i,j)*f4E-am(i,j)*f4W)+prod_a*((EW(i+1,j)-EE(i,j))-dE))/dist_a

                F1(i,j)=((ap(i,j)*f1E-am(i,j)*f1W)+prod_a*((rhoW(i+1,j)-rhoE(i,j))))/dist_a
                F2(i,j)=((ap(i,j)*f2E-am(i,j)*f2W)+prod_a*((mxW(i+1,j)-mxE(i,j))))/dist_a
                F3(i,j)=((ap(i,j)*f3E-am(i,j)*f3W)+prod_a*((myW(i+1,j)-myE(i,j))))/dist_a
                F4(i,j)=((ap(i,j)*f4E-am(i,j)*f4W)+prod_a*((EW(i+1,j)-EE(i,j))))/dist_a
            else
                F1(i,j)=0.5d0*(f1E+f1W)
                F2(i,j)=0.5d0*(f2E+f2W) 
                F3(i,j)=0.5d0*(f3E+f3W)
                F4(i,j)=0.5d0*(f4E+f4W)
            endif  
        enddo
    enddo
    
    do j=0,Ny
        do i=1,Nx
            uN=mxN(i,j)/rhoN(i,j)
            uS=mxS(i,j+1)/rhoS(i,j+1)
            vN=myN(i,j)/rhoN(i,j)
            vS=myS(i,j+1)/rhoS(i,j+1)
            pN=(gamma-1.d0)*(-0.5d0*rhoN(i,j)*(uN**2+vN**2)+EN(i,j))
            pS=(gamma-1.d0)*(-0.5d0*rhoS(i,j+1)*(uS**2+vS**2)+ES(i,j+1))
           
            bp2(i,j)=dmax1(vS+dsqrt(gamma*pS/rhoS(i,j+1)),vN+dsqrt(gamma*pN/rhoN(i,j)),0.d0)
            bm2(i,j)=dmin1(vS-dsqrt(gamma*pS/rhoS(i,j+1)),vN-dsqrt(gamma*pN/rhoN(i,j)),0.d0)

            ! CURH Optimization Block
            select case (scaleSw)
            case('on')
                differ1=dabs(0.5d0*(rhoS(i,j+1)*vS**2-rhoN(i,j)*vN**2)+(pS-pN)/(gamma-1.d0))
                differ2=0.5d0*dabs(rhoS(i,j+1)*uS**2-rhoN(i,j)*uN**2)
                differ=dsqrt(differ1**2+differ2**2)
                if (differ>1.0d-12) then
                    alpha=differ1/differ
                else
                    alpha=0.0d0
                endif
                bp(i,j)=dmax1(vS+alpha*dsqrt(gamma*pS/rhoS(i,j+1)),vN+alpha*dsqrt(gamma*pN/rhoN(i,j)),0.d0)
                bm(i,j)=dmin1(vS-alpha*dsqrt(gamma*pS/rhoS(i,j+1)),vN-alpha*dsqrt(gamma*pN/rhoN(i,j)),0.d0)      
            case ('off')
                bp(i,j)=bp2(i,j)
                bm(i,j)=bm2(i,j)
            end select

            g1N=myN(i,j)
            g2N=myN(i,j)*uN
            g3N=myN(i,j)*vN+pN
            g4N=vN*(EN(i,j)+pN)

            g1S=myS(i,j+1)
            g2S=myS(i,j+1)*uS
            g3S=myS(i,j+1)*vS+pS
            g4S=vS*(ES(i,j+1)+pS)

            ! CURH Optimization Block
            if (shrinkSw=='on') then
                dg1=g1s-g1n
                dg2=g2s-g2n
                dg3=g3s-g3n
                dg4=g4s-g4n

                du1=rhos(i,j+1)-rhon(i,j)
                du2=mxs(i,j+1)-mxn(i,j)
                du3=mys(i,j+1)-myn(i,j)
                du4=es(i,j+1)-en(i,j)
                
                if(du1>0.d0) then
                    du1eps=dmax1(du1,1.d-10)
                else
                    du1eps=dmin1(du1,-1.d-10)
                endif
                if(du2>0.d0) then
                    du2eps=dmax1(du2,1.d-10)
                else
                    du2eps=dmin1(du2,-1.d-10)
                endif
                if(du3>0.d0) then
                    du3eps=dmax1(du3,1.d-10)
                else
                    du3eps=dmin1(du3,-1.d-10)
                endif
                if(du4>0.d0) then
                    du4eps=dmax1(du4,1.d-10)
                else
                    du4eps=dmin1(du4,-1.d-10)
                endif
                ratio1=2.d0*dg1/(du1+du1eps)
                ratio2=2.d0*dg2/(du2+du2eps)
                ratio3=2.d0*dg3/(du3+du3eps)
                ratio4=2.d0*dg4/(du4+du4eps)
                ratiomin=dmin1(ratio1,ratio2,ratio3,ratio4)
                ratiomax=dmax1(ratio1,ratio2,ratio3,ratio4)
                if(ratiomax>0.d0) then
                    bp(i,j)=dmin1(bp(i,j),ratiomax)
                    bm(i,j)=dmax1(bm(i,j),-ratiomax)
                endif
                if(ratiomin<0.d0) then
                    bp(i,j)=dmin1(bp(i,j),-ratiomin)
                    bm(i,j)=dmax1(bm(i,j),ratiomin)
                endif
            endif


            dist_b=bp(i,j)-bm(i,j)
            prod_b=bp(i,j)*bm(i,j)
            if (dist_b>eps_tol) then
                ! rho_star=((bp(i,j)*rhoS(i,j+1)-bm(i,j)*rhoN(i,j))-(g1S-g1N))/dist_b
                ! mx_star=((bp(i,j)*mxS(i,j+1)-bm(i,j)*mxN(i,j))-(g2S-g2N))/dist_b
                ! my_star=((bp(i,j)*myS(i,j+1)-bm(i,j)*myN(i,j))-(g3S-g3N))/dist_b
                ! E_star=((bp(i,j)*ES(i,j+1)-bm(i,j)*EN(i,j))-(g4S-g4N))/dist_b
                
                ! drho=dminmod4((rhoSW(i,j+1)-rho_star),(rho_star-rhoNW(i,j)),(rhoSE(i,j+1)-rho_star),(rho_star-rhoNE(i,j))) !anti-diffusion term  
                ! dmx=dminmod4((mxSW(i,j+1)-mx_star),(mx_star-mxNW(i,j)),(mxSE(i,j+1)-mx_star),(mx_star-mxNE(i,j))) !anti-diffusion term 
                ! dmy=dminmod4((mySW(i,j+1)-my_star),(my_star-myNW(i,j)),(mySE(i,j+1)-my_star),(my_star-myNE(i,j))) !anti-diffusion term 
                ! dE=dminmod4((ESW(i,j+1)-E_star),(E_star-ENW(i,j)),(ESE(i,j+1)-E_star),(E_star-ENE(i,j))) !anti-diffusion term 
    
                ! G1(i,j)=((bp(i,j)*g1N-bm(i,j)*g1S)+prod_b*((rhoS(i,j+1)-rhoN(i,j))-drho))/dist_b
                ! G2(i,j)=((bp(i,j)*g2N-bm(i,j)*g2S)+prod_b*((mxS(i,j+1)-mxN(i,j))-dmx))/dist_b
                ! G3(i,j)=((bp(i,j)*g3N-bm(i,j)*g3S)+prod_b*((myS(i,j+1)-myN(i,j))-dmy))/dist_b
                ! G4(i,j)=((bp(i,j)*g4N-bm(i,j)*g4S)+prod_b*((ES(i,j+1)-EN(i,j))-dE))/dist_b

                G1(i,j)=((bp(i,j)*g1N-bm(i,j)*g1S)+prod_b*((rhoS(i,j+1)-rhoN(i,j))))/dist_b
                G2(i,j)=((bp(i,j)*g2N-bm(i,j)*g2S)+prod_b*((mxS(i,j+1)-mxN(i,j))))/dist_b
                G3(i,j)=((bp(i,j)*g3N-bm(i,j)*g3S)+prod_b*((myS(i,j+1)-myN(i,j))))/dist_b
                G4(i,j)=((bp(i,j)*g4N-bm(i,j)*g4S)+prod_b*((ES(i,j+1)-EN(i,j))))/dist_b
            else
                G1(i,j)=0.5d0*(g1N+g1S)
                G2(i,j)=0.5d0*(g2N+g2S) 
                G3(i,j)=0.5d0*(g3N+g3S)
                G4(i,j)=0.5d0*(g4N+g4S)
            endif 
        enddo
    enddo
    
    !CFL
    if (newdtneed==1) then 
        amax=0.d0
        do j=1,Ny
            do i=0,Nx
                amax=dmax1(amax,ap2(i,j),-am2(i,j))
            enddo
        enddo
        bmax=0.d0
        do j=0,Ny
            do i=1,Nx
                bmax=dmax1(bmax,bp2(i,j),-bm2(i,j)) 
            enddo
        enddo
        dt=0.475d0*dmin1(dx/amax,dy/bmax)
        if (t+dt>Tfinal)  then
            dt=Tfinal-t
        endif
    endif
    
    return
    end subroutine
    
    
