    subroutine initialset_1D(rho,m1,E)
    implicit none
    integer,parameter :: Nmax=10000
    real*8,dimension (0:Nmax,1) :: rho,m1,E,p,u
    real*8,dimension (1:Nmax,1) :: x
    character(20) BoundC,Extype
    integer i,Nx
    real*8 dx,x0,xlength,dt,Tfinal,t,theta,gamma,xW,xE
    common/parameter1/dx,x0,xlength,dt,Tfinal,t,theta,gamma
    common/parameter2/Nx
    common/parameter3/x
    common/parameter4/BoundC,Extype

    do i=1,Nx
        x(i,1)=x0+(i-0.5d0)*dx
    enddo

    select case(Extype)
    case ('MCW')
        !!=========Moving Contact Wave=========
        BoundC='free'
        do i=1,Nx
            if(x(i,1)<0.3d0) then
                rho(i,1)=1.4d0
            else
                rho(i,1)=1.d0
            endif
            m1(i,1)=0.1d0*rho(i,1)
            E(i,1)=1.d0/(gamma-1.d0)+0.5d0*rho(i,1)*0.01d0
            !m1(i,1)=0.d0
            !E(i,1)=1.d0/(gamma-1.d0)
        enddo
    case ('Lax')
        BoundC='free'
        do i=1,Nx
            if(x(i,1)<0.d0) then
                rho(i,1)=0.445d0
                u(i,1)=0.698d0
                p(i,1)=3.528d0
            else
                rho(i,1)=0.5d0
                u(i,1)=0.0d0
                p(i,1)=0.571d0
            endif
            m1(i,1)=rho(i,1)*u(i,1)
            E(i,1)=p(i,1)/(gamma-1.d0)+0.5d0*rho(i,1)*u(i,1)**2
        enddo
    case ('Blast')
        BoundC='solidwall'
        do i=1,Nx
            if(x(i,1)<0.1d0) then
                rho(i,1)=1.0d0
                u(i,1)=0.d0
                p(i,1)=1000.d0
            elseif (x(i,1)<0.9d0) then
                rho(i,1)=1.d0
                u(i,1)=0.0d0
                p(i,1)=0.01d0
            else
                rho(i,1)=1.d0
                u(i,1)=0.0d0
                p(i,1)=100.d0
            endif
            m1(i,1)=rho(i,1)*u(i,1)
            E(i,1)=p(i,1)/(gamma-1.d0)+0.5d0*rho(i,1)*u(i,1)**2
        enddo
    case ('ShuOS')
        BoundC='free'
        do i=1,Nx
            if(x(i,1)<-4.d0) then
                rho(i,1)=3.857143d0
                u(i,1)=2.629369d0
                p(i,1)=10.3333d0
            else
                rho(i,1)=1.d0+0.2d0*dsin(5.d0*x(i,1))
                u(i,1)=0.0d0
                p(i,1)=1.d0
            endif
            m1(i,1)=rho(i,1)*u(i,1)
            E(i,1)=p(i,1)/(gamma-1.d0)+0.5d0*rho(i,1)*u(i,1)**2
        enddo
        case ('test')
        BoundC='free'
        do i=1,Nx
            if(x(i,1)<0.4d0) then
                rho(i,1)=5.99924d0
                u(i,1)=19.5975d0
                p(i,1)=460.894d0
            else
                rho(i,1)=5.99242d0
                u(i,1)=-6.19633d0
                p(i,1)=46.0950d0
            endif
            m1(i,1)=rho(i,1)*u(i,1)
            E(i,1)=p(i,1)/(gamma-1.d0)+0.5d0*m1(i,1)*u(i,1)
        enddo
         case ('SMS')
        BoundC='free'
        do i=1,Nx
            if(x(i,1)<0.1d0) then
                rho(i,1)=3.86d0
                m1(i,1)=-3.1266d0
                E(i,1)=27.0913d0
            else
                rho(i,1)=1.0d0
                m1(i,1)=-3.44d0
                E(i,1)=8.4168d0
            endif
        enddo
    case ('SOD')
        BoundC='free'
        do i=1,Nx
            if(x(i,1)<0.5d0) then
                rho(i,1)=1.d0
                u(i,1)=0.d0
                p(i,1)=1.d0
            else
                rho(i,1)=0.125d0
                u(i,1)=0.d0
                p(i,1)=0.1d0
            endif
            m1(i,1)=rho(i,1)*u(i,1)
            E(i,1)=p(i,1)/(gamma-1.d0)+0.5d0*m1(i,1)*u(i,1)
        enddo
    case ('SPP')
        BoundC='free'
        do i=1,Nx
            if(x(i,1)<0.3d0) then
                rho(i,1)=1.d0
                u(i,1)=0.75d0
                p(i,1)=1.d0
            else
                rho(i,1)=0.125d0
                u(i,1)=0.d0
                p(i,1)=0.1d0
            endif
            m1(i,1)=rho(i,1)*u(i,1)
            E(i,1)=p(i,1)/(gamma-1.d0)+0.5d0*m1(i,1)*u(i,1)
        enddo
    end select

    !select case(Extype)
    !case ('MCW')
    !    !!=========Moving Contact Wave=========
    !    BoundC='free'
    !    do i=1,Nx
    !        xW=x(i,1)-0.5d0*dx
    !        xE=x(i,1)+0.5d0*dx
    !        if(xE<0.3d0) then
    !            rho(i,1)=1.4d0
    !        elseif (xE>0.3d0 .and. xW<0.3d0) then
    !            rho(i,1)=1.d0/dx*(1.4d0*(0.3d0-xW)+1.d0*(xE-0.3d0))
    !        else
    !            rho(i,1)=1.d0
    !        endif
    !        m1(i,1)=0.1d0*rho(i,1)
    !        E(i,1)=1.d0/(gamma-1.d0)+0.5d0*rho(i,1)*0.01d0
    !    enddo
    !case ('SCSR')
    !    BoundC='free'
    !    do i=1,Nx
    !        xW=x(i,1)-0.5d0*dx
    !        xE=x(i,1)+0.5d0*dx
    !        if(xE<0.8d0) then
    !            p(i,1)=1000.d0
    !        elseif (xE>0.8d0 .and. xW<0.8d0) then
    !            p(i,1)=1.d0/dx*(1000.d0*(0.8d0-xW)+0.01d0*(xE-0.8d0))
    !        else
    !            p(i,1)=0.01d0
    !        endif
    !        rho(i,1)=1.d0
    !        u(i,1)=-19.59745d0
    !        m1(i,1)=rho(i,1)*u(i,1)
    !        E(i,1)=p(i,1)/(gamma-1.d0)+0.5d0*rho(i,1)*u(i,1)**2
    !    enddo
    !case ('ShuOS')
    !    !!=========Moving Contact Wave=========
    !    BoundC='free'
    !    do i=1,Nx
    !        xW=x(i,1)-0.5d0*dx
    !        xE=x(i,1)+0.5d0*dx
    !        if (xE<0.0) then
    !            rho(i,1)=3.857143d0
    !            m1(i,1)=rho(i,1)*(-0.920279d0)
    !            E(i,1)=10.33333d0/(gamma-1.d0)+0.5d0*m1(i,1)**2/rho(i,1) !for E
    !        elseif (xW<0.d0 .and. (xE>0.d0 .and. xE<10.d0)) then
    !            rho(i,1)=(3.857143d0*(0.d0-xW)+xE-0.04d0*dcos(5.d0*xE)+0.04d0)/dx!for rho
    !            m1(i,1)=(3.857143d0*(-0.920279d0)*(0.d0-xW)-3.549648d0*(xE-0.04d0*dcos(5.d0*xE)+0.04d0))/dx !for m=rho*u
    !            E(i,1)=1.d0/dx*(10.33333d0/(gamma-1.d0)+0.5d0*3.857143d0*(0.920279d0**2))*(0.d0-xW) &
    !                & +1.d0/dx*(1.d0/(gamma-1.d0)*xE+0.5d0*(3.549648d0**2)*(xE-0.04d0*(dcos(5.d0*xE)-1.d0))) !for E
    !        elseif ((xW>0.d0).and.(xE<10.d0)) then
    !            rho(i,1)=(dx-0.04d0*(dcos(5.d0*xE)-dcos(5.d0*xW)))/dx !for rho
    !            m1(i,1)=-3.549648d0*(dx-0.04d0*(dcos(5.d0*xE)-dcos(5.d0*xW)))/dx !for m=rho*u
    !            E(i,1)=1.d0/(gamma-1.d0)+0.5d0*(-3.549648d0)**2*(dx-0.04d0*(dcos(5.d0*xE)-dcos(5.d0*xW)))/dx !for E
    !        elseif ((xW>0.d0 .and. xW<10.d0).and.(xE>10.d0)) then
    !            rho(i,1)=1.d0-0.04d0*1.d0/dx*(dcos(50.d0)-dcos(5.d0*xW))!for rho
    !            m1(i,1)=(1.d0-0.04d0*1.d0/dx*(dcos(50.d0)-dcos(5.d0*xW)))*(-3.549648d0) !for m=rho*u
    !            E(i,1)=1.d0/(gamma-1.d0)+0.5d0*(-3.549648d0)**2*(1.d0-0.04d0*1.d0/dx*(dcos(50.d0)-dcos(5.d0*xW))) !for E
    !        else
    !            rho(i,1)=1.d0 !for rho
    !            m1(i,1)=-3.549648d0  !for m=rho*u
    !            E(i,1)=1.d0/(gamma-1.d0)+0.5d0*1.d0*((-3.549648d0)**2) !for E
    !        endif
    !    enddo
    !end select

    call extendcell1D(rho,m1,E)

    return
    end subroutine

