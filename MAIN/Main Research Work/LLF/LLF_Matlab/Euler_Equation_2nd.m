function Euler_Equation_2nd
% Burgers_Equation
% ====================================================

x0               = -5      ;
x1               = 15      ;

N                = 201    ;  % number of points;
Gamma            = 1.4    ;
Final_Time       = 5.0   ;
CFL              = 0.45   ;

Movie_Resolution = 0.001    ;
% ====================================================
NV               = 3      ;
% ====================================================
% Initial Condition Setup
% Problem_Case = 1 , Shu-Osher Problem
% Problem_Case = 2 , Sod Problem
% Problem_Case = 3 , Lax Problem
Problem_Case            = 6      ;
[x0, x1, N, Final_Time, Movie_Resolution, fname] = Problem_Setup(Problem_Case);
% ====================================================
% LF Setup : 0 = LLF ; 1 = CIDRA
Case_LF       = 0      ;
CIDRA_Epsilon = 1e-10  ;
CIDRA_String  = 'CIDRA';
if (Case_LF == 0), CIDRA_String  = 'LF'; end
% ====================================================
% WENO Setup
WENO_Order       = 5     ;
WENO_Flag        = 1     ;
WENO_Epsilon     = 1e-12 ;
WENO_Power       = 2     ;
% ====================================================
% Index Setup
r  = (WENO_Order + 1)/2  ;
N0 = 1      	         ;
N2 = N0 + r              ;
N1 = N2 - 1              ;
N3 = N2 + (N-1) - 1      ;
N4 = N3 + 1              ;
N5 = N3 + r              ;
% ====================================================
alpha_LLF = zeros(N5-N0+1,2);
[ x , dx ] = Domain_Setup   ( N       ) ;

% CIDRA_Epsilon = dx^2  ;

Q  =  Initial_Condition (    x )    ;
Q0 =  Q                             ; 
Q  = Boundary_Condition ( Q )       ;

Step = 0; Time = 0; dt = Time_Step ( Q ) ;

% SIMPLE VERSION BELOW (Comment/Uncomment between the * lines)
% **************************BEGIN****SIMPLE VERSION***********************
% while (Time < Final_Time)
%     dt   = Time_Step ( Q )   ;
%     Step = Step + 1          ;
%     Time = Time + dt         ;
%     
%     Q = Runge_Kutta ( Q ) ;
%     
%     fprintf('Plot ... \n') ; plot(x, Q(:,1), 'o-'); grid on
%     title(['t = ', num2str(Time, '%1.2f')]);
% %     axis([-1,1,-(Amplitude+0.5),Amplitude+0.5]); 
%     pause(0.01)
% end
% ***********************END HERE*****************************************

% NOT SO SIMPLE VERSION BELOW (Uncomment between the * lines)
% **************************BEGIN* NOT SO SIMPLE VERSION******************
% ====================================================
disp(['Estimate Number of Time Step : ', num2str(fix((Final_Time-Time)/dt))])
disp(' Total_Step    Step         dt        Time        CPU Time')
disp('--------------------------------------------------------------------------')
formatstr = '  %7d %7d %13.4E %13.4E %13.4E %13.4E %13.4E';
fprintf(formatstr, fix((Final_Time-Time)/dt), Step, dt, Time, 0, 0)
fprintf('\n')
% ====================================================
Next_Save_Time = Time +  Movie_Resolution;
CPU_Begin = cputime;
while (Time < Final_Time)
    dt   = Time_Step ( Q )   ;
    Stability_Check
    Step = Step + 1 ; Time = Time + dt ;

    Save_Indicator = 0;
    if (Time >= Next_Save_Time)
        dt   = dt - (Time-Next_Save_Time);
        Time = Next_Save_Time;
        Save_Indicator = 1;
        Next_Save_Time = Next_Save_Time + Movie_Resolution;
    end

    if Time > Final_Time, dt = dt - (Time - Final_Time); Time = Final_Time; end

    CPU_Start = cputime           ;
    Q         = Runge_Kutta ( Q ) ;
    CPU_End   = cputime           ;
% Save_Indicator = 1;
    if (Save_Indicator)
        subplot(2,1,1)
        if Problem_Case == 6
            Rho = Q(:,1);
            U   = Q(:,2)./Rho;
            E   = Q(:,3);
            
            P   = (Gamma - 1.0)*(E - 0.5*Rho.*U.^2);
            
            Entropy = log(P./Rho.^Gamma);
            plot(x, Entropy,'-o')
            ylim([-0.15,0.15])
            %             save Entropy_LLF.txt Entropy -ascii
        else
            
            fprintf('Plot ... \n') ; plot(x, Q(:,1), '-o'); grid on
        end
        title(['t = ', num2str(Time, '%1.2f')]);
        xlim([x0,x1]);
        %         axis([x0,x1,0,1]);
        subplot(2,1,2)
        fprintf('Plot ... \n') ; plot(x, alpha_LLF(:,1), 'b-o',x, alpha_LLF(:,2), 'r-s'); grid on
        legend('LLF','CIDRA')
        title(['t = ', num2str(Time, '%1.2f')]);
        %         axis([x0,x1,0,1]);
        xlim([x0,x1]);
        
        pause(0.001)
    end

    fprintf(formatstr, fix((Final_Time-Time)/dt), Step, dt, ...
        Time, CPU_End-CPU_Start, CPU_End-CPU_Begin)
    fprintf('\n')
end
% ***********************END HERE****************************************

if Problem_Case == 6
    Rho = Q(:,1);
    U   = Q(:,2)./Rho;
    E   = Q(:,3);

    P   = (Gamma - 1.0)*(E - 0.5*Rho.*U.^2);

    Entropy = log(P./Rho.^Gamma);
    En =[x', Entropy];
    save([fname, '_', CIDRA_String,'_', num2str(N-1),'_Entropy_C2','.txt'], 'En',  '-ascii')
end
% ==================================
Out_Mat = [x', Q(:,1), Q(:,2)./Q(:,1), Q(:,3), alpha_LLF ];
save([fname, '_', CIDRA_String,'_', num2str(N-1),'_SISC.txt'], 'Out_Mat',  '-ascii')
fprintf('Done ... \n');


% ================================================================
    function [x, dx] = Domain_Setup ( N )
        
        dx = (x1-x0)/(N-1);
        x  = x0 + ((N0:N5)-N2)*dx + 1/2*dx  ;
        
    end
% ================================================================
    function Q = Runge_Kutta ( Q )
        
        %  Stage 1
        t = Time                                ;
        
        D_F = D_Fluxes( Q )                     ;
        Q1  = Q + dt * D_F                      ;

        Q1  = Boundary_Condition ( Q1 )     	;

       % Stage 2
        t = Time - dt/2                         ;
        
        D_F = D_Fluxes( Q1 )                    ;
        Q1  = (3.0*Q + Q1 + dt*D_F)/4.0         ;
        
        Q1  = Boundary_Condition ( Q1 )         ;
        
        % Stage 3
        t = Time + dt                           ;
        
        D_F = D_Fluxes( Q1 )                    ;
        Q   = (Q + 2*Q1 + 2*dt*D_F)/3.0         ;

        Q   = Boundary_Condition ( Q )	        ;
    end
% ================================================================
    function dt = Time_Step ( Q )
        Lambda = MaxEigenVector( Q );
        dt = CFL * dx/Lambda;
    end
% ================================================================
    function D_F = D_Fluxes ( Q )
        D_F = zeros(size(Q));
        
        F_Half = WENO_Flux ( Q );
        
        D_F(N2:N3,:) = -(F_Half(N2:N3,:) - F_Half(N2-1:N3-1,:))/dx;
        
    end
% ================================================================
    function F_Half = WENO_Flux ( Q )
        Q_Pc = zeros(     NV,1 ); Q_Mc = zeros(     NV,1 );
        Q_P  = zeros(N5-N0+1,NV); Q_M  = zeros(N5-N0+1,NV);
        Qc   = zeros(      4,NV);
        F_Half  = zeros(N5-N0+1,NV); 
        
        i0 = -1; i1 = 2;
        
        for i = N2-1:N3
            
            Qi     =  Q(i+i0:i+i1,:)     ;
            Qc(:,1) = Qi(:,1)            ;
            Qc(:,2) = Qi(:,2)./Qi(:,1)   ;
            Qc(:,3) = (Gamma-1)*(Qi(:,3)-1/2*Qc(:,1).*Qc(:,2).^2)   ;
            for k = 1:NV
                Q_Pc(k) = C2_P( Qc(:,k) );
                Q_Mc(k) = C2_M( Qc(:,k) );
            end
            
            Q_M(i,1) = Q_Mc(1);
            Q_M(i,2) = Q_Mc(1)*Q_Mc(2);
            Q_M(i,3) = Q_Mc(3)/(Gamma-1) + 1/2*Q_Mc(1)*Q_Mc(2)^2;
            Q_P(i,1) = Q_Pc(1);
            Q_P(i,2) = Q_Pc(1)*Q_Pc(2);
            Q_P(i,3) = Q_Pc(3)/(Gamma-1) + 1/2*Q_Pc(1)*Q_Pc(2)^2;
            [F_Half(i,:),alpha_LLF(i,:)] = LF_Flux( Q_P(i,:), Q_M(i,:) );
        end
    end
% ================================================================
    function [F,alpha] = LF_Flux( U_l, U_r )
        alpha = zeros(1,2);
        F_l = Flux( U_l );
        F_r = Flux( U_r );

        if Case_LF == 0  % LF
            Alpha = LLF( U_r, U_l );
            alpha(1) = Alpha;
        else   % CIDRA
            Alpha = LLF( U_r, U_l );
            alpha(1) = Alpha;
            Alpha = CIDRA( U_r, U_l, F_r, F_l );
            alpha(2) = Alpha;
        end
        F = 1/2*( F_r + F_l - Alpha*( U_r - U_l ) );
    end
%   ================================================================
    function Alpha = LLF( U_r, U_l )
        lambda_max_r = MaxEigenVector( U_r );
        lambda_max_l = MaxEigenVector( U_l );
    	Alpha = max( lambda_max_r, lambda_max_l );
    end


 function Alpha = LLF_M2( U_r, U_l )
        [u_r_abs, c_r] = MaxEigenVector_M( U_r );
        [u_l_abs, c_l] = MaxEigenVector_M( U_l );
        a = abs(u_r_abs/c_r - u_l_abs/c_l);
        d = a;
        t = 1.0;
        if (a<=t), d = (a^2+t^2)/(2*t); end
        d = min(d,1);
        %c = sqrt((Gamma-1)/Gamma);
        c = 1.0;
        Alpha = max( u_r_abs + c*d*c_r, u_l_abs + c*d*c_l );
%     	Alpha = max( u_r_abs, u_l_abs ) + ...
%             sqrt((Gamma-1)/Gamma)*max(c_r, c_l);
    end
%  ================================================================
    function  Alpha = CIDRA(U_r, U_l, F_r, F_l)
        Delta_U = abs( U_r(:,3) - U_l(:,3) );
        u_r = U_r(:,2) ./ U_r(:,1);
        u_l = U_l(:,2) ./ U_l(:,1);
        lambda_min = min( abs(u_l), abs(u_r) );
        lambda_max = LLF_M2( U_r, U_l );
        
        Delta_F = abs(F_r(:,3) - F_l(:,3));
        beta = 2*Delta_F./(Delta_U + max(Delta_U, CIDRA_Epsilon));
        Alpha = MaxMin( lambda_min, lambda_max, beta );
    end
% ================================================================
    function d = MaxMin( a, b, c )
        %  a <= c <= b
        d = max(a, min( b, c ));
    end
% ================================================================
    function Lambda = MaxEigenVector( Q )
        Rho = Q(:,1);
        U   = Q(:,2)./Rho;
        E   = Q(:,3);

        P   = (Gamma - 1.0)*(E - 0.5*Rho.*U.^2);

        C   = sqrt( abs( Gamma*P./Rho ) );
        Lambda = max(abs(U) + C);
    end

  function [U_abs, C] = MaxEigenVector_M( Q )
%         eps = 1e-10;
%         Rho = Q(:,1);
%         Rho2 = Rho.^2 + max(Rho.^2,eps);
%         U   = (2*Q(:,2).*Q(:,1))./Rho2;
%         E   = Q(:,3);
% 
%         P   = (Gamma - 1.0)*(E - 0.5*Rho.*U.^2);
% 
%         C   = sqrt( (2*Gamma*P.*Rho)./Rho2 );
%         Lambda = max(abs(U) + C);
        Rho = Q(:,1);
        U   = Q(:,2)./Rho;
        E   = Q(:,3);

        P   = (Gamma - 1.0)*(E - 0.5*Rho.*U.^2);

        C   = sqrt( abs( Gamma*P./Rho ) );
        U_abs = abs(U);
    end
% ================================================================
    function F = Flux( Q )
        
        F = zeros(size(Q));
        
        Rho = Q(:,1);
        U   = Q(:,2)./Rho;
        E   = Q(:,3);

        P   = (Gamma - 1.0)*(E - 0.5*Rho.*U.^2);

        F(:,1) = Rho.*U;
        F(:,2) = Rho.*U.^2 + P;
        F(:,3) = (E + P).*U;
    end
% ================================================================
    function c = MinMod(a,b)
        c = 1/2*(sign(a)+sign(b))*min(abs(a),abs(b));
    end
% ================================================================
    function fh = C2_P( f )
        a = f(2) - f(1);
        b = f(3) - f(2);
        fh = f(2) + 1/2*MinMod(a,b);
    end
% ================================================================
    function fh = C2_M( f )
        f = f(:);
        fh = C2_P( f(end:-1:1) );
    end
% ================================================================
    function [x0, x1, N, Final_Time, Movie_Resolution, fname] = Problem_Setup(Problem_Case)
        switch Problem_Case
            case (1)
                x0 = -5; 
                x1 = 15;
                N = 401;
                Final_Time = 5.0;
                Movie_Resolution = 0.1;
                fname = 'Shu_Osher';
            case (2)
                x0 = -5;
                x1 =  5;
                N  = 201;
                Final_Time = 2;
                Movie_Resolution = 0.1;
                fname = 'Sod';
            case (3)
                x0 = -5;
                x1 =  5;
                N = 201;
                Final_Time = 1.3;
                Movie_Resolution = 0.1;
                fname = 'Lax';
            case (4)
                x0 = -5;
                x1 =  5;
                N = 401;
                Final_Time = 1;
                Movie_Resolution = 0.1;
                fname = '123';
            case (5)
                x0 = 0;
                x1 = 1;
                N = 1601;
                Final_Time = 0.038;
                Movie_Resolution = 0.001;
                fname = 'Blast';
            case (6)
                x0 = -5;
                x1 = 5;
                N = 1601;
                Final_Time = 5;
                Movie_Resolution = 0.1;
                fname = 'Shock_Entropy';
        	case (7)
                x0 = 0;
                x1 = 1;
                N = 201;
                Final_Time = 2;
                Movie_Resolution = 0.05;
                fname = 'Slowly_moving_Contact';
            case (8)
                x0 = 0;
                x1 = 1;
                N = 201;
                Final_Time = 4;
                Movie_Resolution = 0.05;
                fname = 'Slowly_moving_Shock';
            case (9)
                x0 = 0;
                x1 = 1;
                N = 201;
                Final_Time = 0.1;
                Movie_Resolution = 0.01;
                fname = 'Test1';
            case (10)
                x0 = 0;
                x1 = 1;
                N = 201;
                Final_Time = 0.1;
                Movie_Resolution = 0.01;
                fname = 'Test2';
            case (11)
                x0 = 0;
                x1 = 1;
                N = 2001;
                Final_Time = 0.1;
                Movie_Resolution = 0.01;
                fname = 'Test3';
            case (12)
                x0 = 0;
                x1 = 1;
                N = 201;
                Final_Time = 0.013;
                Movie_Resolution = 0.001;
                fname = 'Test4';
            case (13)
                x0 = 0;
                x1 = 1;
                N = 201;
                Final_Time = 0.05;
                Movie_Resolution = 0.01;
                fname = 'Test5';
            case (14)
                x0 = 0;
                x1 = 1;
                N = 201;
                Final_Time = 0.005;
                Movie_Resolution = 0.0001;
                fname = 'Test6';
            case (15)
                x0 = 0;
                x1 = 1;
                N = 201;
                Final_Time = 0.1;
                Movie_Resolution = 0.01;
                fname = 'Test7';  
          
        end
    end
% ================================================================
     function Q = Initial_Condition( x )
        Q = zeros(N5-N0+1,NV);
        switch  Problem_Case
            case (1) % Shu-Osher problem
                eps = 0.2; xl = -4; k = 5;
                
                Rho = 27/7*ones(size(x)) ; 
                U = 4*sqrt(35)/9*ones(size(x)) ; 
                P = 31/3*ones(size(x));
                
                Rho(x>=xl) = 1 + eps*sin(k*x(x>=xl));
                U(x>=xl) = 0;
                P(x>=xl) = 1;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (2) % SOD
                % eps = 0.2;  k = 5;
                xl = 0;
                
                Rho = 0.125*ones(size(x)) ;
                U = 0*ones(size(x)) ;
                P = 0.1*ones(size(x));
                
                Rho(x>=xl) = 1 ;
                U(x>=xl) = 0;
                P(x>=xl) = 1;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (3) % Lax
                % eps = 0.2;  k = 5;
                xl = 0;
                
                Rho = 0.445*ones(size(x)) ;
                U = 0.698*ones(size(x)) ;
                P = 3.528*ones(size(x));
                
                Rho(x>=xl) = 0.5 ;
                U(x>=xl) = 0;
                P(x>=xl) = 0.5710;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (4) % 123
                % eps = 0.2;  k = 5;
                xl = 0;
                
                Rho = 1*ones(size(x)) ;
                U = -2*ones(size(x)) ;
                P = 0.4*ones(size(x));
                
                Rho(x>=xl) = 1 ;
                U(x>=xl) = 2;
                P(x>=xl) = 0.4;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
             case (5) % Blast Wave
                % eps = 0.2;  k = 5;
                x_1 = 0.1;
                x_2 = 0.9;
                
                Rho = 1*ones(size(x)) ;
                U = 0*ones(size(x)) ;
                P = 1000*ones(size(x));
                
                Rho(x>x_1) = 1 ;
                U(x>x_1) = 0;
                P(x>x_1) = 0.01;
                
                Rho(x>x_2) = 1 ;
                U(x>x_2) = 0;
                P(x>x_2) = 100;
                
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (6) % Shock_Entropy problem
                eps = 1/10; xl = -4.5; k = 20;
                
                Rho = 1.51695*ones(size(x)) ; 
                U = 0.523346*ones(size(x)) ; 
                P = 1.805*ones(size(x));
                
                Rho(x>=xl) = 1 + eps*sin(k*x(x>=xl));
                U(x>=xl) = 0;
                P(x>=xl) = 1;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (7) % Shock_Entropy problem
                xl = 0.3;
                
                Rho =1.4*ones(size(x)) ;
                U =  0.1*ones(size(x)) ;
                P =  1.0*ones(size(x));
                
                Rho(x>=xl) = 1 ;
                U(x>=xl) = 0.1;
                P(x>=xl) = 1;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (8)
                xl = 0.3;
                
                Rho = 3.86*ones(size(x)) ;
                M = -3.1266*ones(size(x)) ;
                E = 27.0913*ones(size(x));
                
                Rho(x>=xl) = 1 ;
                M(x>=xl) = -3.44;
                E(x>=xl) = 8.4168;
            case (9) % Test1
                xl = 0.5;
                
                Rho =1.0*ones(size(x)) ;
                U =  0.0*ones(size(x)) ;
                P =  1.0*ones(size(x));
                
                Rho(x>=xl) = 1 ;
                U(x>=xl) = 0;
                P(x>=xl) = 0.1;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (10) % Test2
                xl = 0.5;
                
                Rho =1.0*ones(size(x)) ;
                U =  0.0*ones(size(x)) ;
                P =  1.0*ones(size(x));
                
                Rho(x>=xl) = 0.125 ;
                U(x>=xl) = 0;
                P(x>=xl) = 0.1;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (11) % Test3
                xl = 0.5;
                
                Rho =1.0*ones(size(x)) ;
                U =  0.0*ones(size(x)) ;
                P =  1.0*ones(size(x));
                
                Rho(x>=xl) = 0.001 ;
                U(x>=xl) = 0;
                P(x>=xl) = 0.8;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (12) % Test4
                xl = 0.5;
                
                Rho =1.0*ones(size(x)) ;
                U =  0.0*ones(size(x)) ;
                P =  0.01*ones(size(x));
                
                Rho(x>=xl) = 1 ;
                U(x>=xl) = 0;
                P(x>=xl) = 1000;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (13) % Test5
                xl = 0.5;
                
                Rho =6.0*ones(size(x)) ;
                U =  8.0*ones(size(x)) ;
                P =  460*ones(size(x));
                
                Rho(x>=xl) = 6 ;
                U(x>=xl) = -6;
                P(x>=xl) = 46;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (14) % Test6
                xl = 0.5;
                
                Rho =600.0*ones(size(x)) ;
                U =  80.0*ones(size(x)) ;
                P =  4600*ones(size(x));
                
                Rho(x>=xl) = 6 ;
                U(x>=xl) = -6;
                P(x>=xl) = 46;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
            case (15) % Test7
                xl = 0.5;
                
                Rho =1.0*ones(size(x)) ;
                U = -2.0*ones(size(x)) ;
                P =  0.4*ones(size(x));
                
                Rho(x>=xl) = 1 ;
                U(x>=xl) = 2;
                P(x>=xl) = 0.4;
                E = P/(Gamma - 1.0) + 0.5*Rho.*U.^2;
        end
        
        if Problem_Case == 8
            Q(:,1) = Rho;
            Q(:,2) = M;
            Q(:,3) = E;
        else
          Q(:,1) = Rho;
          Q(:,2) = Rho.*U;
          Q(:,3) = E;
        end
          
    end
    
% ================================================================
     function f = Boundary_Condition( f )
        switch Problem_Case
            case ({1,6})
                  f(N0:N1,:) = Q0(N0:N1,:);
                  f(N4:N5,:) = Q0(N4:N5,:);
            case (5)
                %    Blast Wave
                L1 = N1-N0; L2 = N5-N4;
                f(N0:N1,:) =  f(N2+L1:-1:N2,:);
                f(N0:N1,2) = -f(N2+L1:-1:N2,2);
                
                f(N4:N5,:) =  f(N3:-1:N3-L2,:);
                f(N4:N5,2) = -f(N3:-1:N3-L2,2);
            otherwise
                for i = N0:N1
                    f(i,:) = f(N2,:);
                end
                for i = N4:N5
                    f(i,:) = f(N3,:);
                end
        end
    end
% ================================================================
    function Stability_Check
        if dt < 1.0e-16, error('=========  STOP ! Unstable  =============='); end
    end
% ================================================================
end
