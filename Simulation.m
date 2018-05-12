classdef Simulation < matlab.System
    % Sinulates gyroscope procession

    % Public, tunable properties
    properties
%         a; b; c; d; da; db; dc; dd;
        r; R0; Ri; H; h; L; mr; p;      % Physical prorerties
        a_0; b_0; c_0; d_0;             % Initial angular offests
        dA_0; dB_0; dC_0; dD_0;         % Initial angular velocities
        X1; X2; X3; X4; X5; X6; X7; X8; % State vector parameters
        X;                              % State vector
        R10; R21; R32; R43;
        H1; H2; H3; H4;

    end

    methods(Access = public)
        
        function init(obj)
            syms t a(t) b(t) c(t) d(t) X1 X2 X3 X4 X5 X6 X7 X8;
            obj.X = [X1 X2 X3 X4 X5 X6 X7 X8];
            obj.X1 = a;
            obj.X2 = b;
            obj.X3 = c;
            obj.X4 = d;
            obj.X5 = diff(a);
            obj.X6 = diff(b);
            obj.X7 = diff(c);
            obj.X8 = diff(d);
            obj.R10 = [cos(obj.X1) -sin(obj.X1) 0; sin(obj.X1) cos(obj.X1) 0; 0 0 1];%Rotation matrices
            obj.R21 = [1 0 0; 0 cos(obj.X2) -sin(obj.X2); 0 sin(obj.X2) cos(obj.X2)];
            obj.R32 = [cos(obj.X3) -sin(obj.X3) 0; sin(obj.X3) cos(obj.X3) 0; 0 0 1];
            obj.R43 = [cos(obj.X4) -sin(obj.X4) 0; sin(obj.X4) cos(obj.X4) 0; 0 0 1];

        end
        
        
        function setPhysicalParams(obj, r_, R0_, Ri_, H_, h_, L_, mr_, p_)
            % Sets the physical parameters of the sysyem
            obj.r = r_;
            obj.R0 = R0_;
            obj.Ri = Ri_;
            obj.H = H_;
            obj.h = h_;
            obj.L = L_;
            obj.mr = mr_;
            obj.p = p_;
        end
        
        
        function setInitialConditions(obj, a0, b0, c0, d0, dA0, dB0, dC0, dD0)
            obj.a_0 = a0; 
            obj.b_0 = b0; 
            obj.c_0 = c0; 
            obj.d_0 = d0;     
            obj.dA_0 = dA0; 
            obj.dB_0 = dB0; 
            obj.dC_0 = dC0; 
            obj.dD_0 = dD0; 
        end
        
        
        function generateEOMs(obj)
            % Generates H1-4 which will be used for numeric integration.
            %Relative angular velocites of system
            g = 9.81;
            syms FGx FGy FGz FOx FOy FOz MGx MGy t;
            syms t  a(t) b(t) c(t) d(t);
            
            fourw43 = [0; 0; diff(d)];
            threew32 = [0; 0; diff(c)];
            twow21 = [diff(b); 0; 0];
            onew10 = [0; 0; diff(a)];


            %Position vector from point O to the centre of mass
            threerOG = [0; 0; obj.L];

            %Absolute angluar velocites of system
            threew3 = threew32 + transpose(obj.R32)*twow21 + ...
                transpose(obj.R32)*transpose(obj.R21)*onew10;
            fourw4 = fourw43 + transpose(obj.R43)*threew3;

            %Moment of inertia of the rotor about centre of mass in frame three
            threeIRG = obj.mr*[(1/12)*(3*obj.Ri^2 + obj.h^2) 0 0; 0 (1/12)*(3*obj.Ri^2 + obj.h^2) 0; 0 0 (1/2)*obj.Ri^2];

            %Moment of inertia of the rotor about centre of mass in frame four
            fourIRG = transpose(obj.R43)*threeIRG*obj.R43;

            %Mass of elements of the system
            mrod = obj.p*pi*(obj.r^2)*obj.H;
            %mtorus = 2*rho*(PI^2)*r^2*R0;
            mframe = mrod + 0.023; %2*mtorus;
            mR = obj.mr - mrod;
            %Gravitational forces
            zeroGR = [0; 0; -mR*g];
            zeroGf = [0; 0; -mframe*g];

            %Reaction moments and forces
            threeMG = [MGx; MGy; 0];
            threeFG = [FGx; FGy; FGz];
            threeFO = [FOx; FOy; FOz];

            %Caclulating the moments of inertia of the frame
            IaG = mrod*[(1/12)*(3*obj.r^2 + obj.H^2) 0 0; 0 (1/12)*(3*obj.r^2 + obj.H^2) 0; 0 0 (1/2)*obj.r^2];
            parallel = [obj.L^2 0 0; 0 obj.L^2 0; 0 0 0];
            IaO = IaG + mrod*parallel;
            IbG = mframe*[(5/8)*obj.r^2 + (1/2)*obj.R0^2 0 0; 0 (5/8)*obj.r^2 + (1/2)*obj.R0^2 0; 0 0 (3/4)*obj.r^2 + obj.R0^2];
            IbO = IbG + mframe*parallel;
            IcG = mframe*[(5/8)*obj.r^2 + (1/2)*obj.R0^2 0 0; 0 (3/4)*obj.r^2 + obj.R0^2 0; 0 0 (5/8)*obj.r^2 + (1/2)*obj.R0^2];
            IcO = IcG + mframe*parallel;

            %Moments of inertia of the frame about point G and O
            IfG = IaG + IbG + IcG;
            IfO = IaO + IbO + IcO;

            %First derivative and second derivative vectors of position vector rOG
            threerOGdot = cross((threew3), (threerOG));
            threerOGdotdot = diff(threerOGdot, t) + cross((threew3), (threerOGdot));

            %Mass Newton-Euler equations for the rotor
            threeFR = threeFG + transpose(obj.R32)*transpose(obj.R21)*transpose(obj.R10)*zeroGR;
            sumthreeFR = mR*threerOGdotdot;

            %Mass Newton-Euler equations for the frame
            threeFf = -threeFG + threeFO + transpose(obj.R32)*transpose(obj.R21)*transpose(obj.R10)*zeroGf;
            sumthreeFf = mframe*threerOGdotdot;

            %Angular momentum of the rotor
            fourhRG = fourIRG*fourw4;
            threehRG = obj.R43*fourhRG;

            %Moment Newton-Euler equations for the rotor
            threeMRG = threeMG;
            sumthreeMRG = diff(threehRG) + cross(formula(threew3), formula(threehRG));

            %Angular mementum of the frame
            threehfO = IfO*threew3;

            %Moment Newton-Euler equations for the frame
            threeMfO = -threeMG + cross((threerOG), (-threeFG + transpose(obj.R32)*transpose(obj.R21)*transpose(obj.R10)*zeroGf));
            sumthreeMfO = diff(threehfO) + cross((threew3), (threehfO));

           %Equations of motion
            g1 = formula(sumthreeMRG);
            g1 = g1(3);

            g2 = formula(sumthreeMfO);
            g2 = g2(3);

            temp1 = formula(threeFR);
            temp2 = formula(sumthreeFR);
            FGx = solve(temp1(1) == temp2(2), FGx);
            FGy = solve(temp1(2) == temp2(2), FGy);

            temp1 = formula(sumthreeMRG);
            MGx = temp1(1);
            MGy = temp1(2);

            temp1 = formula(sumthreeMfO - threeMfO);
            temp1 = subs(temp1);

            g3 = temp1(1);
            g4 = temp1(2);

            g1 = simplify(g1);
            g2 = simplify(g2);
            g3 = simplify(g3);
            g4 = simplify(g4);


            syms A B C D d_A d_B d_C d_D ddA ddB ddC ddD

            G1 = subs(g1, [a b c d diff(a) diff(b) diff(c) diff(d) diff(diff(a)) diff(diff(b)) diff(diff(c)) diff(diff(d))], [A B C D d_A d_B d_C d_D ddA ddB ddC ddD]);
            G2 = subs(g2, [a b c d diff(a) diff(b) diff(c) diff(d) diff(diff(a)) diff(diff(b)) diff(diff(c)) diff(diff(d))], [A B C D d_A d_B d_C d_D ddA ddB ddC ddD]);
            G3 = subs(g3, [a b c d diff(a) diff(b) diff(c) diff(d) diff(diff(a)) diff(diff(b)) diff(diff(c)) diff(diff(d))], [A B C D d_A d_B d_C d_D ddA ddB ddC ddD]);
            G4 = subs(g4, [a b c d diff(a) diff(b) diff(c) diff(d) diff(diff(a)) diff(diff(b)) diff(diff(c)) diff(diff(d))], [A B C D d_A d_B d_C d_D ddA ddB ddC ddD]);

            [X,Y] = equationsToMatrix([G1 == 0; G2 == 0; G3 == 0; G4 == 0], [ddA ddB ddC ddD]);

            %Substitute actual values we measured here
            %t = 1;

            %X = subs(X);
            %Y = subs(Y);

            H = linsolve(X,Y);
            obj.H1 = @(A, B, C, D, d_A, d_B, d_C, d_D) (subs(H(1)));
            obj.H2 = @(A, B, C, D, d_A, d_B, d_C, d_D) (subs(H(2)));
            obj.H3 = @(A, B, C, D, d_A, d_B, d_C, d_D) (subs(H(3)));
            obj.H4 = @(A, B, C, D, d_A, d_B, d_C, d_D) (subs(H(4)));
        end
               
        
        function [alpha, beta, gamma, delta] = simulateGyroscope(obj, fps, tEnd)
           % returns vecors of values for alpha, beta, gamma, and delta
           % corresponding to discrete time states specified by the fps and
           % simulation end time tEnd
           
           ddt = DiffX;     % derivitave calculator for numerical integration
           ddt.giveH(obj.H1, obj.H2, obj.H3, obj.H4);
           
           % initialize output vectors for efficency
           alpha = zeros(1, tEnd*fps);
           beta  = zeros(1, tEnd*fps);
           gamma = zeros(1, tEnd*fps);
           delta = zeros(1, tEnd*fps);
           % generate ode model
           obj.X = [obj.a_0, obj.b_0, obj.c_0, obj.d_0, obj.dA_0, obj.dB_0, obj.dC_0, obj.dD_0];
           [T, Y] = ode45(@ddt.diffX,  [0 tEnd], obj.X);
           
           % produce time vector
           t = linspace(0,tEnd,tEnd*fps);
            
           % compute state at every frame
           for i = 1:numel(t)
               % intrapolate values of A B C and D
               for j = 1:(numel(T)-1)
                   if T(j+1)>t(i)
                           break
                   end
               end
               alpha(i) = Y(j,1) + t(i)*(Y(j+1,5)-Y(j,5))/(T(j+1)-T(j)); %alpha at time step
               beta(i) = Y(j,2) + t(i)*(Y(j+1,6)-Y(j,6))/(T(j+1)-T(j)); %beta at time step
               gamma(i) = Y(j,3) + t(i)*(Y(j+1,7)-Y(j,7))/(T(j+1)-T(j)); %gamma at time step
               delta(i) = Y(j,4) + t(i)*(Y(j+1,8)-Y(j,8))/(T(j+1)-T(j)); %delta at time step
            end
        
            
        end
    end
end
