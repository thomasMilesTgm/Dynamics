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
        
        
        function generateEOMs_G(obj)
            
            syms a(t) b(t) c(t) d(t);
            syms as bs cs ds ad bd  dd add bdd cdd ddd;
            syms r R0 Ri H h L mr p;
            syms m Dx Dy Dz
            syms Fgx Fgy Fgz Mgx Mgy Fox Foy Foz
            cd = sym('cd');
            
            SYMB = [r R0 Ri H h L mr p];
            NUMB = [obj.r, obj.R0, obj.Ri, obj.H, obj.h, obj.L, obj.mr, obj.p];
            
            COMP = [a(t) b(t) c(t) d(t) diff(a(t)) diff(b(t)) diff(c(t)) diff(d(t)) diff(diff(a(t))) diff(diff(b(t))) diff(diff(c(t))) diff(diff(d(t)))];
            READ = [as bs cs ds ad bd cd dd add bdd cdd ddd];
            
            %///////////////////////////Question 23///////////////////////////////////%
            
            IrG_3 = mr*[(3*Ri^2+h^2)/12 0 0; 0 (3*Ri^2+h^2)/12 0; 0 0 Ri^2/2];
            
            IrG_3_num = double(subs(IrG_3,SYMB, NUMB));
            
            %///////////////////////////Question 24///////////////////////////////////%
            
            IfpG_3 = H*pi*r^2*p*[(3*r^2+H^2)/12 0 0; 0 (3*r^2+H^2)/12 0; 0 0 r^2/2];
            
            IfhG_3 = 2*p*pi^2*r^2*R0*[(5/8)*r^2+(1/2)*R0^2 0 0; 0 (5/8)*r^2+(1/2)*R0^2 0; 0 0 (3/4)*r^2+R0^2];
            
            IfvG_3 = 2*p*pi^2*r^2*R0*[(5/8)*r^2+(1/2)*R0^2 0 0; 0 (3/4)*r^2+R0^2 0; 0 0 (5/8)*r^2+(1/2)*R0^2];
            
            IfG_3 = IfpG_3+IfhG_3+IfvG_3;
            
            IfG_3_num = double(subs(IfG_3,SYMB, NUMB));
            
            %///////////////////////////Question 25///////////////////////////////////%
            
            PA_na = m*[Dy^2+Dz^2 -Dx*Dy -Dx*Dz; -Dx*Dy Dx^2+Dz^2 -Dy*Dz; -Dx*Dz -Dy*Dz Dx^2+Dy^2];
            
            PS_na = [m Dx Dy Dz];
            PS_G = [H*pi*r^2*p+4*p*pi^2*r^2*R0 0 0 L];
            
            PA_G = subs(PA_na, PS_na, PS_G);
            
            IfO_3 = IfG_3 + PA_G;
            
            IfO_3_num = double(subs(IfO_3,SYMB, NUMB));
            
            %///////////////////////////Question 26///////////////////////////////////%
            
            R01 = [cos(a(t)) -sin(a(t)) 0 ; sin(a(t)) cos(a(t)) 0; 0 0 1];
            R12 = [1 0 0; 0 cos(b(t)) -sin(b(t)); 0 sin(b(t)) cos(b(t))];
            R23 = [cos(c(t)) -sin(c(t)) 0 ; sin(c(t)) cos(c(t)) 0; 0 0 1];
            R34 = [cos(d(t)) -sin(d(t)) 0 ; sin(d(t)) cos(d(t)) 0; 0 0 1];
            
            w10_1 = [0; 0; diff(a(t))];
            w21_2 = [diff(b(t)); 0; 0];
            w32_3 = [0; 0; diff(c(t))];
            w43_4 = [0; 0; diff(d(t))];
            
            w30_3 = R23.'*R12.'*w10_1 + R23.'*w21_2 + w32_3;
            w40_3 = R23.'*R12.'*w10_1 + R23.'*w21_2 + w32_3 + R34*w43_4;
            
            wd30_3 = diff(w30_3);
            wd40_3 = diff(w40_3);
            
            rOG_3 = [0; 0; L];
            
            rOG_0 = R01*R12*R23*rOG_3;
            
            rddOG_3 = simplify(expand(R23.'*R12.'*R01.'*diff(diff(rOG_0))));
            
            rddOG_3_read = subs(rddOG_3, COMP, READ);
            
            Gr_0 = [0; 0; -9.8*mr];
            
            Gf_0 = [0; 0; -9.8*(H*pi*r^2*p+4*p*pi^2*r^2*R0)];
            
            
            eqs1 = [Fgx; Fgy; Fgz] == mr*rddOG_3 - R23.'*R12.'*R01.'*Gr_0;
            
            eqs2 = [Mgx; Mgy; 0] == simplify(expand(IrG_3*wd40_3+cross(w30_3,(IrG_3*w40_3))));
            
            eqs3 = [Fox; Foy; Foz] == (H*pi*r^2*p+4*p*pi^2*r^2*R0)*rddOG_3 - R23.'*R12.'*R01.'*Gf_0 + [Fgx; Fgy; Fgz];
            
            eqs4 = [0; 0; 0] == simplify(expand([Mgx; Mgy; 0] + cross(rOG_3, [Fgx; Fgy; Fgz]) - cross(rOG_3, R23.'*R12.'*R01.'*Gf_0) + IfO_3*wd30_3 + cross(w30_3, IfO_3*w30_3)));
            
            %/////////////////////////////////////////////////////////////////////////%
            %/////////////////////////////////LAB/////////////////////////////////////%
            %/////////////////////////////////////////////////////////////////////////%
            
            %Save the equations into a large matrix and sub in non diff terms to
            %prevent fuckery
            eqs = [eqs1; eqs2; eqs3; eqs4];
            eqs_READ = subs(eqs, COMP, READ);
            
            %Define variables to be eliminated and convert to matricies
            vars = [Fgx Fgy Fgz Mgx Mgy Fox Foy Foz];
            fprintf ("HELP")
            eliminate(eqs, vars)
            [MAT, SOL] = equationsToMatrix(eqs_READ,vars);
            
            %Use this to reduce the number of equations (type MAT in command to see why
            %this arrangement
            EQN = [0; 0; 0; 0] == [SOL(6); SOL(12); SOL(11)+L*SOL(1)+SOL(5); SOL(10)-L*SOL(2)+SOL(4)];
            
            %Solve the equations and sub back in diff values
            [A B C D] = solve(EQN, [add bdd cdd ddd]);
            syms d_A d_B d_C d_D a_lpha b_eta g_amma d_elta FUCKED;
            INPU = [a_lpha b_eta g_amma d_elta d_A d_B d_C d_D FUCKED FUCKED FUCKED FUCKED];
            %Alpha
            A = subs(A, SYMB, NUMB);
            A = vpa(subs(A, READ, INPU),5)
            %Beta
            B = subs(B, SYMB, NUMB);
            B = vpa(subs(B, READ, INPU),5)
            %Gamma
            C = subs(C, SYMB, NUMB);
            C = vpa(subs(C, READ, INPU),5)
            %Delta
            D = subs(D, SYMB, NUMB);
            D = vpa(subs(D, READ, INPU),5)
            
            obj.H1 = @(a_lpha, b_eta, g_amma, d_elta, d_A, d_B, d_C, d_D)(subs(A));
            obj.H2 = @(a_lpha, b_eta, g_amma, d_elta, d_A, d_B, d_C, d_D)(subs(B));
            obj.H3 = @(a_lpha, b_eta, g_amma, d_elta, d_A, d_B, d_C, d_D)(subs(C));
            obj.H4 = @(a_lpha, b_eta, g_amma, d_elta, d_A, d_B, d_C, d_D)(subs(D));
            
            
        end
            
            
            
        function generateEOMs_N(obj)
            
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
            subs(H(1))
            subs(H(2))
            subs(H(3))
            subs(H(1))
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
               alpha(i) = Y(i,1); %alpha at time step
               beta(i) = Y(i,2); %beta at time step
               gamma(i) = Y(i,3); %gamma at time step
               delta(i) = Y(i,4); %delta at time step
           end       
        end
    end
end
