clc
clear
% Physical Properties
r = 0.002;
R0 = 0.034;
Ri = 0.0285;
H = 0.096;
h = 0.005;
L = 0.048;
mr = 0.045;
p = 7200;

% Initial conditions
a_0 = pi/24;
b_0 = pi/24;
c_0 = pi/24;
d_0 = 0;
dA_0 = -0.1;
dB_0 = 0.1;
dC_0 = pi/2;
dD_0 = 500;

% Simulation timings
FPS = 60;
tEnd = 20;

Sim = Simulation;
Sim.init();
Sim.setPhysicalParams(r, R0, Ri, H, h, L, mr, p);
Sim.setInitialConditions(a_0, b_0, c_0, d_0, dA_0, dB_0, dC_0, dD_0);

% fprintf('Generating EOMs...\n')
% Sim.generateEOMs_N();
% 
% fprintf('Done!\nSimulating...\n')

[a,b,c,d] = Sim.simulateGyroscope(FPS, tEnd);

fprintf('Sim done!\nYou can now run the animation.')
hold on
plot (a);

plot (b);

plot (c);

