clc
clear
% Physical Properties
r = 0.2;
R0 = 3.4;
Ri = 2.85;
H = 9.6;
h = 0.5;
L = 4.8;
mr = 45;
p = 10;

% Initial conditions
a_0 = 0;
b_0 = 0.005;
c_0 = 0.005;
d_0 = 0;
dA_0 = 0;
dB_0 = 0;
dC_0 = 0;
dD_0 = 300;

% Simulation timings
FPS = 50;
tEnd = 20;

Sim = Simulation;
Sim.init();
Sim.setPhysicalParams(r, R0, Ri, H, h, L, mr, p);
Sim.setInitialConditions(a_0, b_0, c_0, d_0, dA_0, dB_0, dC_0, dD_0);

fprintf('Generating EOMs...\n')
Sim.generateEOMs_G();

fprintf('Done!\nSimulating...\n')

[a,b,c,d] = Sim.simulateGyroscope(FPS, tEnd);

fprintf('Sim done!\nYou can now run the animation.')
figure;
hold on
plot (a);
plot (b);
plot (c);
plot (d);