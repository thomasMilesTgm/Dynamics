clc

% Physical Properties
r = 0.002;
R0 = 0.034;
Ri = 0.0285;
H = 0.096;
h = 0.005;
L = 0.048;
mr = 0.045;
p = 4.4;

% Initial conditions
a_0 = 0;
b_0 = 0.0001;
c_0 = 0.0001;
d_0 = 0;
dA_0 = 0;
dB_0 = 0;
dC_0 = 0;
dD_0 = 500;

% Simulation timings
FPS = 10;
tEnd = 10;

Sim = Simulation;
Sim.init();
Sim.setPhysicalParams(r, R0, Ri, H, h, L, mr, p);
Sim.setInitialConditions(a_0, b_0, c_0, d_0, dA_0, dB_0, dC_0, dD_0);

fprintf('Generating EOMs...\n')
Sim.generateEOMs();

fprintf('Done!\nSimulating...\n')

[a,b,c,d] = Sim.simulateGyroscope(FPS, tEnd);

fprintf('Sim done!\nYou can now run the animation.')