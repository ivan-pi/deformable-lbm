%% Compile DLBM

mex -R2018a dlbminit.c -g -Warn all 
mex -R2018a dlbmstep.c -g -Warn all
mex -R2018a dlbmintegrate.c -g -Warn all 

%% Helper functions

ds = 1550; % Density of solids
dw = 998; % Density of water

% Helper functions
nsolid = @(X) ds*dw/(dw + ds*X);
nwater = @(X) nsolid(X)*X;

%% Set parameters

n = 21
L0 = 0.001
X0 = 0.3
Xs = 0.12

D = 1.0e-10

dx = L0/n
dt = 0.1*dx^2/D
csqr = 0.25*dx^2/dt^2

tout = (0:4).*dt

%% Run simulation
tbegin = tic;

    [t,X,us,r] = dlbm(tout,csqr,dt,X0,Xs,D,L0,n);

tfinish = toc(tbegin);

% Print run statistics
total_steps = (tout(end) - tout(1))/dt;
mlups = n * total_steps / (tfinish - tend) / 1000000

% 
%Xavg = trapz(r(:,:) .* X(:,:));
