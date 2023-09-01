function [tout,Xout,usout,rout] = dlbm(tspan,csqr,dt,X0,Xs,D,L0,n)
%
% Inputs:
%   tspan - a row vector, of length >= 2
%   csqr - speed of sound squared
%   dt - time step
%   D - diffusivity
%   L0 - length (or radius of body)
%   X0 - initial concentration (dry basis)
%   Xs - surface concentration (dry basis)
%   n - number of cells in space
%
% Outputs:
%   tout - vector of time
%   Xout - moisture content (dry basis)
%   usout - solid velocity
%   rout - position in slab or cylinder

% Initialize lattice
[X,us,r,N] = dlbminit(int32(n),csqr,dt,X0,L0);

nout = len(tspan);

% Allocate output arrays
tout  = zeros(nout,1);
Xout  = zeros(nout,n);
usout = zeros(nout,n);
rout  = zeros(nout,n);

% Store initial condition
tout(1) = tspan(1);
Xout(1,:)  = X;
usout(1,:) = us;
rout(1,:)  = r;

% Time loop
for k = 2:len(tspan)

    % How many steps should we take?
    nsteps = round( (tspan(k) - tspan(k-1))/dt );

    % Run DLBM
    [X,us,r,N] = dlbmstep(int32(nsteps),D,cs,int32(N),r);

    % Store output
    tout(k) = tout(k-1) + nsteps*dt
    Xout(k,:)  = X; 
    usout(k,:) = us;
    rout(k,:)  = r;
end