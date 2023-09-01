function yp = pde_isothermal(t,y,lag,D,alpha)
    arguments
        t (1,1) double
        y (:,1) double
        lag (:,1) double
        D (1,1) double
        alpha (1,1) double
    end
% PDE_ISOTHERMAL
%
% t     ... time
% y     ... moisture vector and eulerian position vector
% lag   ... solid lagrangian coordinates as column vector
% D     ... diffusivity
% alpha ... volumetric shrinkage coefficient
%


N = length(y)/2;
X = y(1:N);
r = y(N+1:end);

% Preallocate vector of derivative values
yp = zeros(2*N,1);

% Evaluate effective diffusivity
%Dlag = D*((r./lag)./(1 + alpha*X)).^2;
Dlag = D./(1 + alpha*X); % Following Pakowski

%Dlag(1) = D; % correct for limit diffusivity at center

Dc = 0.5*(Dlag(1:end-1) + Dlag(2:end));
lagc = 0.5*(lag(1:end-1) + lag(2:end));

% spatial step
dlag = lag(2) - lag(1);

invdlag2 = 1.0/dlag^2;

% Evaluate central difference
for i = 2:N-1

  yp(i) = invdlag2*(1./lag(i))*(Dc(i)*lagc(i)*(X(i+1) - X(i)) ...
    - Dc(i-1)*lagc(i-1)*(X(i) - X(i-1)));
end  
  
yp(1) = 4*Dlag(1)*(X(2) - X(1))/dlag^2;
yp(N) = 0; % Dirichlet boundary  

% Evaluate radial displacement right-hand side
yp(N+1) = 0;  % Center doesn't move.

for i = 2:N-1
  yp(N+i) = alpha*Dlag(i)*lag(i)*(X(i+1) - X(i-1))/(2*dlag)/r(i);
end

% One-sided difference on cylinder surface
% TODO: Replace with second order

%yp(end) = alpha*Dlag(end)*lag(end)*(X(end)-X(end-1))/dlag;

% The minus is because of the surface normal
dXdlag = (0.5*X(end-2) - 2.0*X(end-1) + 1.5*X(end))/dlag;
yp(end) = alpha * Dlag(end) * lag(end) * dXdlag / r(end);

%yp(end) = alpha * yp(N-1); % We should use the flux somehow here...
