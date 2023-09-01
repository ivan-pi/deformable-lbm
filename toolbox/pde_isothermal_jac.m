function dfdy = pde_isothermal_jac(t,y,lag,D,alpha)
    arguments
        t (1,1) double
        y (:,1) double
        lag (:,1) double
        D (1,1) double
        alpha (1,1) double
    end

N = length(y)/2;

X = y(1:N);
r = y(N+1:end);

Dlag = D./(1 + alpha*X);

Dc = 0.5*(Dlag(1:end-1) + Dlag(2:end));
lagc = 0.5*(lag(1:end-1) + lag(2:end));

dlag = lag(2) - lag(1);
invdlag2 = 1.0/dlag^2;

%
% J11, Gradient of dX/dt w.r.t to X
%

dXdX = zeros(N,3);

% Evaluate central difference
  %yp(i) = invdlag2*(1./lag(i))*(Dc(i)*lagc(i)*(X(i+1) - X(i)) ...
   % - Dc(i-1)*lagc(i-1)*(X(i) - X(i-1))); 

for i = 2:N-1
    dXdX(i-1,1) = invdlag2*(1.0/lag(i))*(Dc(i-1)*lagc(i-1)); 
    dXdX(i,2)   = invdlag2*(1.0/lag(i))*(-Dc(i)*lagc(i) - Dc(i-1)*lagc(i-1));
    dXdX(i+1,3) = invdlag2*(1.0/lag(i))*(Dc(i)*lagc(i)); 
end

% yp(1) = 4*Dlag(1)*(X(2) - X(1))/dlag^2;

dXdX(1,2) = -4*Dlag(1)/dlag^2;
dXdX(2,3) =  4*Dlag(1)/dlag^2;


J11 = spdiags(dXdX,-1:1,N,N);

%
% J12, Gradient of dX/dt w.r.t to r
%

J12 = sparse(N,N);

%
% J21, Gradient of dr/dt w.r.t to X
%
% (--!-- Dlag also depends on X --!--)
%
drdx = zeros(N,4);


% yp(N+i) = alpha*Dlag(i)*lag(i)*(X(i+1) - X(i-1))/(2*dlag)/r(i);

drdx(:,2) = -alpha.*Dlag.* lag ./(2*dlag)./r;
%drdx(:,3) =  alpha.*lag(i)*()
drdx(:,4) =  alpha.*Dlag.* lag ./(2*dlag)./r;

drdx(2,4) = 0;

%dXdlag = (0.5*X(end-2) - 2.0*X(end-1) + 1.5*X(end))/dlag;
%yp(end) = alpha * Dlag(end) * lag(end) * dXdlag / r(end);

drdx(end-2,1) = alpha * Dlag(end) * lag(end) / r(end) * (0.5/dlag);
drdx(end-1,2) = alpha * Dlag(end) * lag(end) / r(end) * (-2.0/dlag);
drdx(end,3)   = alpha * Dlag(end) * lag(end) / r(end) * (1.5/dlag);

%drdx()

J21 = spdiags(drdx,-2:1,N,N);


%
% J22, Gradient of dr/dt w.r.t. to r
%

drdt_dr = zeros(N,1);

% Evaluate radial displacement right-hand side
drdt_dr(1) = 0;  % Center doesn't move.

for i = 2:N-1
    a = alpha*Dlag(i)*lag(i)*(X(i+1) - X(i-1))/(2*dlag);
    drdt_dr(i) = -a/r(i)^2;
end

dXdlag = (0.5*X(end-2) - 2.0*X(end-1) + 1.5*X(end))/dlag;

a = alpha * Dlag(end) * lag(end) * dXdlag;
drdt_dr(end) = -a/r(end)^2;

J22 = spdiags(drdt_dr,0,N,N);

%
% Combine blocks
%

dfdy = [J11 J12; J21 J22];