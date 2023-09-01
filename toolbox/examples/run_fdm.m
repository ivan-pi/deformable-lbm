clear
close all
tic
%%

ds = 1510; % density of solid
dw = 990; % density of water

N = 3001;
radius = 0.001;
X0 = 0.3;
alpha = ds/dw;
D = 1.e-10;

ns0 = dw*ds/(dw + ds*X0);
V0 = pi*radius^2;
dry_V0 = V0/(1 + alpha*X0);
dry_radius = sqrt(dry_V0/pi)

lag = linspace(0,dry_radius,N)';
r0 = linspace(0,radius,N)';

myy0 = [ones(N-1,1)*X0; 0.1];

myr0 = (lag(end)^2 + 2*alpha*cyltrapz(lag,myy0))

%%
y0 = zeros(2*N,1);
y0(1:N) = X0;
y0(N) = 0.1;
y0(N+1:end) = r0;
y0(end) = sqrt(myr0)

tspan = [0 6400];


h = lag(2) - lag(1);
rhs = @(t,y) pde_shrinkage(t,y,h,lag,D,alpha);

%{
opts = odeset('Mass',@mass_isothermal,'MStateDependence','strong',...
  'JPattern',JPat_isothermal(N),'MvPattern',MvPat_isothermal(N), ...
  'RelTol',1e-5,'AbsTol',1e-7,'InitialStep',1);
%}

rhs_jac = @(t,y) pde_isothermal_jac(t,y,lag,D,alpha);

opts = odeset('JPattern', JPat_isothermal(N), ...
    'RelTol',1e-6,'AbsTol',1e-9,'InitialStep',1);

%opts = odeset('Jacobian', rhs_jac, ...
%    'RelTol',1e-6,'AbsTol',1e-9,'InitialStep',1);


sol = ode15s(rhs,tspan,y0,opts);

% plot(lag,sol.y(1:N,:))
plot(lag,sol.y(N+1:end,:),'-')

%%
figure

y = deval(sol,[0,30,60,120,240,480,960,1600,3200,6400]);

plot(y(N+1:end,:),y(1:N,:),'-','LineWidth',1.2)

legend("t = 0","t = 30", "t = 60", "t = 120", "t = 240", "t = 480", ...
    "t = 960", "t = 1600", "t = 3200", "t = 6400")

xlabel("r [m]");
ylabel("X(r), moisture content on a dry basis")
title("Moisture profiles")

%%

xavg = zeros(length(sol.x),1);
for i = 1:length(xavg)
    Rend = sol.y(end,i);
    %xavg(i) = (2./Rend^2)*cyltrapz(sol.y(N+1:end,i),sol.y(1:N,i));
    xavg(i) = (2./dry_radius^2)*cyltrapz(lag,sol.y(1:N,i));
end 

%%
figure
hold on
title("Volumetric shrinkage")

Xfine = linspace(0.1,X0,101);
shrink = (1 + alpha*Xfine)/(1 + alpha*xavg(1));

plot(Xfine,shrink,'-','LineWidth',1.5)
plot(xavg(1:4:end),(sol.y(end,1:4:end)/y0(end)).^2,'o')

xlabel("X_{average}")
ylabel("V/V_0")

%%

figure
hold on
title("Volumetric shrinkage - error")


myvfrac = (sol.y(end,:)'/y0(end)).^2;
vfrac = @(X) (1 + alpha*X)./(1 + alpha*xavg(1));

abs_err = (sol.y(end,1:4:end)'/y0(end)).^2 - vfrac(xavg(1:4:end));
rel_err = abs_err./vfrac(xavg(1:4:end));

plot(xavg(1:4:end),rel_err * 100,'o-')
xlabel("X average")
ylabel("Relative error [%]")

%%
figure
hold on
title("Average moisture content")

yyaxis left
plot(sol.x,xavg,'-','LineWidth',1.5)
xlabel("Time [s]")
ylabel("X [/], moisture dry basis")

yyaxis right
myrfrac = sol.y(end,:)'/y0(end);
rfrac = sqrt(vfrac(xavg));
plot(sol.x,(myrfrac - rfrac)./rfrac*100, 'LineWidth',1.5)
ylabel("Relative error [%] - radius")

%%
str = sprintf('radius, relerr = %.10f, %.10f' , sol.y(end,end), (myrfrac(end) - rfrac(end))./rfrac(end)*100)
disp(str)
toc