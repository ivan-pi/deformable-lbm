function hsum = cyltrapz(r,X)
% CYLTRAPZ Cylindrical integral using trapezoidal method

N = length(X);

hsum = 0;
for i = 1:N-1

    h = r(i+1) - r(i);

    % Gauss integration rule, k = 1
    eta = (2./3.) * (r(i+1)^3 - r(i)^3)/(r(i+1)^2 - r(i)^2);
    w = (r(i+1)^2 - r(i)^2)/(2*h);

    % Linear interpolation to Gauss knot
    Xeta = X(i) + (eta - r(i))*(X(i+1) - X(i))/h;
    hsum = hsum + h*(w*Xeta);

    % Midpoint rule
    %rc = 0.5*(r(i) + r(i+1));
    %Xc = 0.5*(X(i) + X(i+1));
    %hsum = hsum +  h*(rc*Xc);

    % Two-point trapezoidal rule
    %hsum = hsum + h*0.5*(r(i)*X(i) + r(i+1)*X(i+1));
end