function s = cylsimps(r,X)

N = length(X);

s = 0;
for i = 1:2:N-2
  d = r(i+2) - r(i);
  
  h1 = r(i+1) - r(i);
  h2 = r(i+2) - r(i+1);
  
  ta = (d*(2*h1 - h2))/(6*h1);
  tb = d^3/(6*h1*h2);
  tc = -(d*(h1 - 2*h2))/(6*h2);
  
  s = s + ta*X(i)*r(i) + tb*X(i+1)*r(i+1) + tc*X(i+2)*r(i+2);
end

s = 2*s/r(end)^2;