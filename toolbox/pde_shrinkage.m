% pde_shrinkage.c 
%
% This is the right-hand side for an ODE solver
%
% dydt = pde_shrinkage(t,y,h,zeta,D,alpha)
%
%   where
%
%  t ... time
%  y ... the dependent variables [X; r]
%  h ... the grid spacing
%  zeta ... the grid coordinates, with spacing h
%  D ... diffusivity
%  alpha ... volumetric shrinkage coefficient
%
%  X, r, and  zeta are all of the same size N, which is size(y)/2