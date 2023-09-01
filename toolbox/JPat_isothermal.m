function out = JPat_isothermal(M)
% JPAT
%
% N ... the numer of points
%
% The output matrix will be of size (2*N,2*N)

%
% S11, Gradient of dX/dt w.r.t to X
%

S11 = spdiags(ones(M,3),-1:1,M,M); % Finite-difference in space

%
% S12, Gradient of dX/dt w.r.t to r
%

S12 = sparse(M,M);

%
% S21, Gradient of dr/dt w.r.t to X
%

block = zeros(M,4);

% Central difference
block(:,2) = 1;
block(:,4) = 1;

% One-sided difference
block(end-2,1) = 1;
block(end-1,2) = 1;
block(end  ,3) = 1;

% Center, r = 0, no movement
block(2,4) = 0;

S21 = spdiags(block,-2:1,M,M);
%S21 = spdiags(ones(M,4),-2:1,M,M);

%
% S22, Gradient of dr/dt w.r.t to r
%
%S22 = speye(M,M);
S22 = spdiags([0; ones(M-1,1)], 0, M, M);


%
% --- Assemble blocks ---
%
out = [S11 S12; S21 S22];
