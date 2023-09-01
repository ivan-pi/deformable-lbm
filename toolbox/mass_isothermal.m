function M = mass_isothermal(~,y)

N = length(y)/2;

M1 = speye(N);
M2 = sparse(N,N);
M3 = sparse(N,N);

diag = [1; y(N+2:end)];
M4 = spdiags(diag,0,N,N);

M = [M1 M2; M3 M4];