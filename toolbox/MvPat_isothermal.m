function out = MvPat_isothermal(N)
S1 = sparse(N,N);
diag = [0; ones(N-1,1)];
S2 = spdiags(diag,0,N,N);
out = [S1 S1; S1 S2];
end