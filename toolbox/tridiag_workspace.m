function ws = tridiag_workspace(nel,h,D)
%
% Currently only constant diffusivity D, and 
% a constant mesh spacing h are supported

ws.w = zeros(N,4);
ws.ipiv = zeros(N,'int32');

oc_fill_and_factor(s.w,s.ipiv);

end function

%arguments (Output)
%    w (nel*3+1,4) double
%    ipiv (nel*3+1,1) {mustBeNonempty,mustBeUnderlyingType(ipiv,"int32")}
%end