function A = zero2log(A)
% Convert a matrix of parameters including zeros to a logged matrix,
% without log(0) returning inf / nan.
%
%
%

v = spm_vec(A);
i = find(v==0);

v    = log(v); 
v(i) = -1000;

A    = spm_unvec(v,A);

