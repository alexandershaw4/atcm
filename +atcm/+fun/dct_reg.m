function yq = dct_reg(x, y, xq, N)    
% N-dim polynomial regression/spline
%
% AS

%X = apolybasis(x,N);
X = spm_dctmtx(length(x),N);

coff = X\y(:);

yq = spm_dctmtx(length(xq),N)*coff;

end