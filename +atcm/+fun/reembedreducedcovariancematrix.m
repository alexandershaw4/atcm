function cp = reembedreducedcovariancematrix(DCM,CP)

V  = spm_svd(diag(spm_vec(DCM.M.pC)));
cp = V*CP*V';
