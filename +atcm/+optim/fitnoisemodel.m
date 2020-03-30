function [Qp,F] = fitnoisemodel(DCM)

% copy model contents
DM = DCM;

% integration (M.IS) = @gen
DM.M.IS = @gen;

% check c parameter is not fixed!
DM.M.pC.c = [1;1]/8;

% optimse
[Qp,Cp,Eh,F] = spm_nlsi_GN(DM.M,DM.xU,DM.xY);



end

function y = gen(P,M,U)

w     = M.Hz;
y     = {P.c(2) * w.^-exp(P.c(1))};  

end
