function [EP,BETA,R,J,COVB,MSE] = nonlinear_ls_obj(P,DCM)
% System identifcation by nonlinear least squares
%
% [EP,BETA,R,J,COVB,MSE] = nonlinear_ls_obj(P,DCM)
%
% AS

global DD

DD    = DCM;
DD.SP = P;
P     = spm_vec(P);
V     = spm_vec(DCM.M.pC);
ip    = find(V);
cm    = zeros(length(V),length(ip));

% make and store a mapping matrix
for i = 1:length(ip)
    cm(ip(i),i) = 1;
end

% to pass to f(ßx)
DD.P  = P;
DD.V  = V;
DD.cm = cm;

% M  = DCM.M;
% U  = DCM.xU;
% 
% IS = @DCM.M.IS;                               % wraps integrator
% f  = @(x) spm_vec(IS(spm_unvec(x,P),M,U));    % wraps above for 1 vector input call
% e  = @(x) sum( spm_vec(DCM.xY.y) - f(x) ).^2; % objective function


fprintf('Performing non-linear least-squares optimisation\n');

p = ones(length(ip),1);
b = p;

options = statset;
options.Display = 'iter';
options.TolFun  = 1e-6;
options.MaxIter = 128;
options.FunValCheck = 'off';

[BETA,R,J,COVB,MSE] = atcm.optim.nlinfit(p,spm_vec(DCM.xY.y),@fakeDM,b,options);


%[X,F,CP]  = AOf(@fakeDM,p(:),c,DCM.xY.y,niter,12*4,[],-inf);
[~,EP]    = fakeDM(BETA,p);


end

function [y,PP] = fakeDM(B,Px,varargin)
global DD

P    = DD.P;
cm   = DD.cm;

X0 = sum(diag(B.*Px)*cm');
X0(X0==0) = 1;
X0 = full(X0.*exp(P'));
X0 = log(X0);
X0(isinf(X0)) = -1000;

PP = spm_unvec(X0,DD.SP);

IS   = spm_funcheck(DD.M.IS);       % Integrator
y    = IS(PP,DD.M,DD.xU);           % Prediction
y    = spm_vec(y);
y    = real(y);


end