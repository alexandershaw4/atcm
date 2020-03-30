function [ymod,Yreal,pst] = evalERP(DCM,pE,gE)
% quickly integrate and observe an erp model of the form:
%
% { y = g(x,...)
% { x = f(x,...)
%
% with dct residuals, as per
% 
%  y{i} = R*f(x,...)*L'
%
% where L = g(x)
%
% AS

if nargin <3
    gE = DCM.M.gE;
end
if nargin < 2
    pE = DCM.M.pE;
end

% f(x)
%--------------------------------------------------------------
[x,pst,fire] = feval(DCM.M.IS,pE,DCM.M,DCM.xU);

% g(x)
%--------------------------------------------------------------
L   = feval(DCM.M.G,gE,DCM.M);

% R
%--------------------------------------------------------------
R  = DCM.xY.R;
M  = DCM.M;
xY = DCM.xY;

% implement R*y{i}*L'
%--------------------------------------------------------------
% for i = 1:length(y)
%     y{i} = R*y{i}*L';
% end

Ns  = length(DCM.A{1}); 
Nr  = size(DCM.C,1); 
Nt  = length(x);
j   = find(kron( exp(gE.J),               ones(1,Nr)));      % Indices of contributing states
%j    = find(kron( logical(sum(Qg.J)), ones(1,Nr)));

NNs = size(xY.y{1},1);
x0  = ones(NNs,1)*spm_vec(M.x)';         % expansion point for states
for i = 1:Nt
    K{i} = x{i} - x0;                   % centre on expansion point
    y{i} = M.R*K{i}*L'*M.U;             % prediction (sensor space)
    r{i} = M.R*xY.y{i}*M.U - y{i};      % residuals  (sensor space)
    K{i} = K{i}(:,j);                   % Depolarization in sources
end

ymod  = y;
Yreal = spm_unvec( spm_vec(y)+spm_vec(r), y);
pst = xY.pst;
