function [Qp,F] = HierOptim(DCM)
%
%
%
%
%

% Split (variance) parameter set into a series of extrinsic and intrinsics
% for different optimsations
pC = DCM.M.pC;
V  = spm_vec(pC);

% Extrinsic parameters
Ex = [ spm_fieldindices(pC,'A')  ; 
       spm_fieldindices(pC,'AN') ; 
       spm_fieldindices(pC,'B')  ; 
       spm_fieldindices(pC,'BN') ; 
       spm_fieldindices(pC,'C')  ; 
       spm_fieldindices(pC,'L')  ] ;

ExV     = V*0;
ExV(Ex) = V(Ex);
Expc    = spm_unvec(ExV,pC);

% Intrinsic parameters
In      = find( full(V) - full(ExV) );
InV     = V*0;
InV(In) = V(In);
Inpc    = spm_unvec(InV,pC);

% calculate number in each
nE = length(find(spm_vec(Expc)));
nI = length(find(spm_vec(Inpc)));

fprintf('Model contains %d extrinsic and %d intrinsic variable params\n',nE,nI);

% Fit part 1:
M    = DCM.M;
M.pC = Expc;
[Qp,Qg,Cp,Cg,Ce,F,LE] = spm_nlsi_N(M,DCM.xU,DCM.xY);  % EM: inversion

% Fit part 2
M.pE = Qp;
M.pC = Inpc;
[Qp,Qg,Cp,Cg,Ce,F,LE] = spm_nlsi_N(M,DCM.xU,DCM.xY);  % EM: inversion


