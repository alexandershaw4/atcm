function demo_run_tcm_modular()
ns = 1; Reg = build_default_registries();
Comp = compile_tcm(Reg, ns);
% Minimal M, x init
np = numel(Reg.Pops); nk = 1 + numel(Reg.Channels); 
M.x = zeros(ns,np,nk);
%P = struct('BE',log(0.8), 'w_SP_SS_AMPA',log(1), 'w_SP_SS_NMDA',log(1), 's5HT2A',log(0.1));
% Call state eq once
%x0 = spm_vec(M.x);
M.Comp = Comp;
M.Reg = Reg;
M.f = @tcm_modular;

% get parameters
opts.free_mode = 'class-tied';
Pri = setup_priors_extended(opts);

M.labels = Pri.labels;
P = Pri.vec_to_struct(Pri.pE)
C = Pri.vec_to_struct(Pri.pC)
C = spm_unvec(~~spm_vec(P)/8,P);

x0 = spm_vec(zeros(1, numel(Reg.Pops), 1+numel(Reg.Channels)));

[f,~,~] = tcm_modular(x0, [0;0], P, M); 


DCM.M = M;
DCM.M.pE = P;
DCM.M.pC = C;

% fixed point
x = atcm.fun.alexfixed(DCM.M.pE,DCM.M,1e-10,[],[],[],[],1);
DCM.M.x = spm_unvec(x,DCM.M.x);

% transfer function
DCM.M.IS = @atcm.fun.Alex_LaplaceTFwD;
DCM.M.Hz = 1:90;

DCM.M.pE.J = DCM.M.x*0;DCM.M.pE.J(:,:,1) = log([.6 .8 .4 .6 .4 .6 .4 .4]);
DCM.M.pE.J = DCM.M.pE.J(:);
DCM.M.pC.J = DCM.M.pE.J(:)*0;

DCM.M.pE.L = 0;
DCM.M.pC.L = 1/8;

DCM.M.pE.C = 1;
DCM.M.pC.C = 1/8;

DCM.M.pE.d = [0 0];
DCM.M.pC.d = [1 1]./8;

[y,w,G,s] = feval(DCM.M.IS,DCM.M.pE,DCM.M,[]);

end
