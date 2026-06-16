function RunTCM_NGNMM_VL_2026(i)
% RunTCM_NGNMM_VL_2026
%
% Top-level DCM/VL fitting script for the 8-population thalamo-cortical
% next-generation neural mass model based on QIF/theta-neuron mean-field
% equations.

% -------------------------------------------------------------------------
% Data & Design
% -------------------------------------------------------------------------
Data.Datasets     = 'NewSZ.txt';
Data.Design.X     = [];
Data.Design.name  = {'undefined'};
Data.Design.tCode = [1];
Data.Design.Ic    = [1];
Data.Design.Sname = {''};
Data.Prefix       = 'aNGNMM8_';
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);

% -------------------------------------------------------------------------
% Model space
% -------------------------------------------------------------------------
T = 0;
F = (T==1);
B = (T==2);
C = [1]';
L = sparse(1,1);

% fix the iddue with finding the initial states .mat file
%--------------------------------------------------------------------------
p = mfilename('fullpath');
p = fileparts(p);
p = strrep(p,'main_scripts','');
addpath(p);

% -------------------------------------------------------------------------
% Loop over datasets
% -------------------------------------------------------------------------
for i = i

    DCM = [];

    [~, fn, fe] = fileparts(Data.Datasets{i});
    DCM.name = [Data.Prefix fn fe];

    DCM.xY.Dfile = Data.Datasets{i};
    Ns           = length(F);

    DCM.xU.X     = Data.Design.X;
    DCM.xU.name  = Data.Design.name;
    tCode        = Data.Design.tCode;
    DCM.xY.Ic    = Data.Design.Ic;
    DCM.Sname    = Data.Design.Sname;

    if exist(DCM.name,'file')
        fprintf('Skipping model %d/%d - already exists!\n( %s )\n', ...
            i,length(Data.Datasets),DCM.name);
        continue;
    end

    % ---------------------------------------------------------------------
    % Extrinsic connectivity
    % ---------------------------------------------------------------------
    DCM.A{1} = F;
    DCM.A{2} = B;
    DCM.A{3} = L;
    DCM.B{1} = DCM.A{1} | DCM.A{2};

    if ~isempty(DCM.xU.X)
        DCM.B(2:length(DCM.xU.X)) = DCM.B;
    end

    DCM.C = C;

    % ---------------------------------------------------------------------
    % New 8-population next-generation thalamo-cortical model
    % ---------------------------------------------------------------------
    [M_ng, pE_ng, pC_ng, x0_ng] = atcm.make_tc_ngnmm_model(Ns);

    DCM.M.Ns      = M_ng.Ns;
    DCM.M.Npop    = M_ng.Npop;
    DCM.M.Nstates = M_ng.Nstates;
    DCM.M.pop     = M_ng.pop;
    DCM.M.popType = M_ng.popType;
    DCM.M.Kmask   = M_ng.Kmask;
    DCM.M.Erev    = M_ng.Erev;
    DCM.M.alpha0  = M_ng.alpha0;

    % DCM-compatible state space: Ns x Npop x Nstate
    DCM.M.x = x0_ng;

    DCM.M.f  = @atcm.tc_ngnmm_fx;
    DCM.M.IS = @atcm.fun.Alex_LaplaceTFwDNew;
    DCM.options.SpecFun = @atcm.fun.Afft;

    % ---------------------------------------------------------------------
    % Frequency range and data preparation
    % ---------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n', i, length(Data.Datasets));

    fq = [1 90];

    DCM.M.U            = sparse(diag(ones(Ns,1)));
    DCM.options.trials = tCode;
    DCM.options.Tdcm   = [1 1000];
    DCM.options.Fdcm   = fq;
    DCM.options.D      = 1;
    DCM.options.han    = 1;
    DCM.options.h      = 4;
    DCM.options.DoData = 1;
    DCM.options.Fltdcm = fq;
    DCM.options.UseButterband = fq;

    DCM.options.analysis      = 'CSD';
    DCM.xY.modality           = 'LFP';
    DCM.options.spatial       = 'LFP';
    DCM.options.model         = 'tc_ngnmm8';
    DCM.options.Nmodes        = length(DCM.M.U);

    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 0;
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 1/2;

    DCM.xY.name = DCM.Sname;

    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1;

    DCM.xY.y = {abs(DCM.xY.y{:})};
    DCM.xY.y{:} = atcm.fun.agauss_smooth(abs(DCM.xY.y{:}),1)';;

    %DCM.xY.y{:} = atcm.fun.agauss_smooth(abs(DCM.xY.y{:}),1)';

    % ---------------------------------------------------------------------
    % Priors and initial state
    % ---------------------------------------------------------------------
    DCM.M.x  = x0_ng;
    DCM.M.pE = pE_ng;
    DCM.M.pC = pC_ng;

    DCM.M.y  = DCM.xY.y;
    DCM.M.Hz = DCM.xY.Hz;

    assert(isequal(size(DCM.M.x), [Ns DCM.M.Npop DCM.M.Nstates]), ...
        'DCM.M.x must be Ns x Npop x Nstate.');
    assert(numel(DCM.M.pE.J) == numel(DCM.M.x(:)), ...
        'P.J must have one element per state in DCM.M.x(:).');
    assert(numel(DCM.M.pC.J) == numel(DCM.M.x(:)), ...
        'pC.J must have one element per state in DCM.M.x(:).');

    ppE = DCM.M.pE;
    ppC = DCM.M.pC;

    % ---------------------------------------------------------------------
    % Model check before fitting
    % ---------------------------------------------------------------------
    fprintf('--------------- MODEL CHECK ---------------\n');

    [f0,J0,D0] = feval(DCM.M.f, DCM.M.x, 0, DCM.M.pE, DCM.M);
    fprintf('Initial derivative norm: %.6f\n', norm(f0));
    fprintf('State shape: [%d %d %d], vector length: %d\n', ...
        size(DCM.M.x,1), size(DCM.M.x,2), size(DCM.M.x,3), numel(DCM.M.x));

    if any(~isfinite(f0)) || any(~isfinite(J0(:))) || any(~isfinite(D0(:)))
        error('NGNMM model returned non-finite f/J/D at initial state.');
    end


    % ---------------------------------------------------------------------
    % Pre-optimise prior modes before VL
    % ---------------------------------------------------------------------
    % targetBands = [
    %     2   5
    %     6  10
    %     10  14
    %     15  25
    %     30  45
    %     45  65
    %     ];
    % 
    % modeOpts = struct();
    % modeOpts.maxIter = 120;
    % modeOpts.nStarts = 6;
    % modeOpts.display = true;
    % 
    % [pE_mode, modeInfo, modeHist] = atcm.preopt_ngnmm_modes(DCM.M.pE, DCM.M, targetBands, modeOpts);
    % 
    % DCM.M.pE = pE_mode;
    % DCM.M.x  = x0_ng;   % reset before fixed-point search under new priors

    % ---------------------------------------------------------------------
    % Optional fixed-point / operating-point search
    % ---------------------------------------------------------------------
    fprintf('--------------- STATE ESTIMATION ---------------\n');

    x_template = DCM.M.x;

    f_init = feval(DCM.M.f, x_template, 0, DCM.M.pE, DCM.M);
    n_init = norm(f_init(:));

    try
        x_fp = atcm.fun.find_fixed_point_robust(DCM.M.pE, DCM.M, 1e-10);

        % find_fixed_point_robust may return a vector; reshape into DCM.M.x form
        x_fp = spm_unvec(x_fp(:), x_template);

        f_fp = feval(DCM.M.f, x_fp, 0, DCM.M.pE, DCM.M);
        n_fp = norm(f_fp(:));

        fprintf('Initial residual: %.6e\n', n_init);
        fprintf('Candidate operating-point residual: %.6e\n', n_fp);

        if isfinite(n_fp) && n_fp < n_init
            DCM.M.x = x_fp;
            fprintf('Using improved operating point, although residual is not exact.\n');
        else
            DCM.M.x = x_template;
            fprintf('Keeping original x0.\n');
        end

    catch ME
        warning('Fixed point search failed. Continuing with x0.');
        warning('%s', ME.message);
        DCM.M.x = x_template;
    end

    f_check = feval(DCM.M.f, DCM.M.x, 0, DCM.M.pE, DCM.M);
    fprintf('Operating-point derivative norm used for TF: %.6e\n', norm(f_check(:)));

    % Check oscillatory modes
    [f,A,D] = feval(DCM.M.f, DCM.M.x, 0, DCM.M.pE, DCM.M);

    ev = eig(A);
    freq = abs(imag(ev)) ./ (2*pi);
    damp = real(ev);

    T = table(real(ev), imag(ev), freq, 'VariableNames', {'Real','Imag','Hz'});
    disp(sortrows(T,'Hz'));

    % show resonant modes rather than just oscillatory
    ev = eig(A);

    freq = abs(imag(ev)) ./ (2*pi);
    damp = real(ev);
    Q = abs(imag(ev)) ./ max(eps, -2*damp);

    T = table(real(ev), imag(ev), freq, Q, ...
        'VariableNames', {'Real','Imag','Hz','Q'});

    T = sortrows(T,'Hz');
    disp(T);

    % fprintf('\nPotential spectral modes:\n');
    % ix = find(freq > 1 & freq < 90 & real(ev) < 0);
    % [~,ord] = sort(Q(ix),'descend');
    % disp(T(ix(ord),:));

    % ---------------------------------------------------------------------
    % Parameter estimation
    % ---------------------------------------------------------------------
    fprintf('--------------- PARAM ESTIMATION ---------------\n');

    M = aFitDCM(DCM);

    M.aloglikVLthermGC;
    M.update_parameters(M.Ep);

    [y,w,G,s] = feval(DCM.M.IS, spm_unvec(M.Ep,DCM.M.pE), DCM.M, DCM.xU);

    numit = 0;
    while cdist(DCM.xY.y{:}(:)', y{:}(:)') > 0.5 && numit < 8
        numit = numit + 1;

        M.aloglikVLthermGC;
        M.update_parameters(M.Ep);

        [y,w,G,s] = feval(DCM.M.IS, spm_unvec(M.Ep,DCM.M.pE), DCM.M, DCM.xU);
    end

    Qp   = spm_unvec(M.Ep, DCM.M.pE);
    Cp   = M.CP;
    Fval = M.F;

    % ---------------------------------------------------------------------
    % Save outputs
    % ---------------------------------------------------------------------
    DCM.M.pE = ppE;
    DCM.M.pC = ppC;

    DCM.Ep = Qp;
    DCM.Cp = Cp;
    DCM.F  = Fval;

    DCM.M.sim.dt  = 1/600;
    DCM.M.sim.pst = 1000*((0:DCM.M.sim.dt:(2)-DCM.M.sim.dt)');

    [y,w,G,s] = feval(DCM.M.IS, DCM.Ep, DCM.M, DCM.xU);

    DCM.pred   = y;
    DCM.w      = w;
    DCM.G      = G;
    DCM.series = s;

    save(DCM.name, 'DCM');

    close all;
    clear global;

end

end
