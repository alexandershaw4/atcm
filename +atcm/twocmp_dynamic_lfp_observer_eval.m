function C = twocmp_dynamic_lfp_observer_eval(obs,s)
% TWOCMP_DYNAMIC_LFP_OBSERVER_EVAL
% -------------------------------------------------------------------------
% Evaluate a dynamic synaptic-current observer produced by
% build_twocmp_dynamic_lfp_observer_all8.m.
%
%       C(s) = L * [ exp(-s*tauE) C_E - alpha exp(-s*tauI) C_I ]
%
% For spectral/transfer-function code, use s = 1i*2*pi*f.
%
% If s is scalar, C is a matrix.  If s is a vector, C is a cell array with
% one observer matrix per frequency, avoiding accidental huge 3-D arrays.
% -------------------------------------------------------------------------

    if numel(s) > 1
        C = cell(numel(s),1);
        for k = 1:numel(s)
            C{k} = twocmp_dynamic_lfp_observer_eval(obs,s(k));
        end
        return;
    end

    eE = exp(-s*obs.tauE);
    eI = exp(-s*obs.tauI);

    Csrc = eE*obs.CEsrc - obs.alpha*eI*obs.CIsrc;
    C    = obs.L * Csrc;
end
