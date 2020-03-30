function DCM = prepdata(DCM,Ns,tCode)
% Sets options for the empirical data, with new options including
% baselining, filtering and using a set % of trials.

    DCM.M.U            = sparse(diag(ones(Ns,1)));  %... ignore [modes]
    DCM.options.trials = tCode;                     %... trial code [GroupDataLocs]
    DCM.options.Tdcm   = [1 500];                   %... peristimulus time
    DCM.options.Fdcm   = [4 80];                    %... frequency window
    DCM.options.D      = 1;                         %... downsample
    DCM.options.han    = 1;                         %... apply hanning window
    DCM.options.h      = 4;                         %... number of confounds (DCT)
    DCM.options.DoData = 1;                         %... leave on [custom]
    %DCM.options.Bdcm   = [-200 0];                  %... baseline times [new!]
    DCM.options.Fltdcm = [4 80];                    %... bp filter [new!]

    DCM.options.analysis      = 'CSD';              %... analyse type
    DCM.xY.modality           = 'LFP';              %... ECD or LFP data? [LFP]
    DCM.options.spatial       = 'LFP';              %... spatial model [LFP]
    DCM.options.model         = 'tc6';              %... neural model
    DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes

    % Alex additions - 1010 = use aFFT
    DCM.options.UseWelch = 1010;
    DCM.options.Smooth   = 4;
    DCM.options.UseButterband = [4 80];
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);

    DCM.options.DATA = 1 ;
end