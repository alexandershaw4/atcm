function [p,X,M] = batch_analyse_params(list)

X  = atcm.fun.ReadDatasets(list);

for i = 1:length(X)
    load(X{i});
    M{i}  = DCM;
    p(i)  = atcm.fun.unpack_parameters(DCM,DCM.Ep);
end