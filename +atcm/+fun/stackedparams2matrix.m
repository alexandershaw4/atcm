function [matrix,list] = stackedparams2matrix(params)
% e.g.
% for i = 1:num_DCMs
%    params(i) = atcm.get_posteriors_tcm2024(DCM.Ep);
% end
% [matrix,list] = atcm.fun.stackedparams2matrix(params)
%
% AS

% concatenate 3d arrays 
AMPA = tcm_mat2list(cat(3,params.AMPA),'AMPA');
GABAA = tcm_mat2list(cat(3,params.GABAA),'GABAA');
NMDA = tcm_mat2list(cat(3,params.NMDA),'NMDA');
GABAB = tcm_mat2list(cat(3,params.GABAB),'GABAB');

TC = [params.TC]';
CT = [params.CT]';
scale_NMDA = [params.scale_NMDA]';

J = [params.J];

% assemble array
list = []; matrix = [];

f = fieldnames(AMPA);
list = [list; f];
for i = 1:length(f)
    matrix = [matrix AMPA.(f{i})];
end

f = fieldnames(GABAA);
list = [list; f];
for i = 1:length(f)
    matrix = [matrix GABAA.(f{i})];
end

f = fieldnames(NMDA);
list = [list; f];
for i = 1:length(f)
    matrix = [matrix NMDA.(f{i})];
end

f = fieldnames(GABAB);
list = [list; f];
for i = 1:length(f)
    matrix = [matrix GABAB.(f{i})];
end

list = [list; 'TC'];
matrix = [matrix TC];

list = [list; 'CT'];
matrix = [matrix CT];

list = [list; 'scale_NMDA'];
matrix = [matrix scale_NMDA];


% Now remove zero variance parameters
V = var(matrix);

r = find(V==0);

matrix(:,r) = [];
list(r) = [];

end
