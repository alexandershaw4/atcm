function [T,matrix,list] = param_struct_to_table(params)

% concatenate 3d arrays 
AMPA = atcm.fun.tcm_mat2list(cat(3,params.AMPA),'AMPA');
GABAA = atcm.fun.tcm_mat2list(cat(3,params.GABAA),'GABAA');
NMDA = atcm.fun.tcm_mat2list(cat(3,params.NMDA),'NMDA');
GABAB = atcm.fun.tcm_mat2list(cat(3,params.GABAB),'GABAB');

TC = [params.TC]';
CT = [params.CT]';
scale_NMDA = [params.scale_NMDA]';

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

T = array2table(matrix,'VariableNames',list');
