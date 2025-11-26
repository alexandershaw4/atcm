function [P,M,x] = init_tcm_demo()

p = mfilename('fullpath');

[p,f,e] = fileparts(p);

load([p '/tcm_init.mat']);

P = tcm_init.pE;

M = tcm_init.M;

x = M.x;



