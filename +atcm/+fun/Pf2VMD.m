function model = Pf2VMD(y,n)
% returns sqrt(y) and n central frequencies of the n IMFs discovered by
% variational mode decomposition
%
% model = Pf2VMD(y,n)
%
% AS

y = denan(y);

[~,~,data] = atcm.fun.vmd(real(y),'NumIMFs',n);

model = {sqrt(y) data.CentralFrequencies};