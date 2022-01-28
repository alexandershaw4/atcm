function model = Pf2GMM(y)
% returns sqrt(y) and n central frequencies of the n IMFs discovered by
% variational mode decomposition
%
% model = Pf2VMD(y,n)
%
% AS

y = denan(y);

w = 1:length(y);

try
    m = fit(w.',real(y),'Gauss4');
    sv = spm_vec(coeffvalues(m));
    model = {sqrt(y) spm_vec(sort(sv,'descend'))};
catch
    model = {sqrt(y) zeros(12,1)};
end