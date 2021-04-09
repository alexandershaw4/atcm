function con = contributions_dpdy(DCM,Ep)
% Perform a contribution (aka sensitivity) analysis: compute the relative
% effect of each parameter (p) on the model output spectrum (y) as percentage-change 
% from baseline, via finite difference methods - i.e. dp/dy. These are 
% effectively scaled partial numerical derivatives.
%
% AS

% Generate the function
f = @(p) spm_vec(feval(DCM.M.IS,spm_unvec(p,DCM.M.pE),DCM.M,DCM.xU));

% Get x0 / y0
if nargin <2 || isempty(Ep)
    Ep = DCM.M.pE;
end
Ep = spm_vec(Ep);
y0 = f(Ep);

% loop parameters
p = Ep;
v = ~~spm_vec(DCM.M.pC);
c = spm_vec(DCM.M.pC);
n = 0;
N = atcm.fun.generate_pnames(DCM.M.pE);

redmap = zeros(length(p),length(find(v)));
pspec = zeros(length(p),length(y0));

for i = 1:length(p)
    if i > 1; fprintf(repmat('\b',[1 length(txt)])); end
    txt = sprintf('Analysing parameters: %d/%d',i,length(p));
    fprintf(txt);
    if v(i)
        n      = n + 1;
        px     = p;
        px(i)  = px(i) + px(i) * c(i); 
        m(:,i) = f(px);
        
        pspec(i,:) = 100*(m(:,i)-y0)./y0;
        redmap(i,n) = 1;
        
        % also build a library
        par.(N{i}) = pspec(i,:)';
    end
end

con.f = f;
con.p = p;
con.v = v;
con.c = c;
con.m = m;
con.w = DCM.xY.Hz;
con.baseline = y0;
con.RelChangeSpec = pspec;
con.redmap = redmap;
con.par = par;
con.Names = N;
con.pindices = find(spm_vec(DCM.M.pC));

% find the rank of the covariance of the derivative matrix and project
% principal components:
%
% RelCh = con.redmap'*con.RelChangeSpec; % active params only
% N = rank(cov(RelCh')); % cov rank
% [v,D] = eig(cov(RelCh')); % decompose covariance matrix
% DD = diag(D); [~,ord]=sort(DD,'descend'); % sort eigenvalues
% PCV = v(:,ord(1:N))*D(ord(1:N),ord(1:N))*v(:,ord(1:N))'; % project factorised matrix without rank deficiency
%
% Check that the principal components are orthogonal:
% corr((v(:,ord(1:N))'*RelCh)')