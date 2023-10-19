function m = c_oof(w,y,model,varargin)
% 2n order polyfit, or other, with constraint, where constraint is
%
% m    = p(w,y) s.t
% m(1) = y(1) & m(end) = y(end)
%
% AS

if nargin < 3 || isempty(model)
    model = 'poly2';
end

w = w(:);
y = real(y(:));
c = [y(1) y(end)];

warning off;
if ~isempty(varargin)
    p = fit(w,y,model,varargin{:});
else
    p = fit(w,y,model);
end
warning on;

m = p(w);

% constraints
C  = [c(1)./m(1) c(2)./m(end)];
lm = linspace(C(1),C(2),length(m));
m  = m(:).*lm(:);

while any(m>y)
   m  = m + max(y-m);
    m  = [m(2:end); m(end)];
    dy = spm_vec( linspace(c(1),c(2),length(w)) );
    m  = rescale(m) .* dy;
end


end