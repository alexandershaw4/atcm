function FCforward(P,DCM)
% Wraps the integrator function and computes functional connectivity between 
% nodes using amplitude envelope correlations for descrete frequency bands
%
%


% Run the integration on the DCM
[y,w,s,~,~,pst] = feval(DCM.M.IS,P,DCM.M,DCM.xU);

% Get the contributing population membrane potential time series
Ji = find(exp(P.J));
mV = squeeze(s{1}(:,Ji,1,:));
W  = exp(P.J(Ji));

% If this was a single-node model and we've squeezed out dim 1!
if ndims(mV)==2
    mV = shiftdim(mV,-1);
end

L  = exp(P.L);
ns = size(mV,1);

% Generate LFP for each region from weighted combination of cells (W) and
% electrode gain (L):
for i = 1:ns
    LFPs(i,:) = W * squeeze(mV(i,:,:)) * L(i);
end

% Generate filters for If, get envelopes
If = [4 12; 13 30; 40 80];

for i = 1:size(If,1)
    l = If(i,1)./(1./DCM.M.dt*.5);
    u = If(i,2)./(1./DCM.M.dt*.5);   
    [b(i,:),a(i,:)] = butter(2,[l u]);
end

for i = 1:ns
    for j = 1:size(If,1)
        ENV(i,j,:) = abs(hilbert(filter(b(j,:),a(j,:),squeeze(LFPs(i,:)))));
    end
end

for j = 1:size(If,1)
    r{j} = corr(squeeze(ENV(:,j,:))');
end