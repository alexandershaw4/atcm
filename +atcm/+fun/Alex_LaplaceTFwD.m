function [Y,w,G,units,MAG,PHA] = Alex_LaplaceTFwD(P,M,U)
% linearisation and numerical (Laplace) transfer function for a DCM
% witten by Alex Shaw; this version has proper (Laplace domain) handling of
% delays...
%
% takes a dynamical model and observation function;
%
%   dx = f(x,u,P,M)
%    y = g(x,P)
%
% and linearises flow of states [A], inputs [B] and observation [C];
%
%   dx = Ax  + Bu;
%    y = Cx [+ Du];
%
% having computed the linearisation matrices, computes the frequency response 
% at frequencies of interest (vector stored in M.Hz) using the Laplace transform.
%
%   Y(s) = C*inv(sI - A)*BU + C*inv(sI - A)*x0
%
% Usage: [Y,w,G,units] = atcm.fun.Alex_LaplaceTFwD(P,M,U); 
%
% add sub-structure M.sim with M.sim.pst and M.sim.dt to force the routine
% to also recronstruct the time-domain simulation and return in the 4th
% output, since we can access the magnitude and phase data of each component
% from the Laplace transform.
%
% * Update Feb 2024: AS refactored for multiple node models to compute the
% Laplace transform of each region, then compute the cross-spectral density
%
% AS2023

if isnumeric(P)
    P = spm_unvec(P,M.P);
end

if isstruct(P) && isfield(P,'p')
    P = P.p;
end

if isfield(M,'endogenous') && M.endogenous
    Input = 0;
else
    Input = 1;
end

% fixed point search?
if isfield(M,'fixedpoint') && M.fixedpoint == 1
    x = atcm.fun.alexfixed(P,M,1e-10,[],[],1000);
    M.x = spm_unvec(x,M.x);
end

f = @(x,u,varargin) M.f(x,u,P,M);
w = M.Hz;
x0 = M.x(:);
u0 = 1;

[f,A,D] = feval(M.f,M.x,0,P,M);
A = denan(A);
B = spm_diff(M.f,M.x,1,P,M,2);
B = denan(B);
C = exp(P.J);
Ns = size(M.x,1);

% Delay operator
D_exp = arrayfun(@(w) expm(-1i * 2 * pi * w * D), M.Hz, 'UniformOutput', false);

% Loop each node (aka region, source, mode, column ..)
for i = 1:Ns
    win = i:Ns:(length(A));

    AA = A(win,win);
    BB = B(win)*exp(P.C(i));
    
    % if no input to the system these are endogenous fluctations about x0
    if ~Input
        BB = x0(win);
    end

    % we use a static observer model anyway...
    C = exp(P.J(:));
    
    % properly handle delays in the lpalace domain:
    % D(s) = e^−sD = e^−(jω)D
    for j = 1:length(w)
        D_w = D_exp{j};  % Get delay operator at this frequency
        D_w = D_w(win,win);
        %Jm  = D_w * (AA - 1i*2*pi*w(j)*eye(length(AA)));

        CsI = exp(P.d(2))+(1i*2*pi*w(j)*eye(length(AA)));

        Jm  = D_w * (CsI - AA);
        Ym  = (Jm \ BB) + (Jm\x0(win));
        MG(:,j) = Ym;
        Y   = C' * Ym;
        y(j) = Y;
    end

    Y = y;

    % the matlab way (not in use now)
    %G = ss(AA, BB, diag(C), 0);  % Assuming unity output matrix
    G = [];
    Y = (Y);

    if isfield(M,'ham') && M.ham;
        H = hamming(length(w));
        Y = Y(:).*H(:); 
    end

    MAG{i} = (MG);%magnitude;
    PHA{i} = angle(MG)*180/pi;%phase;

    % electrode scaling
    PSD(i,:) = exp(P.L(i))*(Y);

end

% compute cross spectrum from autospectra using complex cobjugate
CSD = zeros(length(w),Ns,Ns);
for i = 1:Ns
    CSD(:,i,i) = PSD(i,:);
    for j = 1:Ns
        if i ~= j
            CSD(:,i,j) = exp(P.Lc(i)) * PSD(i,:) .* conj(PSD(j,:));
            CSD(:,j,i) = CSD(:,i,j);
        end
    end
end

% now that we've computed the CSD from the *complex* data, we can take abs
% and smooth
for i = 1:Ns
    for j = 1:Ns
        %if i < j
            CSD(:,i,j) = atcm.fun.agauss_smooth(abs(CSD(:,i,j)),mean(diff(w))*exp(P.d(1)));
            CSD(:,j,i) = CSD(:,i,j);
        %end
    end
end


Y = {(CSD)};
units = [];

% if continuous-time simluation was requested, compute series
if isfield(M,'sim') && nargout > 3
    pst = M.sim.pst;
    dt  = M.sim.dt;
    
    % remove sim struct and recall top func
    M = rmfield(M,'sim');
    P.J(P.J==-1000)=0;
    [~,~,~,~,MAG,PHA] = atcm.fun.Alex_LaplaceTFwD(P,M,U);

    for k = 1:Ns
        
        mag = squeeze(MAG{k});
        the = squeeze(PHA{k});
    
        for i = 1:size(mag,1)
            for j = 1:size(mag,2)    
                series{k}(i,j,:) = mag(i,j) * sin(2*pi*w(j)*pst/1000 - the(i,j) );
            end
        end
    
        S{k} = squeeze(sum(series{k},2));
        S{k} = S{k} + spm_vec(M.x(k,:,:));
        LFP(k,:) = (C)'*S{k};
    end

    units.series = S;
    units.LFP    = LFP;
    units.dt     = dt;
    units.pst    = pst;
    units.A      = A;
    units.B      = B;
    units.C      = C;
    units.D      = [];
    units.xinit  = M.x(:);
    
    units.mag    = squeeze(mag);
    units.phase  = squeeze(the);
    units.freq   = w(:);

end

% --- New Impulse Response Calculation ---
if isfield(M,'impulse') && M.impulse && nargout > 3
    t = pst;
    impulse_response = zeros(Ns, length(t));

    for k = 1:Ns
        win = k:Ns:(length(A));
        AA = A(win,win);
        BB = B(win);
        CC = C;
        for idx = 1:length(t)
            impulse_response(k,idx) = CC' * expm(AA * t(idx)) * BB;
        end
    end

    units.impulse_response = impulse_response;
    units.impulse_time = t;
end


return;



