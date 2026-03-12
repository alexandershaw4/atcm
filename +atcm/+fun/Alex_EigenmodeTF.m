function [Y,w,modes,units,MAG,PHA] = Alex_EigenmodeTF(P,M,U)
% Modal / eigenmode decomposition of DCM transfer function
% Similar inputs to Alex_LaplaceTFwDNew
%
% Outputs:
%   Y      : cell containing CSD / spectral response
%   w      : frequencies
%   modes  : struct with eigenvalues, eigenvectors, modal contributions
%   units  : linearisation details
%   MAG    : per-state modal magnitudes
%   PHA    : per-state modal phases

if isnumeric(P), P = spm_unvec(P,M.P); end
if isstruct(P) && isfield(P,'p'), P = P.p; end

if isfield(M,'endogenous') && M.endogenous
    Input = 0;
else
    Input = 1;
end

if isfield(M,'fixedpoint') && M.fixedpoint == 1
    x = atcm.fun.alexfixed(P,M,1e-10,[],[],1000);
    M.x = spm_unvec(x,M.x);
end

w   = M.Hz(:);
x0  = M.x(:);
Ns  = size(M.x,1);

[f0,A,D] = feval(M.f,M.x,0,P,M);
A       = full(denan(A));
Bfull   = spm_diff(M.f,M.x,1,P,M,2);
Bfull   = full(denan(Bfull));

Uomega = ones(numel(w),1);
if isfield(M,'external_spectrum')
    Uomega = M.external_spectrum(:);
end

PSD = zeros(Ns,numel(w));
MAG = cell(Ns,1);
PHA = cell(Ns,1);

modes = struct();
modes.region = cell(Ns,1);

for ii = 1:Ns

    win = ii:Ns:length(A);
    n   = numel(win);

    AA  = A(win,win);
    BB  = Bfull(win,:);
    X0  = x0(win);
    Cw  = exp(P.J(win));

    if size(BB,2) > 1
        BB = sum(BB,2);
    end

    drive_scale = 1;
    if isfield(P,'C') && numel(P.C) >= ii
        drive_scale = exp(P.C(ii));
    end

    % Delay-free eigendecomposition first
    [V,L] = eig(AA);
    lam   = diag(L);
    Vinv  = pinv(V);

    modal_resp = complex(zeros(n,numel(w)));
    y          = complex(zeros(1,numel(w)));

    % Modal residues for exogenous case: rk = C*v_k * (w_k^T B)
    % where w_k^T are rows of Vinv
    residues_B = zeros(n,1);
    residues_X0 = zeros(n,1);

    for k = 1:n
        vk = V(:,k);
        wk = Vinv(k,:);
        residues_B(k)  = (Cw(:).' * vk) * (wk * BB);
        residues_X0(k) = (Cw(:).' * vk) * (wk * X0);
    end

    mode_contrib = complex(zeros(n,numel(w)));

    for j = 1:numel(w)
        s = 1i*2*pi*w(j);

        Hj = 0;
        for k = 1:n
            % simple delay-free modal contribution
            if Input
                uk = Uomega(j) * drive_scale;
                mode_contrib(k,j) = residues_B(k) * uk ./ (s - lam(k)) ...
                                  + residues_X0(k)      ./ (s - lam(k));
            else
                mode_contrib(k,j) = residues_X0(k) ./ (s - lam(k));
            end
            Hj = Hj + mode_contrib(k,j);
        end

        y(j) = Hj;

        % state-level modal response
        if Input
            modal_resp(:,j) = V * ((1 ./ (s - lam)) .* (Vinv * (BB * Uomega(j)*drive_scale + X0)));
        else
            modal_resp(:,j) = V * ((1 ./ (s - lam)) .* (Vinv * X0));
        end
    end

    MAG{ii} = abs(modal_resp);
    PHA{ii} = angle(modal_resp) * 180/pi;

    Lgain = 1;
    if isfield(P,'L') && numel(P.L) >= ii
        Lgain = exp(P.L(ii));
    end

    Yloc = y(:);

    H = gradient(gradient(Yloc));
    Yloc = Yloc - (exp(P.d(1))*3)*H;

    PSD(ii,:) = Lgain * Yloc(:).';

    modes.region{ii}.V = V;
    modes.region{ii}.lambda = lam;
    modes.region{ii}.Vinv = Vinv;
    modes.region{ii}.mode_contrib = mode_contrib;
    modes.region{ii}.residues_B = residues_B;
    modes.region{ii}.residues_X0 = residues_X0;
end

CSD = zeros(numel(w),Ns,Ns);
for i = 1:Ns
    CSD(:,i,i) = PSD(i,:).';  % auto
    for j = 1:Ns
        if i ~= j
            Lc = 1;
            if isfield(P,'Lc') && numel(P.Lc) >= i
                Lc = exp(P.Lc(i));
            end

            if Input
                CSD(:,i,j) = Lc * (PSD(i,:).' .* conj(PSD(j,:).'));
            else
                CSD(:,i,j) = Lc * (PSD(i,:).' .* (PSD(j,:).')); 
            end
            CSD(:,j,i) = CSD(:,i,j);
        end
    end
end

% Smooth magnitudes (keeps behaviour of |.| then smooth)
dw = mean(diff(w));
if Ns == 1
    CSD = atcm.fun.agauss_smooth(abs(CSD), dw * exp(P.d(3)));
else
    for i = 1:Ns
        for j = 1:Ns
            CSD(:,i,j) = atcm.fun.agauss_smooth(abs(CSD(:,i,j)), dw * exp(P.d(3)));
            CSD(:,j,i) = CSD(:,i,j);
        end
    end
end

Y     = {(CSD)};


% CSD = zeros(numel(w),Ns,Ns);
% for i = 1:Ns
%     CSD(:,i,i) = PSD(i,:).';
%     for j = 1:Ns
%         if i ~= j
%             CSD(:,i,j) = PSD(i,:).' .* conj(PSD(j,:).');
%             CSD(:,j,i) = CSD(:,i,j);
%         end
%     end
% end
% 
% Y = {CSD};

units = [];
units.x0 = x0;
units.dx = f0;
units.A  = A;
units.B  = Bfull;
units.D  = D;
units.freq = w;