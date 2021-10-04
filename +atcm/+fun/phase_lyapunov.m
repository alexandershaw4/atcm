function [LLE,LL] = phase_lyapunov(P,M,U,ix,label,px,doplot)
% Compute the phase space plot and Lyapunov exponent over a discrete
% period using the Rosenstein method.
%
% LLE = phase_lyapunov(P,M,U,ix,label)
%
% P  = Param structure
% M  = Model spec structure
% U  = inputs and trial spec structure
% ix = indices of the 3 state outputs of interest, e.g.
%      y = f(x,[],P,M);
%      y = y(ix,:);
%
% label = labels of the 3 states of interest
%
%
% AS

% lle = atcm.fun.phase_lyapunov(DCM.Ep,DCM.M,DCM.xU,[
%
% FOR THE 56-STATE THALAMO-CORTICAL MODEL IN atcm.tc_dev:
%                                        L4 L2 L2 L5 L5 L6 TH TH
% states are ordered y(1:8)   = mV     {[ss sp si dp di tp rt rl]}
%                    y(9:16)  = gAMPA  {[ss sp si dp di tp rt rl]}
%                    y(17:24) = gGABAA {[ss sp si dp di tp rt rl]}
%                    y(25:32) = gNMDA  {[ss sp si dp di tp rt rl]}
%                    y(33:40) = gGABAB {[ss sp si dp di tp rt rl]}
%                    y(41:48) = gM     {[ss sp si dp di tp rt rl]}
%                    y(49:56) = gH     {[ss sp si dp di tp rt rl]}

%                                        1  2  3  4  5  6  7  8
%                                        9 10 11 12 13 14 15 16
%                                       17 18 19 20 21 22 23 24
%                                       25 26 27 28 29 30 31 32
%                                       33 34 35 36 37 38 39 40
%                                       41 42 43 44 45 46 47 48
%                                       49 50 51 52 53 54 55 56

% Anonymous function definitions
%--------------------------------------------------------------------------
IS = @atcm.fun.integrator_pst;     % wraps integrator, return only timesers
f  = @(x) IS(spm_unvec(x,P),M,U);  % wraps above for 1 vector input call
g  = @(y) y(ix,:);                 % wraps output for specifc output vectors
gf = @(x) g(f(x));                 % wraps the 2 above into 1 f(x) call

% Allow passing parameter vector also
if nargin >= 6 && ~isempty(px)
    fprintf('Initialised user p-vec\n');
    P = spm_unvec(px,P);
end

% CALL THE FUNCTION
%--------------------------------------------------------------------------
Y = gf(spm_vec(P))';               

Y = Y(181:end,:);

% Plot it
%--------------------------------------------------------------------------
if doplot
    figure('Name','LYAP','Color',[.3 .3 .3],'InvertHardcopy','off','position',[895 170 1134 997]);
    set(gcf, 'MenuBar', 'none');
    set(gcf, 'ToolBar', 'none');

    s  = axes;
    p  = plot3( Y(:,1) , Y(:,2) , Y(:,3)  , 'r' , 'linewidth' , 1 );
    cf = @hot;
    cb = [uint8(cf(length(Y))*255) uint8(ones(length(Y),1))]';
    drawnow;set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cb);
    axis square;
    view(3);
    grid on;

    s(1).YColor = [1 1 1];
    s(1).XColor = [1 1 1];
    s(1).ZColor = [1 1 1];
    s(1).Color  = [.3 .3 .3];

    if nargin < 5 || isempty(label)
        label = {'x' 'y' 'z'};
    end
    xlabel(label{1});
    ylabel(label{2});
    zlabel(label{3});
    rotate3d();
end

% Compute Lyapunov
%--------------------------------------------------------------------------
R2   = 1;
EEM  = real(Y);
rEEM = size(Y,1);
dd   = pdist(EEM,'chebychev');
dd   = squareform(dd);

mad = std(Y(:));
dd  = dd+eye(rEEM)*10*mad;

for k=0:20
    for n=1:rEEM-k
        
        l1=find(0.05*(1/R2)*mad<dd(n,1:end-k)<0.1*(1/R2)*mad)';
        
        u=dd(l1+k,n+k);
        LL(n,1) = log(mean(u));
        
    end
    L(k+1,1)=nanmean(LL);
    K(k+1,1)=k; 
end

lambda=diff(L)./diff(K);

% figure('name','Lyapunov Exponent','NumberTitle','off')
% plot(K,L,'.');
% title(['Lyapunov Exponent'])

%_________________Nonlinear Regression Layapunov Exponents_________________

Lmax=max(L);
L0=L(1);
Lm=L0+0.9*(Lmax-L0);
Ldiff=abs(L-Lm);

Tl= find(Ldiff==min(Ldiff));
x = K(1:Tl);
[betar]=regress(L(1:Tl), [ones(Tl,1) x]);

for iii=1:100
    beta = nlinfit(K(1:Tl),L(1:Tl),@nonlin1,[betar;randn(1,1)]);
    LLE1(iii,1)=beta(2,1);
end
LLE=nanmax(LL);

if doplot
    % appnd title to plot!
    title(sprintf('Maximum Lyapunov Exponent = %d',LLE),'color','w');
    set(findall(gcf,'-property','FontSize'),'FontSize',18);
end

end

function yhat = nonlin1(beta,x)
b1 = beta(1);
b2 = beta(2);
b3 = beta(3);
yhat =b1+b2*x+b3*x./exp(b2*x);
end

