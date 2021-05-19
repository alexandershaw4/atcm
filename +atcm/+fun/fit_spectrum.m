function O = fit_spectrum(w,y,priorF)
% an *approximate* implementation of the fooof algorithm
%
% breaks a spectrum of x-axis values, w, and y-axis values, y;
% (e.g. plot(w,y) )
%
% into,   y = aperiodic + (skew) gmm
%
% where,  - aperiodic is a 2nd order poly with fixed start/end points
%         - spectral peaks are modelled as a n-dim GMM with starting
%           positions priorF
%
% usage:
%        outputs = fit_spectrum(w,y,priorF)
%              w = freq vec
%              y = spectrum
%         priorF = vector of initial guesses for peak locs.
%                  for visual gamma spectrum, try [10 15 40 60]
%
% output struture:
%  out.aperiodic = aperiodic component, generated from out.A*out.b
%  out.error     = fitted model squared error (SSE)
%  out.yinput    = original, inputted spectrum
%  out.yprime    = spectrum fit (i.e. the full model: aperiod+gmm)
%  out.Peaks     = the peak frequency, amplutude, width and skew of
%                  each gaussian bump / peak
%  out.c         = the individual components of the gmm
%  out.priorfit  = the initial (not final) fit of the aperiodic part
%
% AS21
%
% For fooof proper (python) see:
% Donoghue et al, 2020 NatNeuro: https://fooof-tools.github.io/fooof/


w = w(:);
y = smooth(y(:));
yorig = y;


% (1.) Fit aperiodic part as a 2nd order poly with fixed start and end points
%--------------------------------------------------------------------------
iC = [1 length(y)];
N  = 2;
A  = w.^[0:N]; 
C  = A(iC,:);                   % and the constrained points subset
b  = lse(A,y,C,y(iC));
aper = A*b;

% remove aperiodic component
y = y - aper;

% return first guess parts
O.priorfit.aperiodoc = aper;
O.priorfit.A = A;
O.priorfit.b = b;

% (2.) Next fit the remianing data with parameterised skewed gaussians
%--------------------------------------------------------------------------
[err,X1,P,fx] = fitmodel(priorF,A,b,N,w,y,yorig);

% Construct output structure
O.aperiodic = A*spm_vec(X1(1:N+1));
O.error    = err;
O.Peaks    = spm_unvec(X1(N+2:end),P);
O.yprime   = fx(X1);
O.A        = A;
O.b        = spm_vec(X1(1:N+1));
O.w        = w;
O.yinput   = yorig;
O.w        = w;

P=[];
for i = 1:length(priorF)
    P.Freq   = O.Peaks.Freq(i) - w(1);
    P.Amp    = O.Peaks.Amp(i);
    P.Wid    = O.Peaks.Wid(i);
    P.Skew   = O.Peaks.Skew(i);
    O.c(i,:) = atcm.fun.makef(w,P);    
end

O.notes = {'model is aperiod + gmm. aper = A*b and gmm = c'};

end


function [err,X1,P,fx] = fitmodel(priorF,A,b,N,w,y,yorig)

NPK=length(priorF);
for i = 1:length(priorF)
    LOCS(i) = findthenearest(w,priorF(i));
end

PKS = y(LOCS);

% build a skew gauss GMM & fit it
P.Freq = w(LOCS);
P.Amp  = PKS;
P.Wid  = 2*ones(NPK,1);
P.Skew = zeros(NPK,1);

% skew function (see sub funcs)
x0 = (spm_vec(P));
f  = @(x,w) (makegauskew(w,spm_unvec((x),P)));
g  = @(x) sum( (y - f(x,w)).^2 );

opts  = optimset('Display','off');
[X,F] = fminsearch(g,x0,opts);

% save initial fits
O.priorfit.x0 = X;

% build a combination model: slope + gmm == A*b + f(x,w)
px = [b(:); X(:)];
fx = @(x) (A*spm_vec(x(1:N+1))) + f(x(N+2:end),w) ;
gx = @(x) sum( (yorig - fx(x)).^2 );

[X1,F1] = fminsearch(gx,px,opts);

ffx = @(x,w) (A*spm_vec(x(1:N+1))) + f(x(N+2:end),w) ;
[X1,resnorm,residual,~,~,~,J] = lsqcurvefit(ffx,X1,w,yorig,[],[]);
err = sum(residual.^2);

end



function X = makegauskew(w,Fq,Amp,Wid,Skew)
%
% afit.makef(w,Fq,Amp,Wid)
%
% e.g.
% w = 1:100
% S = afit.makef(w,[10 20 50 70],[10 8 5 3],[1 2 5 3]);
% figure;plot(w,S)
%
% AS

if nargin==2 && isstruct(Fq)
   w    = w;
   Skew = Fq.Skew;
   Wid  = Fq.Wid;
   Amp  = Fq.Amp;
   Fq   = Fq.Freq;
end

if length(Fq) > 1
    for i = 1:length(Fq)
        try
            X0 = X0 + makegauskew(w,Fq(i),Amp(i),Wid(i),Skew(i));
        catch
            X0 =      makegauskew(w,Fq(i),Amp(i),Wid(i),Skew(i));
        end
        %X0(i,:) = afit.makef(w,Fq(i),Amp(i),Wid(i));

    end
    %X0 = max(X0); 
    X  = X0;
    return;
end

try Wid ; catch Wid = 2; end
try Amp ; catch Amp = 2; end

% offset negative values in w
mw  = min(w);
X   = 0*w;
f   = findthenearest(Fq,w);
f   = f(1);
w   = w - mw;
X   = X + Amp * exp( -(w-f).^2 / (2*(2*Wid)^2) );
w   = w + mw;

% Apply skew
X = 2*X.*normcdf(Skew*w);
end




function [r,c,V] = findthenearest(srchvalue,srcharray,bias)

if nargin<2
    error('Need two inputs: Search value and search array')
elseif nargin<3
    bias = 0;
end

% find the differences
srcharray = srcharray-srchvalue;

if bias == -1   % only choose values <= to the search value
    
    srcharray(srcharray>0) =inf;
        
elseif bias == 1  % only choose values >= to the search value
    
    srcharray(srcharray<0) =inf;
        
end

% give the correct output
if nargout==1 | nargout==0
    
    if all(isinf(srcharray(:)))
        r = [];
    else
        r = find(abs(srcharray)==min(abs(srcharray(:))));
    end 
        
elseif nargout>1
    if all(isinf(srcharray(:)))
        r = [];c=[];
    else
        [r,c] = find(abs(srcharray)==min(abs(srcharray(:))));
    end
    
    if nargout==3
        V = srcharray(r,c)+srchvalue;
    end
end

end
    
function x = lse(A,b,C,d,solverflag,weights)
% solves A*x = b for X in a least squares sense, given C*x = d
% usage: x = lse(A,b,C,d)
% usage: x = lse(A,b,C,d,solverflag)
% usage: x = lse(A,b,C,d,solverflag,weights)
%
% Minimizes norm(A*x - b),
% subject to C*x = d
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 3.0
% Release date: 1/31/07

% check sizes
[n,p] = size(A);
[r,nrhs] = size(b);
[m,ccols] = size(C);
if n~=r
  error 'A and b are incompatible in size (wrong number of rows)'
elseif ~isempty(C) && (p~=ccols)
  error 'A and C must have the same number of columns'
elseif ~isempty(C) && issparse(C)
  error 'C may not be a sparse matrix'
elseif ~isempty(C) && (m~=size(d,1))
  error 'C and d are incompatible in size (wrong number of rows)'
elseif ~isempty(C) && (size(d,2)~=1)
  error 'd must have only one column'
elseif isempty(C) && ~isempty(d)
  error 'C and d are inconsistent with each other (one was empty)'
elseif ~isempty(C) && isempty(d)
  error 'C and d are inconsistent with each other (one was empty)'
end

% default solver is '\'
if (nargin<5) || isempty(solverflag)
  solverflag = '\';
elseif ~ischar(solverflag)
  error 'If supplied, solverflag must be character'
else
  % solverflag was supplied. Make sure it is legal.
  valid = {'\', 'backslash', 'pinv'};
  ind = strmatch(solverflag,valid);
  if (length(ind)==1)
    solverflag = valid{ind};
  else
    error(['Invalid solverflag: ',solverflag])
  end
end

% default for weights = []
if (nargin<6) || isempty(weights)
  weights = [];
else
  weights = weights(:);
  if (length(weights)~=n) || any(weights<0)
    error 'weights should be empty or a non-negative vector of length n'
  elseif all(weights==0)
    error 'At least some of the weights must be non-zero'
  else
    % weights was supplied. scale it to have mean value = 1
    weights = weights./mean(weights);
    % also sqrt the weights for application as an
    % effective replication factor. remember that
    % least squares will minimize the sum of squares.
    weights = sqrt(weights);
  end
end

% tolerance used on the solve
Ctol = 1.e-13;

if (nargin<3) || isempty(C)
  % solve with \ or pinv as desired.
  switch solverflag
    case {'\' 'backslash'}
      % solve with or without weights
      if isempty(weights)
        x = A\b;
      else
        x = (repmat(weights,1,size(A,2)).*A)\ ...
          (repmat(weights,1,nrhs).*b);
      end
      
    case 'pinv'
      % solve with or without weights
      if isempty(weights)
        ptol = Ctol*norm(A,1);
        x = pinv(A,ptol)*b;
      else
        Aw = repmat(weights,1,size(A,2)).*A;
        ptol = Ctol*norm(Aw,1);
        x = pinv(Aw,ptol)*(repmat(weights,1,nrhs).*b);
      end
  end
  
  % no Constraints, so we are done here.
  return
end

% Which solver do we use?
switch solverflag
  case {'\' 'backslash'}
    % allow a rank deficient equality constraint matrix
    % column pivoted qr to eliminate variables
    [Q,R,E]=qr(C,0);
    
    % get the numerical rank of R (and therefore C)
    if m == 1
%      rdiag = R(1,1);
      rdiag = abs(R(1,1));
    else
      rdiag = abs(diag(R));
    end
    crank = sum((rdiag(1)*Ctol) <= rdiag);
    if crank >= p
      error 'Overly constrained problem.'
    end
    
    % check for consistency in the constraints in
    % the event of rank deficiency in the constraint
    % system
    if crank < m
      k = Q(:,(crank+1):end)'*d;
      if any(k > (Ctol*norm(d)));
        error 'The constraint system is deficient and numerically inconsistent'
      end
    end
    
    % only need the first crank columns of Q
    qpd = Q(:,1:crank)'*d;
    
    % which columns of A (variables) will we eliminate?
    j_subs = E(1:crank);
    % those that remain will be estimated
    j_est = E((crank+1):p);
    
    r1 = R(1:crank,1:crank);
    r2 = R(1:crank,(crank+1):p);
    
    A1 = A(:,j_subs);
    qpd = qpd(1:crank,:);
    
    % note that \ is still ok here, even if pinv
    % is used for the main regression.
    bmod = b-A1*(r1\repmat(qpd,1,nrhs));
    Amod = A(:,j_est)-A1*(r1\r2);
    
    % now solve the reduced problem, with or without weights
    if isempty(weights)
      x2 = Amod\bmod;
    else
      x2 = (repmat(weights,1,size(Amod,2)).*Amod)\ ...
        (repmat(weights,1,nrhs).*bmod);
    end
    
    % recover eliminated unknowns
    x1 = r1\(repmat(qpd,1,nrhs)-r2*x2);
    
    % stuff all estimated parameters into final vector
    x = zeros(p,nrhs);
    x(j_est,:) = x2;
    x(j_subs,:) = x1;
    
  case 'pinv'
    % allow a rank deficient equality constraint matrix
    Ctol = 1e-13;
    % use svd to deal with the variables
    [U,S,V] = svd(C,0);
    
    % get the numerical rank of S (and therefore C)
    if m == 1
      sdiag = S(1,1);
    else
      sdiag = diag(S);
    end
    crank = sum((sdiag(1)*Ctol) <= sdiag);
    if crank >= p
      error 'Overly constrained problem.'
    end
    
    % check for consistency in the constraints in
    % the event of rank deficiency in the constraint
    % system
    if crank < m
      k = U(:,(crank+1):end)'*d;
      if any(k > (Ctol*norm(d)));
        error 'The constraint system is deficient and numerically inconsistent'
      end
    end
    
    % only need the first crank columns of U, and the
    % effectively non-zero diagonal elements of S.
    sinv = diag(S);
    sinv = diag(1./sinv(1:crank));
    % we will use a transformation
    %  Z = V'*X = inv(S)*U'*d
    Z = sinv*U(:,1:crank)'*d;
    
    % Rather than explicitly dropping columns of A, we will
    % work in a reduced space as defined by the svd.
    Atrans = A*V;
    
    % thus, solve (A*V)*Z = b, subject to the constraints Z = supd
    % use pinv for the solution here.
    ptol = Ctol*norm(Atrans(:,(crank+1):end),1);
    if isempty(weights)
      Zsolve = pinv(Atrans(:,(crank+1):end),ptol)* ...
        (b - repmat(Atrans(:,1:crank)*Z(1:crank),1,nrhs));
    else
      w = spdiags(weights,0,n,n);
      Zsolve = pinv(w*Atrans(:,(crank+1):end),ptol)* ...
        (w*(b - repmat(Atrans(:,1:crank)*Z(1:crank),1,nrhs)));
    end
    
    % put it back together in the transformed state
    Z = [repmat(Z(1:crank),1,nrhs);Zsolve];
    
    % untransform back to the original variables
    x = V*Z;
    
end
end


















% M = fit(w,y,'Fourier6');
% 
% for i = 1:6
%     comp(i,:) = eval(['M.a' num2str(i) '*cos(i*w*M.w) + M.b' num2str(i) '*sin(i*w*M.w)']);
% end
% 
% c = M.a0;

% p0 = ones(1,12);
% %p0=rand(size(p0));
% 
% % [PKS,LOCS] = findpeaks(y,'NPeaks',8);
% % 
% % p0([2 5 8 11 14 17 20 23] = w(LOCS);
% % 
% opts = optimoptions('lsqcurvefit');
% opts.OptimalityTolerance = 1e-12;
% opts.MaxFunctionEvaluations = 1e10;
% opts.FunctionTolerance=1e-12;
% opts.MaxIterations=1e8;
% % 
% [X,resnorm,residual,~,~,~,J] = lsqcurvefit(@lfun3,p0,w,y,[],[],opts);
% 
% xb = X;
% for i = 1:6
%     xx = xb(1:3);
%     
%     c(i,:) = xx(1)./((w-xx(2)).^2+xx(3));
%     
%     xb(1:3)=[];
% end
% 
% 
% 
% %f = @(x) lfun3(x,w);
% %g = @(x) sum( (y - f(x)).^2 );
% %[X,F] = fminsearch(g,p0)
% 
%         
% %yp = lfun3(X,w);


% 
% function F = lfun3(p,x)
%     F = p(1)./((x-p(2)).^2+p(3)) + p(4)./((x-p(5)).^2+p(6)) + ...
%         p(7)./((x-p(8)).^2+p(9)) + p(10)./((x-p(11)).^2+p(12));% + ...
%         %p(13)./((x-p(14)).^2+p(15)) + p(16)./((x-p(17)).^2+p(18)) ;;% + ...
%       %  p(19)./((x-p(20)).^2+p(21)) + p(22)./((x-p(23)).^2+p(24)) ;

 