function y = HighResMeanFilt(mIn,n,m)
% Mean window filtering & matrix sampling
% mIn in an n dimensional matrix
% n   is a resampling value (1 = no resamp, .5 = half, 2 = double etc).
% m   is the window kernel size
%
% AS2016 [updt]

if n == 0; n = 1; end

if ndims(mIn) > 2; 
    % handle high dimentions
    mOut = size(mIn);
    mIn  = VecRetainDim(mIn,1);
end


M = imresize(full(mIn),n);
y = smoothmat(M,m,m);


try mOut;
    % reinstate orginal dimentions
    y = spm_unvec(y, ones([mOut(1)*(n) mOut(2)*(n) mOut(3)]));
end

end

function matrixOut = smoothmat(matrixIn,Nr,Nc)
% Smooths 2D array data.  Ignores NaN's.
if nargin < 2, error('Not enough input arguments!'), end

N(1) = Nr; 
N(2) = Nc;

if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

[row,col] = size(matrixIn);
eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);

A = isnan(matrixIn);
matrixIn(A) = 0;

nrmlize = eL*(~A)*eR;
nrmlize(A) = NaN;

matrixOut = eL*matrixIn*eR;
matrixOut = matrixOut./nrmlize;

end

