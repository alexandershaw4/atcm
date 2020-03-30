function vecx2 = VecRetainDim(x,varargin)
% make 2 x n matrix by vectorising around a dimention of x
% AS 2015

try dim = varargin{1}; catch dim = 1; end


sd = 1: ndims(x); 
x  = permute(x,[dim sd(sd~=dim)]);

for i = 1:size(x,1)
    vecx2(i,:)=spm_vec(x(i,:,:));
end
