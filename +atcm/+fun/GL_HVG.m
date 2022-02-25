function GL = GL_HVG(x,w)
% returns the Graph Laplacian of the Horizontal Visibility Graph Adjacency
% matrix:
%    GL = atcm.fun.GL_HVG(x,w)
% AS22

Q  = fast_HVG(x,w);
A  = Q .* ~eye(length(Q));
N  = size(A,1);
GL = speye(N,N) + (A - spdiags(sum(A,2),0,N,N))/4;