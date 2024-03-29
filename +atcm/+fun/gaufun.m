classdef gaufun < handle
    
    properties

    X
    J
    Q

    end

    methods(Static)

        
        function [J,win] = SearchGaussPCA(X,N)
            
            % transpose if req
            if any(size(X)==1) && size(X,1)==1
                X = X';
            end

            % compute N-dim PCs
            for i = 1:N
                J{i} = gaufun.GaussPCA(X,(i));
            end

            % compare
            C = cat(1,J{:});
            for i = 1:N
                %r(i) = corr(X(:),J{i}(:)).^2;
                r(i) = distcorr(C(i,:)',X);

            end
            
            [~,win] = min(r);
                                   
            J = J{win};
            
            
        end

        function J = GaussPCA(X,N)

            if nargin < 2 || isempty(N)
                N = 1;
            end

            for i = 1:size(X,2)

                QM = gaufun.AGenQn(X(:,i),8);
                %QM=(QM+QM')/2;
                [u,s,v] = svd(QM);
                J(i,:) = sum( QM*v(:,1:N), 2);

            end

            J = abs(J);

        end

        function [Q,GL] = AGenQn(x,n)
            % Convert vector into Gaussian process, i.e. a smoothed symmetric matrix form.
            %
            % Q = AGenQn(x,n)
            %
            % Implements function
            %       f = @(x,n) diag(ones(N-abs(x),1),x)/n;
            % for the sequence
            %       M = f(0,1) + f(1,2) + f(-1,2) + f(2,4) + f(-2,4) + f(3,8) + f(-3,8) ...
            %
            % AS2022

            if nargin < 2 || isempty(n)
                n = 3;
            end

            N = length(x);

            f = @(x,n) diag(ones(N-abs(x),1),x)/n;

            % a la, M = f(0,1) + f(1,2) + f(-1,2) + f(2,4) + f(-2,4) + f(3,8) + f(-3,8);

            M = f(0,1);

            for i  = 1:n
                s  = cumprod(repmat(2,[1 i]));
                M  = M + f(i,s(end)) + f(-i,s(end));
            end

            % linear model
            %Q = smooth2(denan(M.\diag(x)),4);
            Q = smooth2(M*diag(x),4);

            if nargout == 2
                %Q  = cdist(x,x) .* Q;
                A  = Q .* ~eye(length(Q));
                N  = size(A,1);
                GL = speye(N,N) + (A - spdiags(sum(A,2),0,N,N))/4;
            end
        end


    end





end