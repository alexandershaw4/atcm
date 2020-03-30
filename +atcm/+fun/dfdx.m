function j = dfdx(DCM,P,order)
% STATE differentiation
%
%
% AS

if nargin <3 || isempty(order)
    order = 1;
end

IS = spm_funcheck(DCM.M.f);
M  = DCM.M;

warning off ;

fx    = spm_vec(feval(IS,spm_vec(M.x),0,P,M));
j     = complex(zeros(length(spm_vec(fx)),length(spm_vec(fx))));
for i = 1:length(fx)
                
        % Compute Jacobi: A(j,:) = ( f(x+h) - f(x-h) ) / (2 * h)
        x0     = fx * 1i;
        x1     = fx * 1i;
        d      = fx(i) * 0.01;
        x0(i)  = x0(i) + d   ;
        x1(i)  = x1(i) - d   ;
                
        f0     = spm_vec(feval(IS,x0,0,P,M));    
        f1     = spm_vec(feval(IS,x1,0,P,M));
        j(i,:) = (f0 - f1) / (2 * d);
        
        if order ==2
        % Alternatively, include curvature
            deriv1 = (f0 - f1) / 2 / d;
            deriv2 = (f0 - 2 * fx + f1) / d ^ 2;
            j(i,:) = deriv1 ./ deriv2;
        end
end

j           = j';
j(isnan(j)) = 0;

warning on;
end