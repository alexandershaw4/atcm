function [j,ip] = jaco(DCM,pE)
% Compute the Jacobian matrix (parameter gradients) for a DCM model using:
%
% j(ip,:) = ( IS( f(x,P(ip)+h) ) - IS( f(x,P(ip)-h) ) ) / (2 * h)
%
% IS is the integrator stored in DCM.M.IS
% f  is the motion equations in DCM.M.f
% P  is a param structure, 2nd input to this function, and as-per DCM.M.pE
% h  is determined as 1% of pE(ip)
%
% j  is the outputted matrix
%
% Calculated for parameters of variance >0 only - determined by DCM.M.pC
%
% AS2019

Ep = spm_vec(pE);
IS = DCM.M.IS;
P  = DCM.M.pE;
V  = spm_vec(DCM.M.pC);
ip = (~~V);

j  = jacf(IS,Ep,P,DCM.M,DCM.xU,ip);

%j(isnan(j)) = 0;

end



function j = jacf(IS,Ep,P,M,U,ip)

% Compute the Jacobian matrix using variable step-size
n  = 0;
warning off ;

fx    = spm_vec(feval(IS,P,M,U));
j     = zeros(length(Ep),length(spm_vec(fx)));
for i = 1:length(Ep)
    if ip(i)
        
        % Print progress
        n = n + 1;
        if n > 1; fprintf(repmat('\b',[1,length(str)])); end
        str  = sprintf('Computing Gradients [ip %d / %d]',n,length(find(ip)));
        fprintf(str);
        
        % Compute Jacobi: A(j,:) = ( f(x+h) - f(x-h) ) / (2 * h)
        P0     = Ep;
        P1     = Ep;
        d      = P0(i) * 0.01;
        P0(i)  = P0(i) + d  ;
        P1(i)  = P1(i) - d  ;
                
        f0     = spm_vec(feval(IS,spm_unvec(P0,P),M,U));
        %f1     = spm_vec(feval(IS,spm_unvec(P1,P),M,U));
        
        %j(i,:) = (f0 - f1) / (2 * d);
        j(i,:) = (f0 - fx) / (d);
        
        % Alternatively, include curvature
        %deriv1 = (f0 - f1) / 2 / d;
        %deriv2 = (f0 - 2 * fx + f1) / d ^ 2;
        %j(i,:) = deriv1 ./ deriv2;
    end
end

warning on;
fprintf('\n');
end