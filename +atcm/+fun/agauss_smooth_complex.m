function y = agauss_smooth_complex(x,n,a)

if nargin < 3 || isempty(a)
    a = 1;
end


yr = agauss_smooth(real(x),n,a);
yi = agauss_smooth(imag(x),n,a);

y = yr + 1i * yi;

end