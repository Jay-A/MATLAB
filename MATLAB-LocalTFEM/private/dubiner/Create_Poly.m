function [f] = Create_Poly(poly_coeffs)
%CREATE_POLY Summary of this function goes here
%   Detailed explanation goes here
n = numel(poly_coeffs);
power = 0;
f = @(x) 0;

for k= n:-1:1
    f = @(x) f(x) + poly_coeffs(k).*x.^power;
    power = power + 1;
end

