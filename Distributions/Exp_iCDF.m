function x = Exp_iCDF(u, lambda)
%Compute the inverse of the cumulative distribution function of Exponential 
%random variable with scale parameter lambda -> default = (lambda = 1)

if nargin < 2
    lambda = 1;
end

x = -log(1-u)/lambda;

end