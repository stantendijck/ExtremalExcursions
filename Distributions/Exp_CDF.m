function F = Exp_CDF(x, lambda)
%Compute the cumulative distribution function of an Exponential random variable
%with scale parameter lambda -> default = (lambda = 1)

if nargin < 2
    lambda = 1;
end

F = 1 - exp(-lambda*max(x,0));

end