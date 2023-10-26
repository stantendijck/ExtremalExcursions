function u = Weibull_Survivor(x, lambda, k)
%Compute the cumulative distribution function of a Weibull random variable
%with parameters (lambda, k) -> default = (lambda = 1, k = 1)

if nargin < 2
    lambda = 1;
end
if nargin < 3
    k = 1;
end

u = exp(-(x/lambda).^k);


end