function x = GEV_iCDF(u, xi, mu, sigma)
%Compute the inverse of the cumulative distribution function of a GPD
%random variable with parameters (xi, mu, sigma) -> default = (mu = 0, sigma = 1)

if nargin < 3
    mu = 0;
end
if nargin < 4
    sigma = 1;
end

x = mu + sigma/xi * ((-log(u)).^(-xi) - 1);

end