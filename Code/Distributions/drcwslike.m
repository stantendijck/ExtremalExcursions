function NLOGL = drcwslike(y, phi, lambda)

if any(lambda <= 0)
    NLOGL = inf;
    return
end

n = length(y(:,1));
sigma2 = lambda(1) + lambda(2) * exp(-lambda(3) * y(:,1));

NLOGL = n/2*log(2*pi) + sum(log(sigma2))/2 + sum((y(:,1) - y(:,2:end)*phi').^2./(2*sigma2));


end