function NLOGL = drchslike3(x, assvar, phi, lambda)

if any(lambda < 0)
    NLOGL = inf;
    return
end

if size(phi,1) > 1
    phi = phi';
end

n = length(x(:,1));
sigma2 = lambda(1) + lambda(2) * exp(-lambda(3) * assvar);

NLOGL = n/2*log(2*pi) + sum(log(sigma2))/2 + sum((x(:,1) - x(:,2:end)*phi').^2./(2*sigma2));


end