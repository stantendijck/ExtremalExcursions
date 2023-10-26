function NLOGL = drchslike2(x, phi, lambda)

if any(lambda <= 0)
    NLOGL = inf;
    return
end

if size(phi,1) > 1
    phi = phi';
end

n = length(x(:,1));
mu = x(:,2:end)*phi';
sigma2 = lambda(1) + lambda(2) * exp(-lambda(3) * x(:,1));

NLOGL = n/2*log(2*pi) + sum(log(sigma2))/2 + sum((x(:,1) - mu).^2./(2*sigma2));


end