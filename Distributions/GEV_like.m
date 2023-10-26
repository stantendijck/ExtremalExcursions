function NLOGL = GEV_like(x, mu, sigma, xi)

if (xi < 0 && max(x) > mu - sigma/xi) || (xi > 0 && min(x) < mu - sigma/xi)
    NLOGL = inf;
    return
end

if sigma < 0
    NLOGL = inf;
    return;
end

n = length(x);
tx = (1 + xi * (x - mu) / sigma).^(-1/xi);
NLOGL = n*log(sigma) - (xi + 1) * sum(log(tx)) + sum(tx);

end