function F = GPD_CDF(x, mu, sigma, xi)

F = 1 - (1 + xi .* (x - mu)./sigma) .^(-1./xi);
I = x < mu;
F(I) = 0;
if any(xi < 0)
    I = xi < 0 & x > mu - sigma./xi;
    F(I) = 1;
end
    

end