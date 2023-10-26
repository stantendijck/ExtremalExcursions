function F = GPD_Survivor(x, mu, sigma, xi)

F = (1 + xi * (x - mu)/sigma) .^(-1/xi);
I = x < mu;
F(I) = 1;
if xi < 0
    I = x > mu - sigma/xi;
    F(I) = 0;
end
    

end