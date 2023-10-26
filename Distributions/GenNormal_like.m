function NLOGL = GenNormal_like(x, m, a, b)

if a < 0 || b < 0
    NLOGL = inf;
    return
end

n = length(x);

NLOGL = -n*log(b/(2*a*gamma(1/b))) + sum((abs(x-m)/a).^b);

end
