function NLOGL = TwoSidedGenNormal_like(x,m,a,b,p)

if any(a < 0) || any(b < 0) || p < 0 || p > 1
    NLOGL = inf;
    return
end

xExc = x(x>=m); nExc = length(xExc);
xBel = x(x<m); nBel = length(xBel);

NLOGL = - nExc * log(p * b(1)/(a(1)*gamma(1/b(1)))) + sum((abs(x(x>=m)-m)/a(1)).^b(1));
NLOGL = NLOGL - nBel*log((1-p) * b(2)/(a(2)*gamma(1/b(2))))  + sum((abs(x(x<m)-m)/a(2)).^b(2));




end