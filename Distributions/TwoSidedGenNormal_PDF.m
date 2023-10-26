function f = TwoSidedGenNormal_PDF(x,m,a,b,p)

if any(a < 0) || any(b < 0) || p < 0 || p > 1
    error('Parameters of TwoSidedGenNormal_PDF do not satisfy constraints');
end

f = nan(size(x));

f(x>=m) = p * b(1)/(a(1)*gamma(1/b(1))) * exp(-(abs(x(x>=m)-m)/a(1)).^b(1));
f(x<m) = (1-p) * b(2)/(a(2)*gamma(1/b(2))) * exp(-(abs(x(x<m)-m)/a(2)).^b(2));




end