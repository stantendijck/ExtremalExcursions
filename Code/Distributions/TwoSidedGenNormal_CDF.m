function F = TwoSidedGenNormal_CDF(x, m, a, b, p)

if nargin < 2
    m = 0;
end
if nargin < 3
    a = 1;
end
if nargin < 4
    b = 2;
end

F = nan(size(x));

F(x>=m) = (1-p) + 2 * p * (GenNormal_CDF(x(x>=m), m, a(1), b(1))-1/2);
F(x<m) = 2 * (1-p) * GenNormal_CDF(x(x<m), m, a(2), b(2));


end