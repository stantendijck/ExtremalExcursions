function F = GenNormal_CDF(x, m, a, b)

if nargin < 2
    m = 0;
end
if nargin < 3
    a = 1;
end
if nargin < 4
    b = 2;
end


F = 1/2 + sign(x-m)/2 .* gammainc(abs((x-m)/a).^b, 1/b);
    

end