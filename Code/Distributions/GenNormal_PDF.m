function u = GenNormal_PDF(x, m, a, b)
%Compute the PDF of a Laplace random variable
%with parameters (mu, b) -> default = (mu = 0, b = 1)

if nargin < 2
    m = 0;
end
if nargin < 3
    a = 1;
end
if nargin < 4
    b = 1;
end

u = b / (2*a*gamma(1/b)) * exp(-(abs(x-m)/a).^b);



end