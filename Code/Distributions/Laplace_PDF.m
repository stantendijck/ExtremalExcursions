function u = Laplace_PDF(x, mu, b)
%Compute the PDF of a Laplace random variable
%with parameters (mu, b) -> default = (mu = 0, b = 1)

if nargin < 2
    mu = 0;
end
if nargin < 3
    b = 1;
end

u = 1/(2*b) * exp(-abs(x-mu)/b);



end