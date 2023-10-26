function y = fromLapl2Exp(x)

I = x < 0;

y = nan(size(x));
y(I) = - log(1 - exp(x(I)) / 2);
y(~I) = log(2) + x(~I);

end