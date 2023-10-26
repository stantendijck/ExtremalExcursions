function NLOGL = HTn_like(Y, X, t, p, logX, MLEs, keefProp)


alp = p(1:2*(t-1));
bet = p(2*t-1);

if any(bet>1)
    NLOGL = inf;
    return
end

% Reparameterization: alp(3) = (mle2 - alpha(2))/mle1 + alpha(3)

xBeta = exp(logX(:,1) .* bet);
newAlp = nan(size(alp));
newAlp(1) = alp(1);
newAlp(2) = (MLEs{end}(1) - newAlp(1))/MLEs{1}(1) + alp(2);
for i = 3:2*(t-1)
    newAlp(i) = -alp(i-1) * MLEs{i-2}(1)/MLEs{i-1}(1) + alp(i);
end

Rsd1 = (Y - X * newAlp) ./ xBeta(:,1);

% true alph
alph =  MLEs{end}(1) + alp(end) * MLEs{end-1}(1);
Rsd2 =  (Y - alph(1) * X(:,1)) ./ xBeta(:,1);

f = keefConstraints(X(:,1), Y, [alph;bet], [0,1], xBeta(:,1), Rsd2, keefProp);
flag_keef = all(f);

if ~flag_keef
    NLOGL = inf;
    return
end

n = size(X,1);
mu2 = sum(Rsd1, 1) / n;
sig2 = sum(abs(Rsd1 - mu2) .^ 2, 1) / (n - 1);
% sig2 = std(Rsd1);

NLOGL = n/2 * log(2*pi) + n/2*log(sig2) + sum(log(xBeta(:,1))) + sum((Rsd1-mu2).^2/(2*sig2), 1);


end