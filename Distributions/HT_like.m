function NLOGL = HT_like(Y, X, p, logX, keefProp)

% p = [alp(1:3),bet(1:2)]
alp = p(1);
bet = p(2);
if any(abs(alp)>1) || any(bet>1)
    NLOGL = inf;
    return
end

xBeta = exp(logX * bet);
Rsd = (Y - alp(1) * X) ./ xBeta(:,1);

f = keefConstraints(X, Y, [alp(1); bet(1)], [0,1], xBeta, Rsd, keefProp);
% f = HeffernanTawn.KeefConstraints(X, Y, [alp(1); bet(1)], [0,1], xBeta);
% if ~all(f1==f)
%     b = 1;
% end
flag_keef = all(f);

if ~flag_keef
    NLOGL = inf;
    return
end

n = size(X,1);

sig2 = sum(abs(Rsd - sum(Rsd, 1) ./ n).^2, 1) ./ (n - 1);
NLOGL = n/2 * log(2*pi) + n/2 * log(sig2) + sum(log(xBeta(:,1))) + sum(Rsd.^2/(2*sig2), 1);


end