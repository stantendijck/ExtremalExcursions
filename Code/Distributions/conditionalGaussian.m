function [condM,condS] = conditionalGaussian(M,S,inds)
% Computes the distribution of Z(inds)|Z(~inds), where Z ~ N(M,S)

n = length(M);
negInds = nan(n-length(inds),1);
cnt = 1;
for i = 1:n
    if ~any(inds == i)
        negInds(cnt) = i;
        cnt = cnt+1;
    end
end

S11 = S(inds,inds);
S12 = S(inds,negInds);
S21 = S(negInds,inds);
S22 = S(negInds,negInds);
S22_inv = S22^(-1);
M1 = M(inds);
M2 = M(negInds);

condM = @(a)( M1 + S12 * S22_inv * (a' - M2) );
condS = S11 - S12 * S22_inv * S21;

end