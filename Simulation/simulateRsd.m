function nextRsd = simulateRsd(prevRsd, allRsd, h)
% Simulate conditional residual given k0 previous residuals
% INPUT:
% prevRsd,      previous residuals: length = k0 < k
% allRsd,       all residuals: size = (n, k)
% h,            bandwidth matrix in kernel density estimation
%
% OUTPUT:
% nextRsd,      a (k x 1) vector containing exactly k residuals conditional
%               on first k0

% Parameters
k0 = length(prevRsd);




flagCopula = false;
if isstruct(allRsd)
    flagCopula = true; % t copula
    k = size(allRsd.rhoHAT,1);
else
    [n,k] = size(allRsd);
end

if ~flagCopula && nargin < 3
    % Silverman's rule of thumb in bandwidth selection
    h = (4/(k+2))^(2/(k+4)) * n^(-2/(k+4)) * diag(var(allRsd));
end

% Pre allocation
nextRsd = nan(k,1);
nextRsd(1:k0) = prevRsd;

%% Simulate the kth residual via kernel density simulation or from copula
if flagCopula
    if any(prevRsd)
        condMu = @(x)(allRsd.rhoHAT(k0+1:end,1:k0) * (allRsd.rhoHAT(1:k0,1:k0) \ x));
        condRho = allRsd.rhoHAT(k0+1:end,k0+1:end) - allRsd.rhoHAT(k0+1:end,1:k0) * (allRsd.rhoHAT(1:k0,1:k0) \ allRsd.rhoHAT(1:k0,k0+1:end));
        d1 = @(x)(x' * (allRsd.rhoHAT(1:k0,1:k0) \ x));
        nextRsdT = condMu(prevRsd)' + mvtrnd(condRho * (allRsd.nuHAT + d1(prevRsd))/(allRsd.nuHAT + k0), allRsd.nuHAT + k0, 1);
        for j = k0+1:length(nextRsd)
            nextRsdU = tcdf(nextRsdT(j-k0), allRsd.nuHAT);
            if isfield(allRsd,'margHAT')
                nextRsd(j) = TwoSidedGenNormal_iCDF(nextRsdU, allRsd.margHAT(j,1), allRsd.margHAT(j,2:3), allRsd.margHAT(j,4:5), allRsd.margHAT(j,6));
            elseif isfield(allRsd,'data')
                nextRsd(j) = ksdensity(allRsd.data(:,j),nextRsdU,'Function','icdf');
            elseif isfield(allRsd,'icdf')
                nextRsd(j) = myinterp1(allRsd.icdf{j}(:,1),allRsd.icdf{j}(:,2),nextRsdU);
            end
        end
    else
        nextRsdT = mvtrnd(allRsd.rhoHAT, allRsd.nuHAT, 1);
        for j = k0+1:length(nextRsd)
            nextRsdU = tcdf(nextRsdT(j-k0), allRsd.nuHAT);
            if isfield(allRsd,'margHAT')
                nextRsd(j) = TwoSidedGenNormal_iCDF(nextRsdU, allRsd.margHAT(j,1), allRsd.margHAT(j,2:3), allRsd.margHAT(j,4:5), allRsd.margHAT(j,6));
            elseif isfield(allRsd,'data')
                nextRsd(j) = ksdensity(allRsd.data(:,j),nextRsdU,'Function','icdf');
            elseif isfield(allRsd,'icdf')
                nextRsd(j) = myinterp1(allRsd.icdf{j}(:,1),allRsd.icdf{j}(:,2),nextRsdU);
            end
        end
    end
else
    if k0 ~= 0
        % Compute weights
        logw = sum(-(nextRsd(1:k0) - allRsd(:,1:k0)').^2 ./ (2*diag(h(1:k0,1:k0))),1)';
        if sum(logw==-inf)==n % if sparsity of data induces problems, increase bandwidth
            nextRsd = simulateRsd(prevRsd, allRsd, 2*diag(h));
            return
        else
            a = max(logw);
            w = exp(-a+logw)/sum(exp(-a+logw));
        end
        
        %% Simulate from the conditional distribution
        % Sample a distribution from the mixture of distribution
        %         try
        index = randsample(n,1,true,w);
        %         catch
        %             de = 1;
        %         end
    else
        index = randsample(n,1);
    end

    % Sample from the sampled distribution
    nextRsd(k0+1:k) = allRsd(index,k0+1:k)' + sqrt(diag(h(k0+1:k,k0+1:k)))*randn(1,1);
end

end