function stormOrg = simulateDrc5(storm, mdldrc, data, MMX, MMY, upcrossingInds, I)
%Transform input storm from standard margins to original margins
% INPUT:
% storm         (n x 2) vector: storm on standard margins
% mdldrc        struct(): directional mdoel
% data          struct(): data
% MM            struct(): MarginalModel


stormOrg = struct();

%% Simulate initial upcrossing associated direction
nStorms = length(data.storms);
if nargin < 5
    upcrossingInds = cellfun(@(x)(x.index(1)),data.storms);
end
assDrc = [data.Xdrc(upcrossingInds),data.Ydrc(upcrossingInds)];

if nargin < 7
    I = randsample(nStorms, 1);
end

initDrc = assDrc(I,:);

%% Transform storm-peak to original margins
% [~,iInit] = max(storms.(method){iOrder}{iStorm}(:,1));
iInit = 1;

hsOrgInit = transform2Org(storm(iInit,1), initDrc(1), MMX);
wsOrgInit = transform2Org(storm(iInit,2), initDrc(2), MMY);

if size(storm,1) ==  1
    stormOrg.data = [hsOrgInit,wsOrgInit];
    stormOrg.drcdata = initDrc;
    return
end

%% Simulate delta direction using empirical model

mysig = @(x,l)(l(1) + l(2)*exp(-l(3)*x));
sigX = mysig(hsOrgInit,mdldrc.X.phat(2:4));

% zX = (mdldrc.X.deltaX(upcrossingInds) - mdldrc.X.phat(1) * mdldrc.X.deltaX(upcrossingInds-1) - mdldrc.X.phat(2) * mdldrc.diffXY(upcrossingInds)) ./ mysig(data.X(upcrossingInds),mdldrc.X.phat(3:5));
zX = (mdldrc.X.deltaX(:,1) - mdldrc.X.phat(1) * mdldrc.X.deltaX(:,2)) ./ mysig(mdldrc.X.assX, mdldrc.X.phat(2:4));
m = sigX * mean(zX) / (1 - mdldrc.X.phat(1));
s = sqrt(sigX^2 * std(zX)^2 / (1 -  mdldrc.X.phat(1)^2));

deltaDrcX = randn(1)*s + m;

%% Calculate next direction
deltaDrcXOrg = quickquantile(mdldrc.theta.X,normcdf(deltaDrcX));
drcX = [initDrc(1),mod(initDrc(1) + deltaDrcXOrg,360)];

%% Transform next HS and WS to original margins

if size(storm,1) >= iInit+1
    hsOrg = transform2Org(storm(iInit+1,1), drcX(end), MMX);
end

%% Next use the directional model to transform the next value to original margins

zXall = (mdldrc.X.deltaX(:,1) - mdldrc.X.phat(1) * mdldrc.X.deltaX(:,2)) ./ mysig(mdldrc.X.assX,mdldrc.X.phat(2:4));

currInd = iInit + 2;
hsOrgVec = [hsOrgInit, hsOrg];
while currInd <= size(storm,1)
    %% Previous values
    latestDrcX = drcX(end);
    latestDeltaDrcX = deltaDrcX;
    latestHS = hsOrgVec(end);
    
    %% Simulate delta direction of Gaussian margins - keep dependence structure
    i = randsample(length(zXall),1);
    deltaDrcX = mdldrc.X.phat(1) * latestDeltaDrcX + mysig(latestHS, mdldrc.X.phat(2:4)) * zXall(i);
    
    %% Transform delta direction to original margins
    deltaDrcXOrg = quickquantile(mdldrc.theta.X,normcdf(deltaDrcX));
    
    %% Update direction vector
    drcX = [drcX, mod(latestDrcX + deltaDrcXOrg,360)];
    
    %% Transform HS and WS back to original margins
    hsOrgVec(end+1) = transform2Org(storm(currInd,1), drcX(end), MMX);
    
    %% Update index
    currInd = currInd + 1;
end

%% Simulate WS

[n,k] = size(mdldrc.Y.rsd);
h = (4/(k+2))^(2/(k+4)) * n^(-2/(k+4)) * diag(var(mdldrc.Y.rsd));

wsOrgVec = nan(size(hsOrgVec));

drcY = nan(size(drcX));
drcY(1) = initDrc(2);
wsOrgVec(1) = transform2Org(storm(1,2), drcY(1), MMY);

prevRsd = mod(drcX(1) - drcY(1) + 180,360) - 180;
nextRsd = simulateRsd(prevRsd,mdldrc.Y.rsd,h);

for i = 2:length(drcX)
    drcY(i) = mod(drcX(i) - nextRsd(2),360);
    wsOrgVec(i) = transform2Org(storm(i,2), drcY(i), MMY);
    nextRsd = simulateRsd(nextRsd(2), mdldrc.Y.rsd,h);
end


stormOrg.data = [hsOrgVec',wsOrgVec'];
stormOrg.drcdata = [drcX',drcY'];

