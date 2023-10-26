function stormsOrgHM = simulateHM2(nStormsSim, trainData, MM_HS, mdl, peaks)

if nargin < 5 || isempty(peaks)
    peaks = rand(nStormsSim,1);
end
%% Defaults
nStorms = length(trainData.storms);
inds = cellfun(@(x)(x.peakInd),trainData.storms);
peakDrcAll = trainData.Xdrc(inds);


%% Fit models
% Linear regression on peak HS and associated WS
expl = [ones(length(inds),1),trainData.X(inds)];
resp = trainData.Y(inds);
beta = expl\resp;
sig = std(expl*beta - resp);

stormsOrgHM = cell(1,1);
stormsOrgHM{1} = cell(nStormsSim,1);
for iStorm = 1:nStormsSim
    %% Simulate peak associated direction
    peakDrcInd = randsample(length(peakDrcAll),1);

    peakDrc = trainData.Xdrc(peakDrcInd);

    %% Simulate peak HS
    u = peaks(iStorm);
    if u > mdl.SP.X.q
        stormpeak = GPD_iCDF((u - mdl.SP.X.q)/(1-mdl.SP.X.q), mdl.SP.X.MLE(3), mdl.SP.X.MLE(1), mdl.SP.X.MLE(2));
    else
%         stormpeak = randsample(mdl.SP.X.data(mdl.SP.X.data < mdl.SP.X.MLE(1)), 1);
        stormpeak = quantile(mdl.SP.X.data(mdl.SP.X.data < mdl.SP.X.MLE(1)), u/mdl.SP.X.q);
    end
    peakHS = transform2Org(stormpeak, peakDrc, MM_HS);
%     peakHS = transform2Org(Laplace_iCDF(peak(iStorm)), peakDrc, MM_HS);

    %% Find closest storms
    distanceMatrix = nan(nStorms,1);

    a = 10; % 5 degrees in direction = 0.5m in HS
    for iStorm2 = 1:nStorms
        inds = trainData.storms{iStorm2}.index;
        [~,I] = max(trainData.X(inds));
        I = I + inds(1) - 1;
        distDrc = abs(mod(trainData.Xdrc(I) - peakDrc + 180,360) - 180);
        distHS = abs(peakHS - trainData.X(I));
        distanceMatrix(iStorm2) = a*distHS + distDrc;
    end

    [~,I] = sort(distanceMatrix);


    %% Sample at random from 10 closest storms
    myStormInd = randsample(I(1:10),1);

    stormInds = trainData.storms{myStormInd}.index;
    orgHS = trainData.X(stormInds);
    stormHS = peakHS/max(orgHS) * orgHS;

    % storm HS drc
    [~,maxind] = max(stormHS);
    peakInd = maxind + stormInds(1) - 1;
    offset = mod(trainData.Xdrc(peakInd) - peakDrc + 180,360) - 180;
    stormHSdrc = mod(trainData.Xdrc(stormInds) - offset,360);

    % storm WS
    orgWS = trainData.Y(stormInds);
    peakWS = beta(1) + beta(2) * peakHS + randn * sig;
%     peakWS = orgWS(maxind) * peakHS/max(orgHS);
    stormWS = peakWS/orgWS(maxind) * orgWS;

    % storm inline WS
    stormWSdrc = mod(trainData.Ydrc(stormInds) - offset,360);
    
    stormsOrgHM{1}{iStorm}.originalMargins = [stormHS,stormWS];
    stormsOrgHM{1}{iStorm}.direction = [stormHSdrc,stormWSdrc];
end

% figure(8); clf;
% col = jet(10);
% for i = 1:10
%     currStormInds = trainData.storms{I(i)}.index;
%     [~,tl] = max(data.X(currStormInds));
%     len = length(currStormInds);
%     subplot(1,2,1);
%     hold on;
%     plot(1-tl:len-tl,data.X(currStormInds),'Color',col(i,:));
%     subplot(1,2,2);
%     hold on;
%     plot(1-tl:len-tl,data.Xdrc(currStormInds),'Color',col(i,:));
% end
% subplot(1,2,1);
% plot(0,peakHS,'k*');
% subplot(1,2,2);
% plot(0,peakDrc,'k*');






end