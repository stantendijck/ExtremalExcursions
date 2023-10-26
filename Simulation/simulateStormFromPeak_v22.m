function storm = simulateStormFromPeak_v22(mdlPeak, mdlAroundPeak, mdlFalling, mdlRising, p, thr, method, rsdModel, stormStart)
% Simulate storm (x>thr + ass.var) starting from the storm peak using
% VAR(p)        if method = 'VAR'
% EVAR(p)       if method = 'EVAR'
% MMEM(p)       if method = 'MMEM'
% MEM(p)        if method = 'MEM'
%
% INPUT
% mdlPeak           struct(): storm peak model
% mdlAroundPeak     struct(): model around storm peak
% mdlFalling        struct(): model for falling part of a storm
% mdlRising         struct(): model for rising part of a storm
% p                 int >= 1: model order
% thr               float:    storm threshold
% method            string:   model identifier for falling and rising parts
% rsdModel          string:   "none" or "copula"
% stormStart        array:    initialisation of stormpeak
%
% OUTPUT
% storm             array:    simulated storm

%% Initialisation
if nargin < 8
    rsdModel = 'none';
end
switch rsdModel
    case 'none'
        flagCopula = false;
    otherwise
        flagCopula = true;
end
switch method
    case 'VAR'
        attr = 'varn';
    case 'EVAR'
        attr = 'typen';
end
if nargin < 9
    stormStart = [];
end
varNms = {'X','Y'};

% PostPeak initialisation
if flagCopula
    switch method
        case {'VAR','EVAR'}
            tDistPostPeak = struct();
            tDistPostPeak.rhoHAT = mdlFalling.rsd.(attr){p}.rhoHAT;
            tDistPostPeak.nuHAT = mdlFalling.rsd.(attr){p}.nuHAT;
            tDistPostPeak.icdf = mdlFalling.rsd.(attr){p}.margHAT;
        case 'MMEM'
            nSizePostPeak = floor(size(mdlFalling.rsd.type1.rhoHAT,1)/2);
            tDistPostPeak.rhoHAT = mdlFalling.rsd.type1.rhoHAT([1:p,nSizePostPeak+1:nSizePostPeak+p],[1:p,nSizePostPeak+1:nSizePostPeak+p]);
            tDistPostPeak.nuHAT = mdlFalling.rsd.type1.nuHAT;
            tDistPostPeak.icdf = mdlFalling.rsd.type1.margHAT([1:p,nSizePostPeak+1:nSizePostPeak+p]);
    end
    switch method
        case {'VAR','EVAR'}
            tDistPrePeak = struct();
            tDistPrePeak.rhoHAT = mdlRising.rsd.(attr){p}.rhoHAT;
            tDistPrePeak.nuHAT = mdlRising.rsd.(attr){p}.nuHAT;
            tDistPrePeak.icdf = mdlRising.rsd.(attr){p}.margHAT;
        case 'MMEM'
            nSizePrePeak = floor(size(mdlRising.rsd.type1.rhoHAT,1)/2);
            tDistPrePeak.rhoHAT = mdlRising.rsd.type1.rhoHAT([1:p,nSizePrePeak+1:nSizePrePeak+p],[1:p,nSizePrePeak+1:nSizePrePeak+p]);
            tDistPrePeak.nuHAT = mdlRising.rsd.type1.nuHAT;
            tDistPrePeak.icdf = mdlRising.rsd.type1.margHAT([1:p,nSizePrePeak+1:nSizePrePeak+p]);
    end
end
switch method
    case 'MMEM'
        alphasPostPeak = [cellfun(@(x)(x.newMLE(1)), mdlFalling.type1.X(1:p));...
            cellfun(@(x)(x.newMLE(1)), mdlFalling.type1.Y(1:p))];
        
        betasPostPeak = [cellfun(@(x)(x.newMLE(2)), mdlFalling.type1.X(1:p));...
            cellfun(@(x)(x.newMLE(2)), mdlFalling.type1.Y(1:p))];
        
        allRsdPostPeak = [cell2mat(cellfun(@(x)(x.MLERsd),mdlFalling.type1.X(1:p),'UniformOutput',false)'),...
            cell2mat(cellfun(@(x)(x.MLERsd),mdlFalling.type1.Y(1:p),'UniformOutput',false)')];
        
        alphasPrePeak = [cellfun(@(x)(x.newMLE(1)), mdlRising.type1.X(1:p));...
            cellfun(@(x)(x.newMLE(1)), mdlRising.type1.Y(1:p))];
        
        betasPrePeak = [cellfun(@(x)(x.newMLE(2)), mdlRising.type1.X(1:p));...
            cellfun(@(x)(x.newMLE(2)), mdlRising.type1.Y(1:p))];
        
        allRsdPrePeak = [cell2mat(cellfun(@(x)(x.MLERsd),mdlRising.type1.X(1:p),'UniformOutput',false)'),...
            cell2mat(cellfun(@(x)(x.MLERsd),mdlRising.type1.Y(1:p),'UniformOutput',false)')];
end

%% Repeat until a valid storm is simulated
flagValidStorm = false;
while ~flagValidStorm
    flagValidStorm = true;
    
    %% Step 1: simulate storm peak
    if any(stormStart); u = stormStart(1,1); else; u = rand; end
    
    if u > mdlPeak.X.q
        stormpeak = GPD_iCDF((u - mdlPeak.X.q)/(1-mdlPeak.X.q), mdlPeak.X.MLE(3), mdlPeak.X.MLE(1), mdlPeak.X.MLE(2));
    else
%         stormpeak = randsample(mdlPeak.X.data(mdlPeak.X.data < mdlPeak.X.MLE(1)), 1);
        stormpeak = quantile(mdlPeak.X.data(mdlPeak.X.data < mdlPeak.X.MLE(1)), u/mdlPeak.X.q);
    end

    %% Step 2: simulate around storm peak

    % Residuals
    if p > 1
        if flagCopula
            nextRsdU = tcdf(mvtrnd(mdlAroundPeak{p}.rsd.rhoHAT, mdlAroundPeak{p}.rsd.nuHAT, 1), mdlAroundPeak{p}.rsd.nuHAT);
            nextRsd = nan(size(nextRsdU));
            for j = 1:length(nextRsdU)
                nextRsd(j) = myinterp1(mdlAroundPeak{p}.rsd.marghat{j}(:,1), mdlAroundPeak{p}.rsd.marghat{j}(:,2),nextRsdU(j));
            end
        else
            nextRsd = simulateRsd([],mdlAroundPeak{p}.rsd.data)';
        end
    else
        nextRsdU = rand;
        nextRsd = myinterp1(mdlAroundPeak{p}.rsd.marghat{1}(:,1), mdlAroundPeak{p}.rsd.marghat{1}(:,2),nextRsdU);
    end
    
    % Associated values
    nextVal = nan(size(nextRsd));
    for j = 1:length(nextRsd)
%         nextRsd(j) = myinterp1(mdlAroundPeak{p}.rsd.marghat{j}(:,1), mdlAroundPeak{p}.rsd.marghat{j}(:,2),nextRsdU(j));
        nextVal(j) = mdlAroundPeak{p}.phat{j}(1) * stormpeak + stormpeak .^ mdlAroundPeak{p}.phat{j}(2) * nextRsd(j);
    end
    storm = reshape([nextVal(1:p-1),stormpeak,nextVal(p:end)],2*p-1,2);
    
    if any(any(isnan(storm)))
        flagValidStorm = false;
    end
    
    % Check if simulated storm is valid and if the storm has ended
    flagFallen = false; flagRisen = false;
    for iObs = p+1:size(storm,1)
        if storm(iObs,1) > stormpeak
            flagValidStorm = false;
            break
        elseif storm(iObs,1) < thr
            flagFallen = true;
            storm = storm(1:iObs-1,:);
            break
        end
    end
    for iObs = p-1:-1:1
        if storm(iObs,1) > stormpeak
            flagValidStorm = false;
            break
        elseif storm(iObs,1) < thr
            flagRisen = true;
            storm = storm(iObs+1:end,:);
            break
        end
    end
    
    %% Step 3: Simulate postpeak trajectory
    
    while flagValidStorm && ~flagFallen    
        switch method
            case {'VAR','EVAR'}
                % Simulate residuals using t copula or multivariate ksdens
                if flagCopula
                    nextRsd = simulateRsd([], tDistPostPeak);
                else
                    nextRsd = simulateRsd([], [mdlFalling.(attr){p}.X{1}.MLERsd, mdlFalling.(attr){p}.Y{1}.MLERsd]);
                end

                % Associated values
                nextVal = nan(length(varNms),1);
                for iVarNme = 1:length(varNms)
                    alp = mdlFalling.(attr){p}.(varNms{iVarNme}){1}.newMLE(1:end-1);
                    bet = mdlFalling.(attr){p}.(varNms{iVarNme}){1}.newMLE(end);
                    dat = reshape(storm(end-p+1:end,:),[],1);
                    nextVal(iVarNme) = dat' * alp + dat(1)^bet * nextRsd(iVarNme);
                end
            case 'MMEM'
                condValue = storm(end-p+1, 1);
                
                prevRsd = nan(2*(p-1), 1);
                for i = 1:p-1
                    prevRsd(i) = (storm(end-p+1+i, 1) - alphasPostPeak(i) * condValue) ./ condValue.^betasPostPeak(i);
                    prevRsd(p-1+i) = (storm(end-p+1+i-1, 2) - alphasPostPeak(p+i) * condValue) ./ condValue.^betasPostPeak(p+i);
                end
                
                prevInds = [1:p-1, p+1:2*p-1]; newInds = [p, 2*p];
                
                if flagCopula
                    ttDist.rhoHAT = tDistPostPeak.rhoHAT([prevInds, newInds],[prevInds, newInds]);
                    ttDist.nuHAT = tDistPostPeak.nuHAT;
                    ttDist.icdf = tDistPostPeak.icdf([prevInds, newInds]);
                    nextRsd = simulateRsd(prevRsd, ttDist);
                else
                    nextRsd = simulateRsd(prevRsd, allRsdPostPeak(:,[prevInds, newInds]));
                end
                
                nextVal = alphasPostPeak(newInds) * condValue + condValue .^ betasPostPeak(newInds) .* nextRsd(end-length(newInds)+1:end);                
        end
        
        % Checks
        if nextVal(1) > stormpeak
            flagValidStorm = false;
        elseif nextVal(1) < thr
            flagFallen = true;
        else
            storm(end+1,:) = nextVal;
        end
    end
    
    %% Step 4: Simulate prepeak trajectory
    storm = storm(end:-1:1,:);
    while flagValidStorm && ~flagRisen
        switch method
            case {'VAR','EVAR'}
                % Simulate residual using t copula or multivariate ksdens
                if flagCopula
                    tDistPostPeak = struct();
                    tDistPostPeak.rhoHAT = mdlRising.rsd.(attr){p}.rhoHAT;
                    tDistPostPeak.nuHAT = mdlRising.rsd.(attr){p}.nuHAT;
                    tDistPostPeak.icdf = mdlRising.rsd.(attr){p}.margHAT;
                    nextRsd = simulateRsd([], tDistPostPeak);
                else
                    nextRsd = simulateRsd([], [mdlRising.(attr){p}.X{1}.MLERsd,mdlRising.(attr){p}.Y{1}.MLERsd]);
                end
                
                % Associated values
                nextVal = nan(length(varNms),1);
                for iVarNme = 1:length(varNms)
                    alp = mdlRising.(attr){p}.(varNms{iVarNme}){1}.newMLE(1:end-1);
                    bet = mdlRising.(attr){p}.(varNms{iVarNme}){1}.newMLE(end);
                    dat = reshape(storm(end-p+1:end,:),[],1);
                    nextVal(iVarNme) = dat' * alp + dat(1)^bet * nextRsd(iVarNme);
                end
            case 'MMEM'
                condValue = storm(end-p+1, 1);
                
                prevRsd = nan(2*(p-1), 1);
                for i = 1:p-1
                    prevRsd(i) = (storm(end-p+1+i, 1) - alphasPrePeak(i) * condValue) ./ condValue.^betasPrePeak(i);
                    prevRsd(p-1+i) = (storm(end-p+1+i-1, 2) - alphasPrePeak(p+i) * condValue) ./ condValue.^betasPrePeak(p+i);
                end
                
                prevInds = [1:p-1, p+1:2*p-1]; newInds = [p, 2*p];
                
                if flagCopula
                    ttDist.rhoHAT = tDistPrePeak.rhoHAT([prevInds, newInds],[prevInds, newInds]);
                    ttDist.nuHAT = tDistPrePeak.nuHAT;
                    ttDist.icdf = tDistPrePeak.icdf([prevInds, newInds]);
                    nextRsd = simulateRsd(prevRsd, ttDist);
                else
                    nextRsd = simulateRsd(prevRsd, allRsdPrePeak(:,[prevInds, newInds]));
                end
                
                nextVal = alphasPrePeak(newInds) * condValue + condValue .^ betasPrePeak(newInds) .* nextRsd(end-length(newInds)+1:end);                
        end
        
        if nextVal(1) > stormpeak
            flagValidStorm = false;
        elseif nextVal(1) < thr
            flagRisen = true;
        else
            storm(end+1,:) = nextVal;
        end
    end

    storm = storm(end:-1:1,:);
end
    



end