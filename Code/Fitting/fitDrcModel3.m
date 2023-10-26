function mdldrc = fitDrcModel3(data)

mdldrc = struct();

% only storm indices!
nStorms = length(data.storms);
X = cell(nStorms,1); Y = cell(nStorms,1);
assX1 = cell(nStorms,1); assY1 = cell(nStorms,1);
drcX = cell(nStorms,1); drcY = cell(nStorms,1);
diffDrcX = cell(nStorms,1); diffDrcY = cell(nStorms,1);
for iStorm = 1:nStorms
    I = data.storms{iStorm}.index;
    X{iStorm} = data.X(I);
    Y{iStorm} = data.Y(I);
    
    drcX{iStorm} = data.Xdrc(I);
    drcY{iStorm} = data.Ydrc(I);
    
    if length(X{iStorm}) > 1
        assX1{iStorm} = X{iStorm}(1:end-1);
        assY1{iStorm} = Y{iStorm}(1:end-1);
    end
    diffDrcX{iStorm} = mod(diff(data.Xdrc(I))+180,360)-180;
    diffDrcY{iStorm} = mod(diff(data.Ydrc(I))+180,360)-180;
end


gaussianDDX = norminv(epit(cell2mat(diffDrcX)));
gaussianDDY = norminv(epit(cell2mat(diffDrcY)));
vecOfLengths = cellfun(@length,diffDrcX);

gaussianDDX_cell = cell(nStorms,1);
gaussianDDY_cell = cell(nStorms,1);

assX = cell(nStorms,1); assY = cell(nStorms,1);
stackedDataX = cell(nStorms,1); stackedDataY = cell(nStorms,1);
cnt = 1;
for iStorm = 1:nStorms
    if vecOfLengths(iStorm) > 0
        gaussianDDX_cell{iStorm} = gaussianDDX(cnt:cnt+vecOfLengths(iStorm)-1);
        gaussianDDY_cell{iStorm} = gaussianDDY(cnt:cnt+vecOfLengths(iStorm)-1);
        cnt = cnt + vecOfLengths(iStorm);
    end

    if vecOfLengths(iStorm) >= 2
        assX{iStorm} = X{iStorm}(2:end-1);
        assY{iStorm} = Y{iStorm}(2:end-1);
        
        stackedDataX{iStorm} = [gaussianDDX_cell{iStorm}(2:end),gaussianDDX_cell{iStorm}(1:end-1)];
        stackedDataY{iStorm} = [gaussianDDY_cell{iStorm}(2:end),gaussianDDY_cell{iStorm}(1:end-1)];
        
%         stackedDataOrgX{iStorm} = [gaussianDDX_cell{iStorm}(2:end),gaussianDDX_cell{iStorm}(1:end-1)];
%         stackedDataOrgY{iStorm} = [gaussianDDY_cell{iStorm}(2:end),gaussianDDY_cell{iStorm}(1:end-1)];
    end
end

assX = cell2mat(assX);
assY = cell2mat(assY);
stackedDataX = cell2mat(stackedDataX);
stackedDataY = cell2mat(stackedDataY);


%%
opts = optimset('MaxFunEvals',10000,'MaxIter',10000);
p0 = [1,1,1,1];
like = @(data, p)(drchslike3(data{1}, data{2}, p(1), p(2:4)));
phatX = fminsearch(@(p)like({stackedDataX, assX}, p), p0, opts);


MCMC = adptMCMC({stackedDataX, assX},like,phatX',1e3, 0.01);
MCMC = adptMCMC({stackedDataX, assX},like,MCMC.MLE,1e4, 0.01);
% MCMC = adptMCMC({stackedDataX, assX},like,MCMC.MLE,1e4, 0.001);
% MCMC = adptMCMC({stackedDataX, assX},like,MCMC.MLE,1e4, 0.001);
% pltMCMC(MCMC);

mdldrc.X.phat = MCMC.MLE;


%%

rsd = cell(nStorms,1);
for iStorm = 1:nStorms
    if length(drcX{iStorm})>1
        rsd{iStorm} = [drcX{iStorm}(1:end-1)-drcY{iStorm}(1:end-1),drcX{iStorm}(2:end) - drcY{iStorm}(2:end)];
    end
end

rsd = mod(cell2mat(rsd)+180,360)-180;
mdldrc.Y.rsd = rsd;


mdldrc.X.deltaX = stackedDataX;
mdldrc.Y.deltaY = stackedDataY;
mdldrc.X.assX = assX;
mdldrc.Y.assY = assY;

mX = cell2mat(diffDrcX);
mY = cell2mat(diffDrcY);
% mdldrc.theta.X = sort(mX(cell2mat(assX1)>8));
% mdldrc.theta.Y = sort(mY(cell2mat(assX1)>8));
mdldrc.theta.X = sort(mX);
mdldrc.theta.Y = sort(mY);


end