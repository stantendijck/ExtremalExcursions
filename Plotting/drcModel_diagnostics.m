%% Data
nStorms = length(data.storms);
X = cell(nStorms,1);
Y = cell(nStorms,1);
assX = cell(nStorms,1);
assY = cell(nStorms,1);
diffDrcX = cell(nStorms,1);
diffDrcY = cell(nStorms,1);
for iStorm = 1:nStorms
    I = data.storms{iStorm}.index;
    X{iStorm} = data.X(I);
    Y{iStorm} = data.Y(I);
    
    if length(I) > 1
        assX{iStorm} = X{iStorm}(1:end-1);
        assY{iStorm} = Y{iStorm}(1:end-1);
        diffDrcX{iStorm} = mod(diff(data.Xdrc(I))+180,360)-180;
        diffDrcY{iStorm} = mod(diff(data.Ydrc(I))+180,360)-180;
    end
end

drcX = cell(nStorms,1);
for iStorm = 1:nStorms
    I = data.storms{iStorm}.index;
    drcX{iStorm} = data.Xdrc(I);
    for i = 2:length(drcX{iStorm})
        while drcX{iStorm}(i) - drcX{iStorm}(i-1) > 180
            drcX{iStorm}(i) = drcX{iStorm}(i) - 360;
        end
        while drcX{iStorm}(i) - drcX{iStorm}(i-1) < -180
            drcX{iStorm}(i) = drcX{iStorm}(i) + 360;
        end
    end    
end
drcY = cell(nStorms,1);
for iStorm = 1:nStorms
    I = data.storms{iStorm}.index;
    drcY{iStorm} = data.Ydrc(I);
    if drcY{iStorm}(1) > 260 && drcX{iStorm}(1) < 100
        drcY{iStorm}(1) = drcY{iStorm}(1) - 360;
    elseif drcY{iStorm}(1) < 100 && drcX{iStorm}(1) > 260
        drcY{iStorm}(1) = drcY{iStorm}(1) + 360;
    end
    for i = 2:length(drcY{iStorm})
        while drcY{iStorm}(i) - drcY{iStorm}(i-1) > 180
            drcY{iStorm}(i) = drcY{iStorm}(i) - 360;
        end
        while drcY{iStorm}(i) - drcY{iStorm}(i-1) < -180
            drcY{iStorm}(i) = drcY{iStorm}(i) + 360;
        end
    end    
end

figure(1); clf;
hold on;
iCnt = 0;
for iStorm = 1:nStorms
    len = length(diffDrcX{iStorm});
    plot(iCnt+1:iCnt+len,diffDrcX{iStorm});
    plot(iCnt+1:iCnt+len,diffDrcY{iStorm});
    iCnt = iCnt + len + 3;
end
xlim([0 150]);

drcXmat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),drcX,'UniformOutput',false));
drcYmat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),drcY,'UniformOutput',false));
Xmat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),X,'UniformOutput',false));
Ymat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),Y,'UniformOutput',false));

figure(2); clf;
ax1 = subplot(2,1,1);
plot(1:length(drcXmat),drcXmat,'k-')
hold on;
plot(1:length(drcYmat),drcYmat,'r-')
xlim([0 150])
ax2 = subplot(2,1,2);
plot(1:length(drcXmat),Xmat,'k-')
hold on;
plot(1:length(drcYmat),Ymat,'r-')
xlim([0 150])
linkaxes([ax1,ax2],'x');

%%

% hold on;
iCnt = 0;
avg = cell(nStorms,1);
for iStorm = 1:nStorms
%     len = length(drcX{iStorm});
%     if len >= 10
%         plot(iCnt+1:iCnt+len,drcX{iStorm},'k-');
%         plot(iCnt+1:iCnt+len,drcY{iStorm},'r-');
        
        
        
        windDrc = drcY{iStorm};
        newWindDrc = nan(size(windDrc));
        for i = 1:length(windDrc)
            newWindDrc(i) = mean(windDrc(max(1,-3+i):min(end,1+i)));
        end
%         plot(iCnt+1:iCnt+len,windDrc,'m-');
        
        deltaDrc = mod(newWindDrc - drcX{iStorm} + 180,360) - 180;
        avg{iStorm} = drcX{iStorm} + deltaDrc/2;
        residuals{iStorm} = [windDrc-avg{iStorm},drcX{iStorm}-avg{iStorm}];
%         plot(iCnt+1:iCnt+len,avg{iStorm},'r-');
%         iCnt = iCnt + len + 3;
%     end
end
% xlim([0 150]);

figure(3); clf;
subplot(1,2,1);
plot(cell2mat(assX),cell2mat(diffDrcX),'k.');
xlabel('HS'); ylabel('DeltaDrc-HS');
subplot(1,2,2);
plot(cell2mat(assY),cell2mat(diffDrcY),'k.');
xlabel('WS'); ylabel('DeltaDrc-WS');


figure(20); clf;
hold on;
iCnt = 0;
for iStorm = 1:nStorms
    len = length(drcX{iStorm});
    %         residuals{iStorm} = [windDrc-avg,drcX{iStorm}-avg];
    plot(iCnt+1:iCnt+len,residuals{iStorm}(:,1),'m-');
    plot(iCnt+1:iCnt+len,residuals{iStorm}(:,2),'k-');
    iCnt = iCnt + len + 3;
    %     end
end
xlim([0 150]);

%%
deltaAvg = cell(nStorms,1);
assHS = cell(nStorms,1);
assWS = cell(nStorms,1);
for iStorm = 1:nStorms
    if length(X{iStorm}) >= 2
        deltaAvg{iStorm} = mod(diff(avg{iStorm})+180,360)-180;
        assHS{iStorm} = X{iStorm}(1:end-1);
        assWS{iStorm} = Y{iStorm}(1:end-1);
    end
end
deltaAvg = cell2mat(deltaAvg);
assHS = cell2mat(assHS);
assWS = cell2mat(assWS);

%%
grdHS = linspace(4,15,15);
grdWS = linspace(5,30,15);
grdMidPointsHS = (grdHS(1:end-1)+grdHS(2:end))/2;
grdMidPointsWS = (grdWS(1:end-1)+grdWS(2:end))/2;
std1 = zeros(length(grdMidPointsHS),1);
std2 = zeros(length(grdMidPointsHS),1);
for i = 1:length(grdMidPointsHS)
    std1(i) = std(deltaAvg(assHS>grdHS(i)));
    std2(i) = std(deltaAvg(assWS>grdWS(i)));
end

figure(21); clf;
subplot(1,2,1);
plot(assHS,deltaAvg,'k.');
hold on;
plot(grdMidPointsHS,std1);
% ylim([-30 30]);
subplot(1,2,2);
plot(assWS,deltaAvg,'k.');
hold on;
plot(grdMidPointsWS,std2);
% ylim([-50 50]);



    




%%
% Calculate curvature of HS-Drc
% figure(2); clf;
% hold on;
% iCnt = 0;
curvatureHSDrc = 0; curvatureWindDrc = 0; curvatureNewWindDrc = 0;
for iStorm = 1:nStorms
    curvatureHSDrc = curvatureHSDrc + sum(abs(drcX{iStorm}(1:end-2) - 2*drcX{iStorm}(2:end-1) + drcX{iStorm}(3:end)));
    windDrc = drcY{iStorm};
    newWindDrc = nan(size(windDrc));
    for i = 1:length(windDrc)
        newWindDrc(i) = mean(windDrc(max(1,-3+i):min(end,1+i)));
    end
    curvatureWindDrc = curvatureWindDrc + sum(abs(drcY{iStorm}(1:end-2) - 2*drcY{iStorm}(2:end-1) + drcY{iStorm}(3:end)));
    curvatureNewWindDrc = curvatureNewWindDrc + sum(abs(newWindDrc(1:end-2) - 2*newWindDrc(2:end-1) + newWindDrc(3:end)));
end
fprintf('Curvature Wave direction: %.2f\nCurvature wind direction: %.2f\nCurvature new wind direction: %.2f\n',curvatureHSDrc,curvatureWindDrc,curvatureNewWindDrc)




%%
laggedDataY = cell(nStorms,1);
for iStorm = 1:nStorms
    if length(diffDrcY{iStorm}) > 1
        laggedDataY{iStorm} = [diffDrcY{iStorm}(2:end),diffDrcY{iStorm}(1:end-1)];
    else
        laggedDataY{iStorm} = zeros(0,2);
    end
end
Y = cell2mat(laggedDataY);

laggedDataX = cell(nStorms,1);
for iStorm = 1:nStorms
    if length(diffDrcX{iStorm}) > 1
        laggedDataX{iStorm} = [diffDrcX{iStorm}(2:end),diffDrcX{iStorm}(1:end-1)];
    else
        laggedDataX{iStorm} = zeros(0,2);
    end
end
X = cell2mat(laggedDataX);

% figure(3); clf;
% plot(X(:,1),X(:,2),'k.');
% figure(4); clf;
% plot(X(:,1),Y(:,1),'k.');


%% Simulation
sim = stormsOrg.EVAR{3};
nStorms = length(sim);
% nStorms = 1000;

diffDrcXSim = cell(nStorms,1);
for iStorm = 1:nStorms
    diffDrcXSim{iStorm} = mod(diff(sim{iStorm}.direction(:,2))+180,360)-180;
end

diffDrcYSim = cell(nStorms,1);
for iStorm = 1:nStorms
    diffDrcYSim{iStorm} = mod(diff(sim{iStorm}.direction(:,1))+180,360)-180;
end

drcXsim = cell(nStorms,1);
for iStorm = 1:nStorms
    drcXsim{iStorm} = sim{iStorm}.direction(:,1);
    for i = 2:length(drcXsim{iStorm})
        while drcXsim{iStorm}(i) - drcXsim{iStorm}(i-1) > 180
            drcXsim{iStorm}(i) = drcXsim{iStorm}(i) - 360;
        end
        while drcXsim{iStorm}(i) - drcXsim{iStorm}(i-1) < -180
            drcXsim{iStorm}(i) = drcXsim{iStorm}(i) + 360;
        end
    end    
end
drcYsim = cell(nStorms,1);
for iStorm = 1:nStorms
    drcYsim{iStorm} = sim{iStorm}.direction(:,2);
    if drcYsim{iStorm}(1) > 260 && drcXsim{iStorm}(1) < 100
        drcYsim{iStorm}(1) = drcYsim{iStorm}(1) - 360;
    elseif drcYsim{iStorm}(1) < 100 && drcXsim{iStorm}(1) > 260
        drcYsim{iStorm}(1) = drcYsim{iStorm}(1) + 360;
    end
    for i = 2:length(drcYsim{iStorm})
        while drcYsim{iStorm}(i) - drcYsim{iStorm}(i-1) > 180
            drcYsim{iStorm}(i) = drcYsim{iStorm}(i) - 360;
        end
        while drcYsim{iStorm}(i) - drcYsim{iStorm}(i-1) < -180
            drcYsim{iStorm}(i) = drcYsim{iStorm}(i) + 360;
        end
    end    
end

drcXsimmat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),drcXsim,'UniformOutput',false));
drcYsimmat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),drcYsim,'UniformOutput',false));
Xsimmat = cell2mat(cellfun(@(x)([x.originalMargins(:,1);nan;nan;nan]),sim,'UniformOutput',false));
Ysimmat = cell2mat(cellfun(@(x)([x.originalMargins(:,2);nan;nan;nan]),sim,'UniformOutput',false));

figure(6); clf;
ax3 = subplot(2,1,1);
plot(1:length(drcXsimmat),drcXsimmat,'k-')
hold on;
plot(1:length(drcYsimmat),drcYsimmat,'r-')
xlim([0 150])
ax4 = subplot(2,1,2);
plot(1:length(Xsimmat),Xsimmat,'k-')
hold on;
plot(1:length(Ysimmat),Ysimmat,'r-')
xlim([0 150])
linkaxes([ax3,ax4],'x');




figure(6); clf;
hold on;
iCnt = 0;
for iStorm = 1:nStorms
    len = length(drcXsim{iStorm});
%     if len >= 10
%         fprintf('%d\n',iStorm);
        plot(iCnt+1:iCnt+len,drcXsim{iStorm},'k-');
        plot(iCnt+1:iCnt+len,drcYsim{iStorm},'r-');
        iCnt = iCnt + len + 3;
%     end
end
xlim([0 150]);
title('Wave drc: ARIMA; Wind drc: conditional residuals');



% figure(6); clf;
% iCnt = 0;
% for iStorm = 1:nStorms
%     len = length(diffDrcXSim{iStorm});
%     plot(iCnt+1:iCnt+len,diffDrcXSim{iStorm});
%     hold on;
%     iCnt = iCnt + len + 3;
% end
% xlim([0 150]);

%%
laggedDataY = cell(nStorms,1);
for iStorm = 1:nStorms
    if length(diffDrcYSim{iStorm}) > 1
        laggedDataY{iStorm} = [diffDrcYSim{iStorm}(2:end),diffDrcYSim{iStorm}(1:end-1)];
    else
        laggedDataY{iStorm} = zeros(0,2);
    end
end
Ysim = cell2mat(laggedDataY);

laggedDataX = cell(nStorms,1);
for iStorm = 1:nStorms
    if length(diffDrcXSim{iStorm}) > 1
        laggedDataX{iStorm} = [diffDrcXSim{iStorm}(2:end),diffDrcXSim{iStorm}(1:end-1)];
    else
        laggedDataX{iStorm} = zeros(0,2);
    end
end
Xsim = cell2mat(laggedDataX);

figure(7); clf;
plot(Xsim(:,1),Xsim(:,2),'k.');
figure(8); clf;
plot(Xsim(:,1),Ysim(:,1),'k.');


%% Linear regression

XsimMat = [ones(size(Xsim(:,1))),Xsim(:,1)];
betaSim =  (XsimMat' * XsimMat)^(-1) * XsimMat' * Xsim(:,2);

XdataMat = [ones(size(X(:,1))),X(:,1)];
betaData =  (XdataMat' * XdataMat)^(-1) * XdataMat' * X(:,2);

%%
binsX = cell(nStorms,1); binsY = cell(nStorms,1);
for iStorm = 1:nStorms
    binsX{iStorm} = BinAlc(data.Xdrc(data.storms{iStorm}.index), MM.X{2}.DrcEdg);
    binsY{iStorm} = BinAlc(data.Ydrc(data.storms{iStorm}.index), MM.Y{2}.DrcEdg);
end

figure(10); clf;
hold on;
iCnt = 0;
cntrY = 0; cntrX = 0;cntrXY = 0;
for iStorm = 1:nStorms 
    len = length(binsX{iStorm});
    plot(iCnt+1:iCnt+len,binsX{iStorm},'k-');
    plot(iCnt+1:iCnt+len,binsY{iStorm},'r-');
    iCnt = iCnt + len + 3;
    if all(binsX{iStorm}==binsX{iStorm}(1))
        cntrX = cntrX + 1;
    end
    if all(binsY{iStorm}==binsY{iStorm}(1))
        cntrY = cntrY + 1;
    end
    if all(binsX{iStorm}==binsX{iStorm}(1)) && all(binsY{iStorm}==binsY{iStorm}(1))
        cntrXY = cntrXY + 1;
    end
end
fprintf('Percentage of storms for which wave direction stays in same bin: %.2f\n',cntrX/nStorms*100);
fprintf('Percentage of storms for which wind direction stays in same bin: %.2f\n',cntrY/nStorms*100);
fprintf('Percentage of storms for direction stays in the same bin: %.2f\n',cntrXY/nStorms*100);

matX = zeros(8,8);
for iStorm = 1:nStorms
    for iObs = 2:length(binsX{iStorm})
        iX = binsX{iStorm}(iObs-1);
        iY = binsX{iStorm}(iObs);
        matX(iX,iY) = matX(iX,iY) + 1;
    end
end
matX./sum(matX,2)

matXYX = zeros(8,3,8);
for iStorm = 1:nStorms
    for iObs = 2:length(binsX{iStorm})
        iX = binsX{iStorm}(iObs-1);
        iY = binsY{iStorm}(iObs-1);
        iZ = binsX{iStorm}(iObs);
        matXYX(iX,iY,iZ) = matXYX(iX,iY,iZ) + 1;
    end
end
matXYX./sum(matXYX,3)

matXYY = zeros(8,3,3);
for iStorm = 1:nStorms
    for iObs = 2:length(binsX{iStorm})
        iX = binsX{iStorm}(iObs-1);
        iY = binsY{iStorm}(iObs-1);
        iZ = binsY{iStorm}(iObs);
        matXYY(iX,iY,iZ) = matXYY(iX,iY,iZ) + 1;
    end
end
matXYY./sum(matXYY,3)

matXYY = zeros(8,3,8,3);
for iStorm = 1:nStorms
    for iObs = 2:length(binsX{iStorm})
        iX = binsX{iStorm}(iObs-1);
        iY = binsY{iStorm}(iObs-1);
        iZ = binsX{iStorm}(iObs);
        iW = binsY{iStorm}(iObs);
        matXYY(iX,iY,iZ,iW) = matXYY(iX,iY,iZ,iW) + 1;
    end
end
matXYY./sum(matXYY,3)

matY = zeros(3,3);
for iStorm = 1:nStorms
    for iObs = 2:length(binsY{iStorm})
        iX = binsY{iStorm}(iObs-1);
        iY = binsY{iStorm}(iObs);
        matY(iX,iY) = matY(iX,iY) + 1;
    end
end
matY./sum(matY,2)

%%

