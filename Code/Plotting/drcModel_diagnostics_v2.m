function drcModel_diagnostics_v2(sim, data)
%%
% sim = stormsOrg.EVAR{3};

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


drcXmat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),drcX,'UniformOutput',false));
drcYmat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),drcY,'UniformOutput',false));
Xmat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),X,'UniformOutput',false));
Ymat = cell2mat(cellfun(@(x)([x;nan;nan;nan]),Y,'UniformOutput',false));

figure(8); clf;
ax1 = subplot(2,1,1);
plot(1:length(drcXmat),drcXmat,'k-')
hold on;
plot(1:length(drcYmat),drcYmat,'r-')
xlim([0 150])
ylabel('Direction');
xlabel('t');
title('Data (red = wind, black = wave)');
ax2 = subplot(2,1,2);
plot(1:length(drcXmat),Xmat,'k-')
hold on;
plot(1:length(drcYmat),Ymat,'r-')
xlim([0 150])
ylabel('Size');
xlabel('t');
linkaxes([ax1,ax2],'x');

%% Simulation
nStorms = length(sim);
% nStorms = 1000;

Xsim = cell(nStorms,1);
Ysim = cell(nStorms,1);
diffDrcXSim = cell(nStorms,1);
diffDrcYSim = cell(nStorms,1);
for iStorm = 1:nStorms
    Xsim{iStorm} = sim{iStorm}.originalMargins(:,1);
    Ysim{iStorm} = sim{iStorm}.originalMargins(:,2);
    diffDrcXSim{iStorm} = mod(diff(sim{iStorm}.direction(:,2))+180,360)-180;
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

%%
figure(9); clf;
ax3 = subplot(2,1,1);
plot(1:length(drcXsimmat),drcXsimmat,'k-')
hold on;
plot(1:length(drcYsimmat),drcYsimmat,'r-')
xlim([0 150])
% xlim([31700 31850])
ylabel('Direction');
xlabel('t');
title('Simulation (red = wind, black = wave)');
ax4 = subplot(2,1,2);
plot(1:length(Xsimmat),Xsimmat,'k-')
hold on;
plot(1:length(Ysimmat),Ysimmat,'r-')
xlim([0 150])
% xlim([31700 31850])
ylabel('Size');
xlabel('t');
linkaxes([ax3,ax4],'x');

%% check: rate of peak directions

peakDrcInd = cellfun(@(x)(find(x==max(x))),X,'UniformOutput',false);
peakDrcHS = cellfun(@(x,y)(x(y)),drcX,peakDrcInd);

peakDrcIndsim = cellfun(@(x)(find(x==max(x))),Xsim,'UniformOutput',false);
peakDrcsimHS = cellfun(@(x,y)(x(y)),drcXsim,peakDrcIndsim);

peakDrcInd = cellfun(@(x)(find(x==max(x))),Y,'UniformOutput',false);
peakDrcWS = cellfun(@(x,y)(x(y)),drcY,peakDrcInd);

peakDrcIndsim = cellfun(@(x)(find(x==max(x))),Ysim,'UniformOutput',false);
peakDrcsimWS = cellfun(@(x,y)(x(y)),drcYsim,peakDrcIndsim);


figure(10); clf;
subplot(2,2,1);
histogram(mod(peakDrcsimHS,360),'NumBins',36,'Normalization','pdf','FaceColor','black');
hold on;
histogram(mod(peakDrcHS,360),'NumBins',36,'Normalization','pdf','FaceColor','red');
title('Peak associated direction H_S');
xlabel('Direction');

subplot(2,2,2);
plot(mod(peakDrcsimHS,360),mod(peakDrcsimWS,360),'k.');
hold on;
plot(mod(peakDrcHS,360),mod(peakDrcWS,360),'r.');
xlabel('Drc-HS'); ylabel('Drc-WS');

subplot(2,2,3);
plot(mod(peakDrcsimWS,360),mod(peakDrcsimHS,360),'k.');
hold on;
plot(mod(peakDrcWS,360),mod(peakDrcHS,360),'r.');
ylabel('Drc-HS'); xlabel('Drc-WS');

subplot(2,2,4); 
histogram(mod(peakDrcsimWS,360),'NumBins',36,'Normalization','pdf','FaceColor','black');
hold on;
histogram(mod(peakDrcWS,360),'NumBins',36,'Normalization','pdf','FaceColor','red');
xlabel('Direction');
title('Peak associated direction W_S');

%% check: hs - ws

dsim = mod(drcXsimmat - drcYsimmat + 180,360) - 180;
ddata = mod(drcXmat - drcYmat + 180,360) - 180;

figure(11); clf;
% subplot(1,2,1);
subplot(1,2,1);
plot(Xsimmat,dsim,'k.');
hold on;
plot(Xmat,ddata,'r.');
title('(H_S, Drc(H_S) - Drc(W_S))');
xlabel('H_S');
ylabel('Drc(H_S) - Drc(W_S)');
legend({'sim','data'});


subplot(1,2,2);
plot(Ysimmat,dsim,'k.');
hold on;
plot(Ymat,ddata,'r.');
title('(W_S, Drc(H_S) - Drc(W_S))');
xlabel('W_S');
ylabel('Drc(H_S) - Drc(W_S)');
legend({'sim','data'});

%% Distributions of Hs, Ws and (Hs,Ws)
figure(12); clf;
subplot(2,2,1);
h = histogram(Xmat,'Normalization','pdf','FaceColor','black');
hold on;
histogram(Xsimmat,'Normalization','pdf','BinWidth',h.BinWidth,'FaceColor','red');
xlabel('H_S');
legend({'data','sim'});

subplot(2,2,2);
plot(Xsimmat,Ysimmat,'k.');
hold on;
plot(Xmat,Ymat,'r.');
xlabel('H_S');
ylabel('W_S');
legend({'sim','data'});

subplot(2,2,3);
plot(Ysimmat,Xsimmat,'k.');
hold on;
plot(Ymat,Xmat,'r.');
xlabel('H_S');
ylabel('W_S');
legend({'sim','data'});

subplot(2,2,4);
h = histogram(Ymat,'Normalization','pdf','FaceColor','black');
hold on;
histogram(Ysimmat,'Normalization','pdf','BinWidth',h.BinWidth,'FaceColor','red');
xlabel('W_S');
legend({'data','sim'});

%% Lagged distributions

figure(21); clf;

subplot(2,3,1);
plot(Xsimmat(1:end-1),Xsimmat(2:end),'k.')
hold on;
plot(Xmat(1:end-1),Xmat(2:end),'r.')
xlabel('H_S(t)'); ylabel('H_S(t+1)');

subplot(2,3,2);
plot(Xsimmat(1:end-2),Xsimmat(3:end),'k.')
hold on;
plot(Xmat(1:end-2),Xmat(3:end),'r.')
xlabel('H_S(t)'); ylabel('H_S(t+2)');

subplot(2,3,3);
plot(Xsimmat(1:end-3),Xsimmat(4:end),'k.')
hold on;
plot(Xmat(1:end-3),Xmat(4:end),'r.')
xlabel('H_S(t)'); ylabel('H_S(t+3)');

subplot(2,3,4);
plot(Ysimmat(1:end-1),Ysimmat(2:end),'k.')
hold on;
plot(Ymat(1:end-1),Ymat(2:end),'r.')
xlabel('W_S(t)'); ylabel('W_S(t+1)');

subplot(2,3,5);
plot(Ysimmat(1:end-2),Ysimmat(3:end),'k.')
hold on;
plot(Ymat(1:end-2),Ymat(3:end),'r.')
xlabel('W_S(t)'); ylabel('W_S(t+2)');

subplot(2,3,6);
plot(Ysimmat(1:end-3),Ysimmat(4:end),'k.')
hold on;
plot(Ymat(1:end-3),Ymat(4:end),'r.')
xlabel('W_S(t)'); ylabel('W_S(t+3)');



