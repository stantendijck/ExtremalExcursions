%% Illustration peak period, post-peak and pre-peak
rng(1234);

t = -4:5;
X = [1.2,2.6,2.4,3,3.3,3.1,2.8,1.6,2,1.4];
Y = [2,1.7,2.5,2.8,2.7,1.8,1.75,1.4,0.8,0.3];

t2 = -2:3;
X2 = [1.4,1.2,2,1.6,1.8,1.1];
Y2 = [1.9,2,1.4,1.2,1.6,0.9];

figure(10); clf;
tl = tiledlayout(1,2,'tilespacing','compact','padding','compact');
nexttile;
plot(t,X,'k-o');
hold on;
plot(t,Y,'r-o');
xlim([-5 7]);
plot(xlim,[1.5,1.5],'b--');
grid on;
% text('peak');

an = annotation('doublearrow',[0.13,0.235],[0.3,0.3]);
c = an.Color;
an.Color = 'black';
text(-2.3,1,'Pre peak');

an = annotation('doublearrow',[0.235,0.372],[0.3,0.3]);
c = an.Color;
an.Color = 'black';
text(1.3,1,'Post peak');

an = annotation('doublearrow',[0.17,0.3],[0.25,0.25]);
c = an.Color;
an.Color = 'black';
text(-0.5,0.78,'Peak');
xlabel('t');
ylabel('size');


nexttile;
plot(t2,X2,'k-o');
hold on;
plot(t2,Y2,'r-o');
xlim([-3 4]);
plot(xlim,[1.5,1.5],'b--');
grid on;
ylim([0 2.2]);



an = annotation('doublearrow',[0.715,0.725],[0.3,0.3]);
c = an.Color;
an.Color = 'black';
text(-0.5,0.63,'Pre peak');

an = annotation('doublearrow',[0.725,0.835],[0.3,0.3]);
c = an.Color;
an.Color = 'black';
text(0.6,0.63,'Post peak');

an = annotation('doublearrow',[0.605,0.835],[0.25,0.25]);
c = an.Color;
an.Color = 'black';
text(-0.2,0.5,'Peak');

xlabel('t');
ylabel('size');

legend({'Y_1','Y_2'})

savePics('figurespaper/illustration.pdf',saveOn,'beamer');

%% Matrix plot: Hs, Ws
I = cell2mat(cellfun(@(x)(x.index'),data.storms,'UniformOutput',false));

figure(11); clf;

t = tiledlayout(3,4,'TileSpacing','none','Padding','compact');


%% HS
colValue = 0;
for i = 1:4
    nexttile;
    
    if i == 1
        histogram(data.lap.X,'BinWidth',0.8,'Normalization','pdf');
        set(gca,'view',[90 -90],'ytick',0:0.1:0.4,'yticklabel',[]);
        xlabel('H^L_{S,i}');
        ylabel('Histogram');
    else
        plot(data.lap.X(i:end),data.lap.X(1:end-i+1),'.','Color',colValue*ones(1,3));
        xlabel(sprintf('H^L_{S,i+%d}',i));
        set(gca,'xticklabel',[],'yticklabel',[]);
        ylim([-14 14]);
    end
    %         set(gca,'xticklabel',[]);
    %         ylabel('H^L_{S,t}');
    xlim([-14 14]); %
    grid on;
end

%% HS,WS
for i = 1:4
    nexttile;
    plot(data.lap.Y(i:end),data.lap.X(1:end-i+1),'.','Color',colValue*ones(1,3));
    
    if i == 1
        xlabel('W^L_{s,i}');
        ylabel('H^L_{S,i}');
        set(gca,'xticklabel',[]);
    else
        xlabel(sprintf('W^L_{s,i+%d}',i));
        set(gca,'xticklabel',[],'yticklabel',[]);
    end
    xlim([-14 14]); ylim([-14 14]);
    grid on;
end

%% WS

for i = 1:4
    nexttile;
    if i == 1
        histogram(data.lap.Y,'BinWidth',0.8,'Normalization','pdf');
        set(gca,'view',[90 -90],'ytick',0:0.1:0.4,'yticklabel',[]);
        xlabel('W^L_{s,i}');
        ylabel('Histogram');
    else
        plot(data.lap.Y(i:end),data.lap.Y(1:end-i+1),'.','Color',colValue*ones(1,3));
        xlabel(sprintf('W^L_{s,i+%d}',i));
        set(gca,'xticklabel',[],'yticklabel',[]);
        ylim([-14 14]);
    end
    xlim([-14 14]);
    grid on;
end
savePics('figurespaper/matrixplot.pdf',saveOn,'beamer');

%% Slice of time-series

[~,I] = sort(cellfun(@(x)(-data.X(x.peakInd)),testData.storms));

start= 1;
newDrcX = nan(length(data.Xdrc),1);
newDrcY = nan(length(data.Xdrc),1);
changeX = nan(length(data.Xdrc)-1,1);
changeY = nan(length(data.Xdrc)-1,1);
newDrcX(1) = data.Xdrc(start);
newDrcY(1) = data.Ydrc(start);
for i = 1:length(data.Xdrc)-1
    changeX(i) = mod(data.Xdrc(start+i) - data.Xdrc(start+i-1) + 180,360) - 180;
    newDrcX(i+1) = newDrcX(i) + changeX(i);
    changeY(i) = mod(data.Ydrc(start+i) - data.Ydrc(start+i-1) + 180,360) - 180;
    newDrcY(i+1) = newDrcY(i) + changeY(i);
    while abs(newDrcY(i+1) - newDrcX(i+1)) > 180
        sign = (newDrcY(i+1) - newDrcX(i+1)) ./ abs(newDrcY(i+1) - newDrcX(i+1));
        newDrcY(i+1) = newDrcY(i+1) - sign*360;
    end
end


wt = datetime(data.t,'ConvertFrom','datenum');

figure(12); clf;
t = tiledlayout(3,4,'TileSpacing','none','Padding','compact');

inds = [1, round(length(I)*0.05), round(length(I)*0.1), round(length(I)*0.15)];
for i = 1:4
    myInd = I(inds(i));
    pk = wt(testData.storms{myInd}.peakInd);
    
    
    nexttile(i);
    plot(wt,data.X,'k-');
    hold on;
    plot(wt,data.Y,'r-');
    xlim([pk-30/8 pk+30/8]);
    if i == 4, legend({'H_{S,i}','W_{s,i}'}); end
    if i == 1, ylabel('size'); else; set(gca,'yticklabel',[]); end
    set(gca,'xticklabel',[]);
    grid on;
    ylim([0 30]);
    
    
    
    nexttile(4+i);
    plot(wt,data.lap.X,'k-');
    hold on;
    plot(wt,data.lap.Y,'r-');
    xlim([pk-30/8 pk+30/8]);
    if i == 4, legend({'H_{S,i}^L','W_{s,i}^L'}); end
    if i == 1, ylabel('size'); else; set(gca,'yticklabel',[]); end
    set(gca,'xticklabel',[]);
    grid on;
    ylim([-5 11]);
    
    
    newDrcY = newDrcY - 360*floor(newDrcX(testData.storms{myInd}.peakInd)/360);
    newDrcX = newDrcX - 360*floor(newDrcX(testData.storms{myInd}.peakInd)/360);
    
    nexttile(8+i);
    xl = start+(0:length(data.Xdrc)-1);
    plot(wt,newDrcX,'k-',wt,newDrcY,'r-');
    if i == 1, ylabel('degree'); else; set(gca,'yticklabel',[]); end
    if i == 4, legend({'\theta_i^H','\theta_i^W'},'Location','SouthEast'); end
    xlim([pk-30/8 pk+30/8]);
    ylim([-80 360]);
    grid on;
end
savePics('figurespaper/dataslice.pdf',saveOn,'beamer');


%% Correlation
figure(13); clf;
iV = 0:30;
corrxy = arrayfun(@(x)(corr(data.lap.X(1:end-x),data.lap.Y(1+x:end))),iV);
corrx = arrayfun(@(x)(corr(data.lap.X(1:end-x),data.lap.X(1+x:end))),iV);
corry = arrayfun(@(x)(corr(data.lap.Y(1:end-x),data.lap.Y(1+x:end))),iV);

t = tiledlayout(1,2,'TileSpacing','none','Padding','compact');

nexttile;
plot(iV*3,corrx,'k-',iV*3,corry,'b-',iV*3,corrxy,'r-');
xlabel('lag (hours)'); ylabel('Correlation');
% title('');
legend({'$\mathrm{Corr}\left(H^L_{S,i},H^L_{s,i+\mathrm{lag}}\right)$',...
    '$\mathrm{Corr}\left(W^L_{S,i},W^L_{s,i+\mathrm{lag}}\right)$',...
    '$\mathrm{Corr}\left(H^L_{S,i},W^L_{s,i+\mathrm{lag}}\right)$'},'interpreter','latex');
grid on;
xlim([0 iV(end)*3+1.5]);
ylim([-0.2 1]);

iV = 0:5;
% corrxy = arrayfun(@(x)(corr(mod(diff(newDrcX(1:end-x))+180,360)-180,mod(diff(newDrcY(1+x:end))+180,360)-180)),iV);
% corrx = arrayfun(@(x)(corr(mod(diff(newDrcX(1:end-x))+180,360)-180,mod(diff(newDrcX(1+x:end))+180,360)-180)),iV);
% corry = arrayfun(@(x)(corr(mod(diff(newDrcY(1:end-x))+180,360)-180,mod(diff(newDrcY(1+x:end))+180,360)-180)),iV);
corrxy = arrayfun(@(x)(corr(mod(diff(data.Xdrc(1:end-x))+180,360)-180,mod(diff(data.Ydrc(1+x:end))+180,360)-180)),iV);
corrx = arrayfun(@(x)(corr(mod(diff(data.Xdrc(1:end-x))+180,360)-180,mod(diff(data.Xdrc(1+x:end))+180,360)-180)),iV);
corry = arrayfun(@(x)(corr(mod(diff(data.Ydrc(1:end-x))+180,360)-180,mod(diff(data.Ydrc(1+x:end))+180,360)-180)),iV);

gammat = mod(data.Ydrc - data.Xdrc + 180,360) - 180;
corrgamma = arrayfun(@(x)(corr(mod(diff(gammat(1:end-x))+180,360)-180,mod(diff(gammat(1+x:end))+180,360)-180)),iV);
% plot();%,'Color',[0.6,0.6,0.6])

nexttile;
plot(iV*3,corrx,'k-',iV*3,corry,'b-',iV*3,corrxy,'r-',iV*3,corrgamma,'g-');
xlabel('lag (hours)')
xlim([0 iV(end)*3+0.5]);
ylim([-0.2 1]);
set(gca,'yticklabel',[]);
% title('');
grid on;
legend({'$\mathrm{Corr}\left(\Delta\theta^H_{S,i},\Delta\theta^H_{s,i+\mathrm{lag}}\right)$',...
    '$\mathrm{Corr}\left(\Delta\theta^W_{S,i},\Delta\theta^W_{s,i+\mathrm{lag}}\right)$',...
    '$\mathrm{Corr}\left(\Delta\theta^H_{S,i},\Delta\theta^W_{s,i+\mathrm{lag}}\right)$',...
    '$\mathrm{Corr}\left(\gamma_i,\gamma_{i+\mathrm{lag}}\right)$'},'interpreter','latex');

savePics('figurespaper/correlation.pdf',saveOn,'beamer',10,10*1/3);

%% Trajectory plots
meths = {'EVAR','MMEM','VAR','HM'};
for imeth = 1 %1:4
    meth = meths{imeth};
    switch meth
        case 'HM'
            orders = 1;
        otherwise
            orders = 1:6;
    end
    % meth = 'EVAR'; order = 4;
    for iorder = 4 %1:length(orders)
        AH = zeros(0,70);
        AW = zeros(0,70);
        order = orders(iorder);
        stormpeak = 12;
        % compareX = true;
        compareX = true;
        plotAll = true;
        methcompare = 'data';
        cnt = 0;
        for iStorm = 1:opts.nStormsSim
            switch meth
                case 'data'
                    if iStorm > length(data.storms)
                        continue
                    end
                    currStormX = data.X(data.storms{iStorm}.index);
                    currStormY = data.Y(data.storms{iStorm}.index);
                case 'HM'
                    currStormX = stormsOrg.(meth){1}{iStorm}.originalMargins(:,1);
                    currStormY = stormsOrg.(meth){1}{iStorm}.originalMargins(:,2);
                otherwise
                    currStormX = stormsOrg.(meth){order}{iStorm}.originalMargins(:,1);
                    currStormY = stormsOrg.(meth){order}{iStorm}.originalMargins(:,2);
            end
            [M,I] = max(currStormX);
            L = length(currStormX);
            if abs(M - stormpeak) > 0.5
                continue
            else
                if compareX
                    myVar = currStormX;
                end
                myVarH = currStormX;
                myVarW = currStormY;
                
                cnt = cnt + 1;
                %         plot(-I:-I+L-1,currStormX,'r-','LineWidth',1);
                %         plot(-I:-I+L-1,currStormY,'k-','LineWidth',1);
                AH(end+1,35-I:34-I+L) = myVarH;
                AH(end,1:34-I) = -inf;
                AH(end,35-I+L:end) = -inf;
                
                AW(end+1,35-I:34-I+L) = myVarW;
                AW(end,1:34-I) = -inf;
                AW(end,35-I+L:end) = -inf;
            end
        end
        
        AHc = zeros(0,70);
        AWc = zeros(0,70);
        cnt = 0;
        for iStorm = 1:opts.nStormsSim
            switch methcompare
                case 'data'
                    if iStorm > length(data.storms)
                        continue
                    end
                    currStormX = data.X(data.storms{iStorm}.index);
                    currStormY = data.Y(data.storms{iStorm}.index);
                otherwise
                    currStormX = stormsOrg.(methcompare){order}{iStorm}.originalMargins(:,1);
                    currStormY = stormsOrg.(methcompare){order}{iStorm}.originalMargins(:,2);
            end
            [M,I] = max(currStormX);
            L = length(currStormX);
            
            if abs(M - stormpeak) > 0.5
                continue
            else
                if compareX
                    myVar = currStormX;
                end
                myVarH = currStormX;
                myVarW = currStormY;
                
                cnt = cnt + 1;
                %         plot(-I:-I+L-1,currStormX,'r-','LineWidth',1);
                %         plot(-I:-I+L-1,currStormY,'k-','LineWidth',1);
                AHc(end+1,35-I:34-I+L) = myVarH;
                AHc(end,1:35-I) = -inf;
                AHc(end,35-I+L:end) = -inf;
                AWc(end+1,35-I:34-I+L) = myVarW;
                AWc(end,1:35-I) = -inf;
                AWc(end,35-I+L:end) = -inf;
            end
        end
        
        col = {'r','m'};
        
        % Plot model trajectories HS
        % [~,I] = min(mean(isinf(AH),2));
        lowerQuantile = 0.1; upperQuantile = 0.9;
        
        
        xl = (-33:36)*3;
        
        %
        figure(33); clf;
        t = tiledlayout(3,3,'TileSpacing','none','Padding','compact');
        % t = tiledlayout(6,6,'TileSpacing','none','Padding','compact');
        
        nexttile(1);
        % nexttile([2 2]);
        plot(xl,AH','k-');
        hold on;
        % plot(xl,AH(I,:),'r-','LineWidth',2); hold on;
        ylabel('H_S');
        grid on;
        title(sprintf('%s(%d)',meth,order));
        axis tight;
        xlim([-50 70]);
        xx = xlim;
        % set(gca,'xticklabel',{[]});
        ylim([4 stormpeak+1]);
        set(gca,'xticklabel',{[]});
        % xlim([-45 60]);
        
        
        nexttile(2);
        % nexttile([2 2]);
        plot(xl,AHc','r-'); hold on;
        % plot(xl,AHc(I,:),'r-','LineWidth',2); hold on;
        grid on;
        title('Data');
        xlim([xx(1) xx(2)]);
        ylim([4 stormpeak+1]);
        set(gca,'yticklabel',{[]});
        set(gca,'xticklabel',{[]});
        
        
        nexttile(3);
        % nexttile([2 2]);
        plot(xl,quantile(AH,0.5),'k-','LineWidth',1);
        hold on;
        plot(xl,quantile(AH,[lowerQuantile;upperQuantile]),'k--','LineWidth',1);
        plot(xl,quantile(AHc,0.5),'r-','LineWidth',1);
        plot(xl,quantile(AHc,[lowerQuantile;upperQuantile]),'r--','LineWidth',1);
        grid on
        % ylim([0 30]);
        % axis tight;
        xlim([-32 45]);
        ylim([4 stormpeak+1]);
        title('Summary');
        % xlim([xx(1) xx(2)]);
        set(gca,'yticklabel',{[]});
        set(gca,'xticklabel',{[]});
        
        % Plot model trajectories Ws
        nexttile(4);
        % nexttile([2 2]);
        plot(xl,AW','k-'); hold on;
        % plot(xl,AW(I,:),'r-','LineWidth',2); hold on;
        ylabel('W_s');
        grid on;
        xlim([xx(1) xx(2)]);
        ylim([0 30]);
        
        % Plot data trajectories HS
        % [~,I] = min(mean(isinf(AHc),2));
        
        
        % Plot data trajectories Ws
        nexttile(5);
        % nexttile([2 2]);
        plot(xl,AWc','r-'); hold on;
        % plot(xl,AWc(I,:),'r-','LineWidth',2); hold on;
        grid on;
        xlim([xx(1) xx(2)]);
        ylim([0 30]);
        set(gca,'yticklabel',{[]});
        
        % Plot summaries of HS trajectories
        lowerQuantile = 0.1; upperQuantile = 0.9;
        
        
        
        % Plot summaries of Ws trajectories
        nexttile(6);
        % nexttile([2 2]);
        plot(xl,quantile(AW,0.5),'k-','LineWidth',1);
        hold on;
        plot(xl,quantile(AW,[lowerQuantile;upperQuantile]),'k--','LineWidth',1);
        plot(xl,quantile(AWc,0.5),'r-','LineWidth',1);
        plot(xl,quantile(AWc,[lowerQuantile;upperQuantile]),'r--','LineWidth',1);
        % xlim([xx(1) xx(2)]);
        grid on
        % title(meth)
        % ylim([0 1.2]);
        % axis tight;
        ylim([0 30]);
        xlim([-32 45]);
        % title('Summary');
        % xlim([xx(1) xx(2)]);
        set(gca,'yticklabel',{[]});
        % title(sprintf('Compare %s (black) vs %s (blue)',meth,methcompare));
        
        %
        % figure(34); clf;
        nexttile(7,[1 3]);
        % nexttile(26,[2 4]);
        plot(xl,1-mean(isinf(AH)),'k-','LineWidth',1);
        hold on;
        plot(xl,1-mean(isinf(AHc)),'r-','LineWidth',1)
        grid on
        xlim([xx(1) xx(2)]);
        % title('Existence');
        ylabel('Probability');
        xlabel('Hours (relative to the peak)');
        
        savePics(sprintf('figurespaper/trajectories%s%d.pdf',meth,order),saveOn,'beamer');
    end
end

%% Transform HM trajectories onto standard amargins
for iStorm = 1:opts.nStormsSim
    currStorm = stormsOrg.HM{1}{iStorm};
    [currStormX, currStormY] = MarginsHM(currStorm,[MM_HS,MM_WS]);
    stormsOrg.HM{1}{iStorm}.standardMargins = [currStormX, currStormY];
end


%% Chi(u) estimates / plot 3

figure(37); clf;
lags = [1,4];
nlags = length(lags);

t = tiledlayout(nlags,3,'TileSpacing','none','Padding','compact');

% lags = [14,8];
ul = 3:0.5:6;

colors = {'k','r','b','g'};
linestyle = {'-','--','-.',':','-'};

for ilag = 1:length(lags)
    lag = lags(ilag);
    meths = {'data','EVAR','MMEM','HM'};%,'VAR','HM'};
    for imeth = 1:4
        meth = meths{imeth};
        switch meth
            case {'HM','data'}
                orders = 1;
            otherwise
                %                 orders = 1:6;
                orders = [1,4];
        end
        % meth = 'EVAR'; order = 4;
        for iorder = 1:min(4,length(orders))
            
            order = orders(iorder);
            cnt = 0;
            chiHW = zeros(2,length(ul));
            chiW = zeros(2,length(ul));
            chiH = zeros(2,length(ul));
            
            for iu = 1:length(ul)
                u = ul(iu);
                for iStorm = 1:opts.nStormsSim
                    switch meth
                        case 'data'
                            if iStorm > length(data.storms)
                                continue
                            end
                            currStormX = data.lap.X(data.storms{iStorm}.index);
                            currStormY = data.lap.Y(data.storms{iStorm}.index);
                        case 'HM'
                            currStormX = stormsOrg.(meth){order}{iStorm}.standardMargins(:,1);
                            currStormY = stormsOrg.(meth){order}{iStorm}.standardMargins(:,2);
                            %                             [currStormX, currStormY] = MarginsHM(currStorm,[MM_HS,MM_WS]);
                        otherwise
                            currStormX = stormsOrg.(meth){order}{iStorm}.standardMargins(:,1);
                            currStormY = stormsOrg.(meth){order}{iStorm}.standardMargins(:,2);
                    end
                    
                    for i = 1:length(currStormX)
                        if i + lag > length(currStormX)
                            flagX = false;
                            flagY = false;
                        else
                            if currStormX(i+lag) > u
                                flagX = true;
                            else
                                flagX = false;
                            end
                            if currStormY(i+lag) > u
                                flagY = true;
                            else
                                flagY = false;
                            end
                        end
                        
                        chiH(1,iu) = chiH(1,iu) + (currStormX(i) > u & flagX);
                        chiH(2,iu) = chiH(2,iu) + (currStormX(i) > u & ~flagX);
                        
                        chiW(1,iu) = chiW(1,iu) + (currStormY(i) > u & flagY);
                        chiW(2,iu) = chiW(2,iu) + (currStormY(i) > u & ~flagY);
                        
                        chiHW(1,iu) = chiHW(1,iu) + (currStormX(i) > u & flagY);
                        chiHW(2,iu) = chiHW(2,iu) + (currStormX(i) > u & ~flagY);
                    end
                end
            end
            
            if ~strcmp(meth,'data')
                nexttile((ilag-1)*3+1);
                %                 nexttile;
                title(sprintf('$\\chi_H(u,%d)$',lag), 'interpreter', 'latex');
                plot(ul,chiH(1,:)./sum(chiH),'Color',colors{imeth},'LineStyle',linestyle{iorder},'LineWidth',1.5);
                hold on;
                grid on;
                %             ylim([0 .4]);
                
                nexttile((ilag-1)*3+2);
                %                 nexttile;
                title(sprintf('$\\chi_{HW}(u,%d)$',lag), 'interpreter', 'latex');
                
                plot(ul,chiHW(1,:)./sum(chiHW),'Color',colors{imeth},'LineStyle',linestyle{iorder},'LineWidth',1.5);
                hold on;
                grid on;
                %             ylim([0 .3]);
                
                nexttile((ilag-1)*3+3);
                %                 nexttile;
                title(sprintf('$\\chi_W(u,%d)$',lag), 'interpreter', 'latex');
                
                plot(ul,chiW(1,:)./sum(chiW),'Color',colors{imeth},'LineStyle',linestyle{iorder},'LineWidth',1.5);
                hold on;
                grid on;
                %             ylim([0 .3]);
            end
            if strcmp(meth,'data') % then add confidence bounds
                n = chiH(1,:) + chiH(2,:);
                p = chiH(1,:) ./ sum(chiH);
                std = sqrt(p.*(1-p)./n);
                
                nexttile((ilag-1)*3+1);
                %                 nexttile;
                patch([ul,ul(end:-1:1)],[p + 2*std,max(0,p(end:-1:1)-2*std(end:-1:1))],colors{imeth},'FaceAlpha',0.15,'EdgeAlpha',0.2);
                hold on;
                n = chiHW(1,:) + chiHW(2,:);
                p = chiHW(1,:) ./ sum(chiHW);
                std = sqrt(p.*(1-p)./n);
                
                nexttile((ilag-1)*3+2);
                %                 nexttile;
                p1 = patch([ul,ul(end:-1:1)],[p + 2*std,max(0,p(end:-1:1)-2*std(end:-1:1))],colors{imeth},'FaceAlpha',0.15,'EdgeAlpha',0.2);
                hold on;
                n = chiW(1,:) + chiW(2,:);
                p = chiW(1,:) ./ sum(chiW);
                std = sqrt(p.*(1-p)./n);
                
                nexttile((ilag-1)*3+3);
                %                 nexttile;
                patch([ul,ul(end:-1:1)],[p + 2*std,max(0,p(end:-1:1)-2*std(end:-1:1))],colors{imeth},'FaceAlpha',0.15,'EdgeAlpha',0.2);
                hold on;
            end
            
            
        end
    end
    %     savePics(sprintf('figurespaper/trajectories%s%d.pdf',meth,order),saveOn,'beamer');
end
savePics('modelfit3_1.pdf',saveOn,'paper',10,10*2/3)


