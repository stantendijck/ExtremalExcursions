%% Plots for simulation study
d = @(x,y)(abs(x-y)./max(x,y));
average = @(x)(mean(x));

%% Response picture

variable = impact;
p2 = linspace(0.97,0.999,20);
q = cell(length(p2),1);
qHM = cell(length(p2),1);
qEmp = cell(length(p2),1);
qEmpTrain = cell(length(p2),1);

funs = {'max','sum'};
% HsThr = [5,7];
% c = [0.4, 0.2];
MdlPrm = {[0.2, 0.3], [6, 7]};
HsThr = MdlPrm{2};
c = MdlPrm{1};
% HsThr = [6,8];
% c = [0.3,0.15];

figure(30); clf;
t = tiledlayout(2,2,'TileSpacing','compact','Padding','none');
for iPrm = 1:2
    for iFun = 1:2
        nexttile
        
        clear v;
        nOrd = opts.maxLag-1;
        for ip = 1:length(p2)
            nMeth = 3;
            q{ip} = nan(nBootstraps,nMeth,nOrd);
            qEmp{ip} = nan(nBootstraps,1);
            qEmpTrain{ip} = nan(nBootstraps,1);
            qHM{ip} = nan(nBootstraps,1);
            for iBt = 1:nBootstraps
                for iMeth = 1:nMeth
                    for iOrd = 1:nOrd
                        A = variable{iBt}{iPrm,iFun}.models{iMeth,iOrd};
                        q{ip}(iBt,iMeth,iOrd) = quantile(A,p2(ip));
                        
                    end
                end
                A = variable{iBt}{iPrm,iFun}.test{1};
                qEmp{ip}(iBt) = quantile(A,p2(ip));
                
                A = variable{iBt}{iPrm,iFun}.train{1};
                qEmpTrain{ip}(iBt) = quantile(A,p2(ip));
                
                A = variable{iBt}{iPrm,iFun}.HM{1};
                qHM{ip}(iBt) = quantile(A, p2(ip));
                
            end
            
        end
        
        %w = cell2mat(arrayfun(@(y)(mean(cell2mat(arrayfun(@(x)(d(q{x}(:,ceil(y/6),y-(ceil(y/6)-1)*6),qEmp{x})),1:length(p2),'UniformOutput',false)),2)),1:18,'UniformOutput',false));
        %w = mean(cell2mat(arrayfun(@(x)(d(q{x}(:,ceil(y/nOrd),y-(ceil(y/nOrd)-1)*nOrd),qEmp{x})),1:length(p2),'UniformOutput',false)),2);
        w = cell2mat(arrayfun(@(y)(mean(cell2mat(arrayfun(@(x)(d(q{x}(:,ceil(y/nOrd),y-(ceil(y/nOrd)-1)*nOrd),qEmp{x})),1:length(p2),'UniformOutput',false)),2)),1:(3*nOrd),'UniformOutput',false));
        wHM = mean(cell2mat(arrayfun(@(x)(d(qHM{x},qEmp{x})),1:length(p2),'UniformOutput',false)),2);
        wTrain = mean(cell2mat(arrayfun(@(x)(d(qEmpTrain{x},qEmp{x})),1:length(p2),'UniformOutput',false)),2);
        wTot = [wHM,w,wTrain];
        mv = mean(wTot,1);
        
        B = 1000;
        mvB = nan(B,nOrd*3+2);
        mvB(1,:) = mv;
        for iB = 2:B
            mvB(iB,:) = mean(wTot(randsample(nBootstraps,nBootstraps,true),:),1);
        end
        qv = quantile(mvB,[0.1 0.9]);
        
        plot(1:nOrd,mv(1)*ones(1,nOrd),'k--','LineWidth',2);
        hold on;
        plot(1:nOrd,mv(2:(1+nOrd)),'r--','LineWidth',2);
        plot(1:nOrd,mv((2+nOrd):(1+2*nOrd)),'m:','LineWidth',2);
        plot(1:nOrd,mv((2+2*nOrd):(1+3*nOrd)),'b--','LineWidth',2);
        
        xconf = [1:nOrd,nOrd:-1:1] ;
        yconf = [qv(1,1)*ones(1,nOrd) qv(2,1)*ones(1,nOrd)];
        p = patch(xconf,yconf,'k','FaceAlpha',0.2,'EdgeAlpha',0);
        
        yconf = [qv(1,2:(1+nOrd)) qv(2,(1+nOrd):-1:2)];
        p = patch(xconf,yconf,'r','FaceAlpha',0.2,'EdgeAlpha',0);
        
        yconf = [qv(1,(2+2*nOrd):(1+3*nOrd)) qv(2,(1+3*nOrd):-1:(2+2*nOrd))];
        p = patch(xconf,yconf,'b','FaceAlpha',0.2,'EdgeAlpha',0);
        
        title(sprintf('R^{%s}(%d, %.1f)',funs{iFun},HsThr(iPrm),c(iPrm)));
        grid on
        
        if iPrm == 2 && iFun == 1
            xlabel('Model order');
        end
        set(gca,'xtick',1:6);
        
    end
end
legend({'HM','EVAR','EVAR_0','MMEM'});
savePics('figurespaper/response.pdf',saveOn,'beamer');


