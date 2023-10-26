function SCV = output_diagnostics(data, methods, stormsOrg, indMethod, indOrder)
%%
sim = stormsOrg.(methods{indMethod}){indOrder};
maxModelOrder = length(stormsOrg.EVAR);
nStorms = length(sim);
% indMethod = 1; indOrder = 3;

%% Storm trajectory

% xL = linspace(0,360,1000);
% yL = arrayfun(@(x)(transform2Org(Laplace_iCDF(0.95), x, MM_HS)),xL);

% Let's say storm has impact if Hs > 10 and additional impact if WS > 20
% impact: 
% if Hs < 8:                Impact = 0
% if Hs > 8 and WS < 20:    Impact = A + B*(Hs - 10)^2
% if Hs > 8 and WS > 20:    Impact = A + B*(Hs - 10)^2 + C*(Ws - 20)^2

A = 20;
B = 1;
C = 1/5;
HsThr = 10;
WsThr = 20;

figure(13); clf;

impact = zeros(nStorms,length(methods),maxModelOrder);
col = jet(length(methods)*maxModelOrder);
leg = cell(0,0);
for iMethod = 1:length(methods)
    for iOrder = 1:maxModelOrder
        nStorms = length(stormsOrg.(methods{iMethod}){iOrder});
        for iStorm = 1:nStorms
            currStorm = stormsOrg.(methods{iMethod}){iOrder}{iStorm}.originalMargins;
        %     currStorm = [data.X(data.storms{iStorm}.index),data.Y(data.storms{iStorm}.index)];
            I = currStorm(:,1) > 10;
            if any(I)
                impact(iStorm,iMethod,iOrder) = sum(A + B*(currStorm(I,1) - HsThr).^2 + C*max(0,currStorm(I,2) - WsThr).^2);
            end
        end
        [F,X] = ecdf(impact(:,iMethod,iOrder),'Function','survivor');
        plot(X,F,'Color',col((iMethod-1)*maxModelOrder + iOrder,:)); hold on;
        leg{end+1} = sprintf('%s(%d)',methods{iMethod},iOrder);
    end
end
nStorms = length(data.storms);
impactEmp = zeros(nStorms,1);
for iStorm = 1:nStorms
    currStorm = [data.X(data.storms{iStorm}.index),data.Y(data.storms{iStorm}.index)];
    I = currStorm(:,1) > 10;
    if any(I)
        impactEmp(iStorm) = sum(A + B*(currStorm(I,1) - HsThr).^2 + C*max(0,currStorm(I,2) - WsThr).^2);
    end
end
[F,X] = ecdf(impactEmp,'Function','survivor');
plot(X,F,'k-','LineWidth',2);
leg{end+1} = sprintf('data');
legend(leg);
ylim([0 0.1]);
set(gca,'YScale','log')

%
muGPD = A; % needs to be at least A
like = @(data,p)(GPD_like(data,muGPD,p(1),p(2)));
p0 =  [1;0.1];
mydata = impactEmp(impactEmp > muGPD);
mcmc1 = adptMCMC(mydata,like,p0,1e4);
pltMCMC(mcmc1,14,{'\sigma','\xi'})
% phat = fminsearch(@(p)(like(mydata,p)),mcmc1.MLE);



like = @(data,p)(GPD_like(data,muGPD,p(1),p(2)));
p0 =  [1;0.1];
mydata = impact(impact(:,indMethod,indOrder)>muGPD,indMethod,indOrder);
mcmc2 = adptMCMC(mydata,like,p0,1e4);
pltMCMC(mcmc2,17,{'\sigma','\xi'})
% fminsearch(@(p)(like(mydata,p)),p0)
fprintf('MLE from data: xi = %.2f\n',mcmc1.MLE(2));
fprintf('MLE from model: xi = %.2f\n',mcmc2.MLE(2));

%%
% rho = @(p,a,b)((a-b).*(p-(a<b)));
% 
% estimates = reshape(quantile(impact,pc),length(methods),maxModelOrder)';
% SCV = nan(maxModelOrder,length(methods));
% for iMethod = 1:length(methods)
%     for iOrder = 1:maxModelOrder
% %         estimates(iMethod,iOrder)
%         SCV(iOrder,iMethod) = mean(rho(pc,estimates(iOrder,iMethod),impactEmp));
%     end
% end


end


