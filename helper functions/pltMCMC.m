function pltMCMC(mcmc,figCntr,varNms)

if nargin < 2
    figCntr = 10;
end

if nargin < 3
    varNms = cell(size(mcmc.p,1),1);
    for i = 1:size(mcmc.p,1)
        varNms{i} = sprintf('P_%d',i);
    end
end

[~,I] = min(mcmc.like);

figure(figCntr);clf;
% plot(mcmc.p(1,mcmc.burnIn+mcmc.adpItr:end)');
plot(mcmc.p(:,mcmc.burnIn+mcmc.adpItr:end)');
hold on;
% plot(mcmc.alpha2)
% plot(mcmc.p(3,mcmc.burnIn+mcmc.adpItr:end)');
plot(I-mcmc.burnIn-mcmc.adpItr+1,mcmc.p(:,I)','r*');
legend(varNms);

make_it_tight = length(varNms) == 5;
subplot2 = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.18 0.1], [0.05 0.03]);




figure(figCntr+2);clf;
subplot(1,2,1);
plot(mcmc.like - min(mcmc.like) + 1); set(gca,'YScale','log');
hold on;
plot((mcmc.burnIn+mcmc.adpItr)*ones(1,2),ylim,'k--');
plot(I,mcmc.like(I) - min(mcmc.like) + 1,'r*');
ylabel('(scaled) like');
subplot(1,2,2);
plot(mcmc.accept_rate);
hold on;
plot((mcmc.burnIn+mcmc.adpItr)*ones(1,2),ylim,'k--');
ylabel('accept rate');

figure(figCntr+1);clf;
cnt = 1;
for i = 1:size(mcmc.p,1)
    for j = 1:size(mcmc.p,1)
        if make_it_tight
            subplot2(size(mcmc.p,1),size(mcmc.p,1),cnt);
        else
            subplot(size(mcmc.p,1),size(mcmc.p,1),cnt);
        end
        if i ~= j
            plot(mcmc.p(j,mcmc.burnIn+mcmc.adpItr:end),mcmc.p(i,mcmc.burnIn+mcmc.adpItr:end),'k.');
            hold on;
            plot(mcmc.MLE(j),mcmc.MLE(i),'r*');
            xlabel(varNms{j});
            ylabel(varNms{i});
        else
            histogram(mcmc.p(i,mcmc.burnIn+mcmc.adpItr:end),'FaceColor','black');
            xlabel(varNms{i});
        end
        cnt = cnt + 1;
    end
end

end