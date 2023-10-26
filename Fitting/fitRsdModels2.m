function mdl = fitRsdModels2(mdl, printOn)
% Fit residual models for input mdl - struct

if nargin < 2
    printOn = true;
end

if printOn
    fprintf('Fitting marginal and copula models for BS sample\n');
end

for i = 1:length(mdl)
    expl = mdl{i}.data.explanatory;
    mdl{i}.rsd.data = nan(length(expl), length(mdl{i}.phat));
    mdl{i}.rsd.marghat = cell(length(mdl{i}.phat), 1);
    dataUnif = nan(length(expl), length(mdl{i}.phat));
    for j = 1:length(mdl{i}.phat)
        %% Model parameters
        p = mdl{i}.phat{j};
        resp = mdl{i}.data.response(:,j);
        
        %% Calculate residuals
        mdl{i}.rsd.data(:,j) = (resp - p(1) * expl) ./ (expl .^ p(2));
            
        %% Marginal model
        XI = linspace(min(mdl{i}.rsd.data(:,j))-1.5,max(mdl{i}.rsd.data(:,j))+1.5,1000)';
        f = ksdensity(mdl{i}.rsd.data(:,j), XI, 'Function', 'cdf');
        
        begin = find(f==0);
        if any(begin); begin = begin(end); else; begin = 1; end
        
        einde = find(f==1);
        if any(einde); einde = einde(1);else; einde = length(f);end
        
        mdl{i}.rsd.marghat{j} = [f(begin:einde), XI(begin:einde)];
        
        %% Transform to uniform margins
        dataUnif(:,j) = interp1(XI(begin:einde), f(begin:einde), mdl{i}.rsd.data(:,j));
    end
    
    if size(dataUnif,2) > 1
%         [rhoHAT, nuHAT] = copulafit('t', dataUnif, 'Method', 'ApproximateML');

        %% Save the parameters
%         mdl{i}.rsd.rhoHAT = rhoHAT;
%         mdl{i}.rsd.nuHAT = nuHAT;
    end
end

end






