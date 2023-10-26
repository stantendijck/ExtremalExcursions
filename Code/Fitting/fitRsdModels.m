function rsd = fitRsdModels(mdl, printOn, dens)
% Fit residual models for input mdl - struct

if nargin < 2
    printOn = true;
end
if nargin < 3
    dens = 'TGN';
end

rsd = struct();
if printOn
    fprintf('Fitting marginal and copula models for BS sample\n');
end

%% Collect type 1 residuals
type1Xrsd = cell2mat(cellfun(@(x)(x.MLERsd),mdl.type1.X,'UniformOutput',false)');
type1Yrsd = cell2mat(cellfun(@(x)(x.MLERsd),mdl.type1.Y,'UniformOutput',false)');
nX = size(type1Xrsd,2);
nY = size(type1Yrsd,2);

%% Fit marginal models and transform residuals to uniform margins
type1XrsdU = nan(size(type1Xrsd));
type1YrsdU = nan(size(type1Yrsd));
phat = struct();

switch dens
    case 'TGN'
        phat.type1 = zeros(nX + nY, 6);
    otherwise
        phat.type1 = cell(nX + nY, 1);
end

for jX = 1:nX
    switch dens
        case 'TGN'
            p0 = fminsearch(@(p)(GenNormal_like(type1Xrsd(:,jX),p(1),p(2),p(3))),[0,1,1]);
            phat.type1(jX,:) = fminsearch(@(p)(TwoSidedGenNormal_like(type1Xrsd(:,jX),p(1),p(2:3),p(4:5),p(6))),[p0([1,2,2,3,3]),0.5]);
            type1XrsdU(:,jX) = TwoSidedGenNormal_CDF(type1Xrsd(:,jX),phat.type1(jX,1),phat.type1(jX,2:3),phat.type1(jX,4:5),phat.type1(jX,6));
        otherwise
            XI = linspace(min(type1Xrsd(:,jX))-0.5,max(type1Xrsd(:,jX))+0.5,1000)';
            f = ksdensity(type1Xrsd(:,jX), XI, 'Function', 'cdf');
            begin = find(f==0);
            if any(begin)
                begin = begin(end);
            else
                begin = 1;
            end
            einde = find(f==1);
            if any(einde)
                einde = einde(1);
            else
                einde = length(f);
            end
            phat.type1{jX} = [f(begin:einde), XI(begin:einde)];
            type1XrsdU(:,jX) = interp1(XI(begin:einde), f(begin:einde), type1Xrsd(:,jX));
%             type1XrsdU(:,jX) = ksdensity(type1Xrsd(:,jX), type1Xrsd(:,jX), 'Function', 'cdf');
    end
end
for jY = 1:nY
    switch dens
        case 'TGN'
            p0 = fminsearch(@(p)(GenNormal_like(type1Yrsd(:,jY),p(1),p(2),p(3))),[0,1,1]);
            phat.type1(nX+jY,:) = fminsearch(@(p)(TwoSidedGenNormal_like(type1Yrsd(:,jY),p(1),p(2:3),p(4:5),p(6))),[p0([1,2,2,3,3]),0.5]);
            type1YrsdU(:,jY) = TwoSidedGenNormal_CDF(type1Yrsd(:,jY),phat.type1(nX+jY,1),phat.type1(nX+jY,2:3),phat.type1(nX+jY,4:5),phat.type1(nX+jY,6));
        otherwise
            XI = linspace(min(type1Yrsd(:,jY))-0.5,max(type1Yrsd(:,jY))+0.5,1000)';
            f = ksdensity(type1Yrsd(:,jY), XI, 'Function', 'cdf');
            begin = find(f==0);
            if any(begin)
                begin = begin(end);
            else
                begin = 1;
            end
            einde = find(f==1);
            if any(einde)
                einde = einde(1);
            else
                einde = length(f);
            end
            phat.type1{nX+jY} = [f(begin:einde), XI(begin:einde)];
            type1YrsdU(:,jY) = interp1(XI(begin:einde), f(begin:einde), type1Yrsd(:,jY));
%             type1YrsdU(:,jY) = ksdensity(type1Yrsd(:,jY), type1Yrsd(:,jY), 'Function', 'cdf');
    end
end

%% Fit Copula model
% [rhoHAT,nuHAT] = copulafit('t', [type1XrsdU, type1YrsdU], 'Method', 'ApproximateML');

%% Save the parameters
% rsd.type1.rhoHAT = rhoHAT;
% rsd.type1.nuHAT = nuHAT;
rsd.type1.margHAT = phat.type1;

%% Type n
nTypeN = length(mdl.typen);
rsd.typen = cell(nTypeN, 1);
phat.typen = cell(nTypeN, 1);
for iN = 1:nTypeN
    %% Collect type n residuals
    typeNXrsd = cell2mat(cellfun(@(x)(x.MLERsd),mdl.typen{iN}.X,'UniformOutput',false)');
    typeNYrsd = cell2mat(cellfun(@(x)(x.MLERsd),mdl.typen{iN}.Y,'UniformOutput',false)');
    nX = size(typeNXrsd,2);
    nY = size(typeNYrsd,2);
    
    %% Fit marginal models and transform residuals to uniform margins
    typeNXrsdU = nan(size(typeNXrsd));
    typeNYrsdU = nan(size(typeNYrsd));
    
    switch dens
        case 'TGN'
            phat.typen{iN} = zeros(nX + nY, 6);
        otherwise
            phat.typen{iN} = cell(nX + nY, 1);
    end
    
    for jX = 1:nX
        switch dens
            case 'TGN'
                p0 = fminsearch(@(p)(GenNormal_like(typeNXrsd(:,jX),p(1),p(2),p(3))),[0,1,1]);
                phat.typen{iN}(jX,:) = fminsearch(@(p)(TwoSidedGenNormal_like(typeNXrsd(:,jX),p(1),p(2:3),p(4:5),p(6))),[p0([1,2,2,3,3]),0.5]);
                typeNXrsdU(:,jX) = TwoSidedGenNormal_CDF(typeNXrsd(:,jX),phat.typen{iN}(jX,1),phat.typen{iN}(jX,2:3),phat.typen{iN}(jX,4:5),phat.typen{iN}(jX,6));
            otherwise
                XI = linspace(min(typeNXrsd(:,jX))-0.5,max(typeNXrsd(:,jX))+0.5,1000)';
                f = ksdensity(typeNXrsd(:,jX), XI, 'Function', 'cdf');
                begin = find(f==0);
                if any(begin)
                    begin = begin(end);
                else
                    begin = 1;
                end
                einde = find(f==1);
                if any(einde)
                    einde = einde(1);
                else
                    einde = length(f);
                end
                phat.typen{iN}{jX} = [f(begin:einde), XI(begin:einde)];
                typeNXrsdU(:,jX) = interp1(XI(begin:einde), f(begin:einde), typeNXrsd(:,jX));
%                 typeNXrsdU(:,jX) = ksdensity(typeNXrsd(:,jX),typeNXrsd(:,jX),'Function','cdf');
        end
    end
    
    for jY = 1:nY
        switch dens
            case 'TGN'
                p0 = fminsearch(@(p)(GenNormal_like(typeNYrsd(:,jY),p(1),p(2),p(3))),[0,1,1]);
                phat.typen{iN}(nX+jY,:) = fminsearch(@(p)(TwoSidedGenNormal_like(typeNYrsd(:,jY),p(1),p(2:3),p(4:5),p(6))),[p0([1,2,2,3,3]),0.5]);
                typeNYrsdU(:,jY) = TwoSidedGenNormal_CDF(typeNYrsd(:,jY),phat.typen{iN}(nX+jY,1),phat.typen{iN}(nX+jY,2:3),phat.typen{iN}(nX+jY,4:5),phat.typen{iN}(nX+jY,6));
            otherwise
                XI = linspace(min(typeNYrsd(:,jY))-0.5,max(typeNYrsd(:,jY))+0.5,1000)';
                f = ksdensity(typeNYrsd(:,jY), XI, 'Function', 'cdf');
                begin = find(f==0);
                if any(begin)
                    begin = begin(end);
                else
                    begin = 1;
                end
                einde = find(f==1);
                if any(einde)
                    einde = einde(1);
                else
                    einde = length(f);
                end
                
                phat.typen{iN}{nX+jY} = [f(begin:einde), XI(begin:einde)];
                typeNYrsdU(:,jY) = interp1(XI(begin:einde), f(begin:einde), typeNYrsd(:,jY));
%                 typeNYrsdU(:,jY) = ksdensity(typeNYrsd(:,jY),typeNYrsd(:,jY),'Function','cdf');
        end
    end
    
    %% Fit Copula model
%     [rhoHAT,nuHAT] = copulafit('t',[typeNXrsdU,typeNYrsdU],'Method', 'ApproximateML');
    
    %% Save the parameters
%     rsd.typen{iN}.rhoHAT = rhoHAT;
%     rsd.typen{iN}.nuHAT = nuHAT;
    rsd.typen{iN}.margHAT = phat.typen{iN};
    
end

%% Var n
nVarN = length(mdl.varn);
rsd.varn = cell(nVarN, 1);
phat.varn = cell(nVarN, 1);
for iN = 1:nVarN
    %% Collect type n residuals
    varNXrsd = cell2mat(cellfun(@(x)(x.MLERsd),mdl.varn{iN}.X,'UniformOutput',false)');
    varNYrsd = cell2mat(cellfun(@(x)(x.MLERsd),mdl.varn{iN}.Y,'UniformOutput',false)');
    nX = size(varNXrsd,2);
    nY = size(varNYrsd,2);
    
    %% Fit marginal models and transform residuals to uniform margins
    varNXrsdU = nan(size(varNXrsd));
    varNYrsdU = nan(size(varNYrsd));
    
    switch dens
        case 'TGN'
            phat.typen{iN} = zeros(nX + nY, 6);
        otherwise
            phat.typen{iN} = cell(nX + nY, 1);
    end
    
    for jX = 1:nX
        switch dens
            case 'TGN'
                p0 = fminsearch(@(p)(GenNormal_like(varNXrsd(:,jX),p(1),p(2),p(3))),[0,1,1]);
                phat.varn{iN}(jX,:) = fminsearch(@(p)(TwoSidedGenNormal_like(varNXrsd(:,jX),p(1),p(2:3),p(4:5),p(6))),[p0([1,2,2,3,3]),0.5]);
                varNXrsdU(:,jX) = TwoSidedGenNormal_CDF(varNXrsd(:,jX),phat.varn{iN}(jX,1),phat.varn{iN}(jX,2:3),phat.varn{iN}(jX,4:5),phat.varn{iN}(jX,6));
            otherwise
                XI = linspace(min(varNXrsd(:,jX))-0.5,max(varNXrsd(:,jX))+0.5,1000)';
                f = ksdensity(varNXrsd(:,jX), XI, 'Function', 'cdf');
                begin = find(f==0);
                if any(begin)
                    begin = begin(end);
                else
                    begin = 1;
                end
                einde = find(f==1);
                if any(einde)
                    einde = einde(1);
                else
                    einde = length(f);
                end
                phat.typen{iN}{jX} = [f(begin:einde), XI(begin:einde)];
                varNXrsdU(:,jX) = interp1(XI(begin:einde), f(begin:einde), varNXrsd(:,jX));
%                 varNXrsdU(:,jX) = ksdensity(varNXrsd(:,jX),varNXrsd(:,jX),'Function','cdf');
        end
    end
    
    for jY = 1:nY
        switch dens
            case 'TGN'
                p0 = fminsearch(@(p)(GenNormal_like(varNYrsd(:,jY),p(1),p(2),p(3))),[0,1,1]);
                phat.varn{iN}(nX+jY,:) = fminsearch(@(p)(TwoSidedGenNormal_like(varNYrsd(:,jY),p(1),p(2:3),p(4:5),p(6))),[p0([1,2,2,3,3]),0.5]);
                varNYrsdU(:,jY) = TwoSidedGenNormal_CDF(varNYrsd(:,jY),phat.varn{iN}(nX+jY,1),phat.varn{iN}(nX+jY,2:3),phat.varn{iN}(nX+jY,4:5),phat.varn{iN}(nX+jY,6));
            otherwise
                XI = linspace(min(varNYrsd(:,jY))-0.5,max(varNYrsd(:,jY))+0.5,1000)';
                f = ksdensity(varNYrsd(:,jY), XI, 'Function', 'cdf');
                begin = find(f==0);
                if any(begin)
                    begin = begin(end);
                else
                    begin = 1;
                end
                einde = find(f==1);
                if any(einde)
                    einde = einde(1);
                else
                    einde = length(f);
                end
                phat.typen{iN}{nX+jY} = [f(begin:einde), XI(begin:einde)];
                varNYrsdU(:,jY) = interp1(XI(begin:einde), f(begin:einde), varNYrsd(:,jY));
%                 varNYrsdU(:,jY) = ksdensity(varNYrsd(:,jY),varNYrsd(:,jY),'Function','cdf');
        end
    end
    
    %% Fit Copula model
%     [rhoHAT,nuHAT] = copulafit('t',[varNXrsdU,varNYrsdU],'Method', 'ApproximateML');
    
    %% Save the parameters
%     rsd.varn{iN}.rhoHAT = rhoHAT;
%     rsd.varn{iN}.nuHAT = nuHAT;
    rsd.varn{iN}.margHAT = phat.varn{iN};
    
end


end






