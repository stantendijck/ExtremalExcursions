function [HSL,WSL]=MarginsHM(currStorm, MM)
%[YGmbl,YUnif]=Margins(MM)



MM_HS = MM(1);
MM_WS = MM(2);

%% Model Non Exceedances:
cdf = cell(MM_HS.Bn.nBin,1);

IExc = MM_WS.Y > MM_WS.Thr(MM_WS.Bn.A,1);
dataNonExcHS = MM_HS.Y(~IExc);
for i = 1:MM_HS.Bn.nBin
    binI = MM_HS.Bn.A(~IExc) == i;
    y = dataNonExcHS(binI);
    [f,x] = ecdf(y);
    f(1) = []; x(1) = [];
    cdf{i}.HS = griddedInterpolant(x,f);
end

IExc = MM_WS.Y > MM_WS.Thr(MM_WS.Bn.A,1);
dataNonExcWS = MM_WS.Y(~IExc);
for i = 1:MM_WS.Bn.nBin
    binI = MM_WS.Bn.A(~IExc) == i;
    y = dataNonExcWS(binI);
    [f,x] = ecdf(y);
    f(1) = []; x(1) = [];
    cdf{i}.WS = griddedInterpolant(x,f);
end

HS = currStorm.originalMargins(:,1);
WS = currStorm.originalMargins(:,2);

HSdrc = currStorm.direction(:,1);
WSdrc = currStorm.direction(:,2);

n = length(HS);

% first find out if above QR threshold

% first find which knots are closest
HSL = zeros(n,1);
WSL = zeros(n,1);
for i = 1:n
    tau = MM_HS.NEP(1);
    I = BinAlc(HSdrc(i),MM_HS.Bn.Edg{1});
    threshold = MM_HS.Thr(I,1);
    if HS(i) > threshold % use exceedance model
        % get gp-parameters
        xi = MM_HS.Shp(1);
        sigma = MM_HS.Scl(I,1);
        
        
        HSL(i) = Laplace_iCDF(tau + (1-tau) * GPD_CDF(HS(i),threshold,sigma,xi));
    else
        HSL(i) = Laplace_iCDF(tau * cdf{I}.HS(HS(i)));
    end
    
    tau = MM_WS.NEP(1);
    I = BinAlc(WSdrc(i),MM_WS.Bn.Edg{1});
    threshold = MM_WS.Thr(I,1);
    if WS(i) > threshold % use exceedance model
        % get gp-parameter
        xi = MM_WS.Shp(1);
        sigma = MM_WS.Scl(I,1);
        
        
        WSL(i) = Laplace_iCDF(tau + (1-tau) * GPD_CDF(WS(i),threshold,sigma,xi));
    else
        WSL(i) = Laplace_iCDF(tau * cdf{I}.WS(WS(i)));
    end
end
end
