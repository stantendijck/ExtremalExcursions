function MM = marginal(X, Xdrc)

MM = struct();

%% Quantile regression


% Bn=CovariateBinning(X,Edg,IsPrd,Y,RspLbl,CvrLbl);
nBin = 20;
Edg = {(0:(360/nBin):359.99)'};
% Edg = [Edg(1:end-1),Edg(2:end)];
% Edg = mat2cell(Edg,ones(11,1),2);

IsPrd = 1;


Dat.X = Xdrc;
Dat.Y = X;
Dat.IsPrd = IsPrd;
Dat.CvrLbl = {'Direction'};
Dat.RspLbl = {'Hs'};

Bn=CovariateBinning(Dat.X, Edg, Dat.IsPrd, Dat.Y, Dat.RspLbl, Dat.CvrLbl);

iDmn = 1;
NEP = 0.9;
MM = MarginalModel(Dat,iDmn,NEP,Bn);


% Bss=Basis("Direction",1,16,32,3,1);
% Bss = Bss.MakeAll(Xdrc);
% [Dmn, BssOpts] = Bss.to_Bars(Bss.A, "Drc");
% 
% Dst = Distribution.Laplace(X, Dmn, BssOpts(1));
% 
% Dst.Tau = 0.9;
% Dst = Dst.maxLike;
% 
% %% Split data into exceedances and non-exceedances
% 
% % [~,ind] = min(abs(bss-Dat.X'));
% IExc = X > Dst.Prm.Prd(Dst.Dmn.A);
% 
% figure; clf;
% plot(Xdrc(~IExc),X(~IExc),'.','Color',[0.6,0.6,0.6]);
% hold on;
% plot(Xdrc(IExc),X(IExc),'k.');
% plot(Dst.Dmn.X, Dst.Prm.Prd,'r-','LineWidth',2);
% title('Quantile Regression');
% ylabel('H_S');
% xlabel('Direction');
% 
% 
% %% model exceedances
% BssExc = Basis("Direction",1,16,32,3,1);
% BssExc = BssExc.MakeAll(Xdrc(IExc));
% [DmnExc, BssOptsExc] = BssExc.to_Bars(BssExc.A, "Drc");
% 
% dataExc = X(IExc) - Dst.Prm.Prd(DmnExc.A);
% 
% DistGP = Distribution.GeneralisedParetoOrtho(dataExc, DmnExc, BssOptsExc);
% 
% DistGP.checkCffFlag = false;
% 
% % Cross validation set up
% CV.nCV=10;
% CV.nCVreps=1;
% 
% Kpp = 0.05;
% LmbSet = 10.^(linspace(-3,3,10)');
% CV.SmthSet=cat(3, LmbSet, LmbSet*Kpp); 
% 
% DistGP = DistGP.maxLikeCV(CV);
% 
% 
% %%
% figure; clf;
% subplot(1,2,1);
% plot(Xdrc(IExc), dataExc,'k.');
% subplot(2,2,2);
% plot(DistGP.Dmn.X, DistGP.Prm(1).Prd,'r-'); title(DistGP.Prm(1).Lbl); hold on;
% subplot(2,2,4);
% plot(DistGP.Dmn.X, DistGP.Prm(2).Prd,'r-'); title(DistGP.Prm(2).Lbl); hold on;
% 
% 
% %% Transform exceedances to Laplace margins
% mu = 0;
% xi = DistGP.Prm(1).Prd(DistGP.Dmn.A);
% sig = DistGP.Prm(2).Prd(DistGP.Dmn.A) ./ (1 + xi);
% newDataExc = Laplace_iCDF(GPD_CDF(dataExc, mu, sig, xi)*(1-Dst.Tau) + Dst.Tau);
% 
% 
% %%
% figure; clf;
% subplot(1,2,1);
% plot(Xdrc(IExc), dataExc,'k.');
% subplot(1,2,2);
% plot(Xdrc(IExc), newDataExc,'k.');
% 
% %% non Exceedances
% dataNonExc = X(~IExc);
% dataNonExcDrc = Xdrc(~IExc);
% 
% MdlNonExc = cell(Dst.Dmn.N_Bin,1);
% for i = 1:Dst.Dmn.N_Bin
%     binI = Dst.Dmn.A(~IExc) == i;
%     y = dataNonExc(binI);
%     [f,x] = ecdf(y);
%     f(1) = []; x(1) = [];
%     icdf = griddedInterpolant(f, x);
%     cdf = griddedInterpolant(x,f);
%     MdlNonExc{i}.icdf = icdf;
%     MdlNonExc{i}.cdf = cdf;
%     MdlNonExc{i}.margins.Y = Laplace_iCDF(cdf(y) * Dst.Tau - eps);
%     MdlNonExc{i}.margins.X = dataNonExcDrc(binI);
% end
% 
% 
% %% Laplace margins
% orgData = nan(size(X));
% orgData(IExc) = newDataExc;
% for i = 1:Dst.Dmn.N_Bin
%     orgData(~IExc&Dst.Dmn.A==i) = MdlNonExc{i}.margins.Y;
% end
% figure(3);clf; plot(Xdrc,orgData,'k.');
% title('Data on standard Laplace margins');
% xlabel('Direction');
% 
% 
% %%
% MM.DrcX = Xdrc;
% MM.Y = X;
% MM.margins = orgData;
% 
% MM.QR.mdl = Dst;
% 
% MM.Exc = struct();
% MM.Exc.mdlname = 'GeneralisedPareto';
% MM.Exc.mdl = DistGP;
% 
% MM.NonExc = struct();
% MM.NonExc.mdlname = 'Empirical';
% MM.NonExc.mdl = MdlNonExc;


