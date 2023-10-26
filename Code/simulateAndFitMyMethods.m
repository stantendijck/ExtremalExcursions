function [stormsOrg, mdl, mdldrc] = simulateAndFitMyMethods(trainData, opts, MM_HS, MM_WS, onlyFit)

if nargin < 5
    onlyFit = false;
end

fprintf("Fitting directional models\n");
mdldrc = fitDrcModel3(trainData);

%% Model Fit
mdl = mymodel(trainData, opts.maxLag+1);
if onlyFit
    stormsOrg = [];
    return;
end

%% Simulation on Laplace Margins
%Independent Sims
% storms = mysim(mdl, opts.nStormsSim, opts.methods, opts.simulationType, opts.rsdModel, trainData);

%Dependent Sims:
% For each sim, HS-Peak on Laplace margins; and upcrossing direction are the same
peaks = rand(opts.nStormsSim,1);

storms = mysim2(mdl, opts.nStormsSim, opts.methods, opts.simulationType, opts.rsdModel, trainData, peaks);

A = length(trainData.storms);
B = max(cellfun(@(x)(max(cellfun(@length,storms.(x)))),opts.methods));

fprintf("Simulating direction\n");
I = randsample(A, B, true);
stormsOrg = simulateDrcAll5(storms, mdldrc, trainData, MM_HS, MM_WS, I); % ARIMA for wave, but residual model for wind

% Historical Matching
fprintf('Simulating from historical  matching\n');
stormsOrg.HM = simulateHM2(opts.nStormsSim, trainData, MM_HS, mdl, peaks);

end