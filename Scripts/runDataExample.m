%% Code for fitting models and simulating on original margins
clear;
clc;
close all;
set(0,'DefaultFigureWindowStyle','docked')
% %% Add Paths
addpath(genpath('../Code'));

%% User settings
opts = struct;
% Choose which models to fit
opts.maxLag = 5;                      % maxLag >= 1 to consider. NB: maxLag = 5 => models up to order 4
opts.methods = {'EVAR','VAR','MMEM'}; % pick from 'EVAR', 'VAR' (= EVAR0), 'MMEM'

% Choose between fitting model on all data or running a simulation study
% switch by setting: runCase = "SimulationStudy" or "All"
% 
% In simulation study, you have to set additionally:
% opts.nStormsSim = 20000;
% nBootstraps = 50;
% NB: this takes a very long time
%
% In running the model just once on the data without bootstrapping, set instead:
% opts.nStormsSim = 2000;
% nBootstraps = 1;
% This should run relatively quickly.
opts.nStormsSim = 0; 
nBootstraps = 1;
runCase = "All";                      % Type: "All" or "SimulationStudy"

if ~exist('Figures','dir')
    mkdir('Figures');
end


%% Default settings
parOn = false;
nCores = 4;
if parOn && isempty(gcp('nocreate'))
    h = parcluster;
    h.NumWorkers = nCores;
    parpool(h.NumWorkers);
end

%% Load data
filenme = '../Data/data.mat';
data = getData(filenme);

%% Marginal Models: HS
MM_HS = marginal(data.X, data.Xdrc);
[data.lap.X, ~] = MM_HS.Margins(1);

%% Marginal Model: WS
MM_WS = marginal(data.Y, data.Ydrc);
[data.lap.Y, ~] = MM_WS.Margins(1);

%% Storm data
q = 0.95;
data = ProcessData(data, q);
nStorms = length(data.storms);

%% Split up data into train and test
rng(123456);

% Defaults
opts.simulationType = 'peak';         % simulate from peak onwards, no other options supported
opts.rsdModel = 'none';               % 'none': empirical copula for residuals, no other options supported

% Pre allocation
SCV = cell(nBootstraps, 1);

impact=cell(nBootstraps, 1);

tic();
% Folds
myMdlFits = cell(nBootstraps,1);
mdldrc = cell(nBootstraps,1);
for iBt = 1:nBootstraps
    fprintf('Number of bootstraps %d/%d\n',iBt,nBootstraps);
    
    switch runCase
        case "All"
            newInds = randsample(nStorms,nStorms,true);
        otherwise
            newInds = randsample(nStorms,nStorms);
    end
    if iBt == 1
        newInds = (1:nStorms)';
    end
    
    switch runCase
        case "SimulationStudy"
            q = round(0.25*nStorms);
            trainInds = newInds(1:q);
            testInds = newInds(q+1:end);
        case "All"
            trainInds = newInds;
            testInds = newInds;
    end

    trainData = data; trainData.storms = data.storms(trainInds);
    testData = data; testData.storms = data.storms(testInds);
    
    %% Fitting and simulating from MMEM, EVAR and HM
    [stormsOrg, myMdlFits{iBt}, mdldrc{iBt}] = simulateAndFitMyMethods(trainData, opts, MM_HS, MM_WS);
    
    %% Output diagnostics
    impact{iBt} = outSCV4(testData, stormsOrg, trainData);

end
toc();
%%

% nB = 200, Rsd models turned off - used for generateTablesForPaper
% save('ParameterEstimatesWorkspace.mat','-v7.3')

%%
saveOn = false;
mdl = myMdlFits{1};         % pick first bootstrap for generating plots for paper
generatePlotsForPaper;

%%
generatePlotsFromSimulationStudy;

%% Generate Tables: uncertainty = Bootstrap uncertainty
s = 'R'; % change from 'R' to 'F' to go from PrePeak to PostPeak
generateTablesForPaper;

%% Plot diagnostics model fit

% Residuals vs Explanatory variables
figure(1); clf;
for i = 1:8
    subplot(3,3,i);
    expl = myMdlFits{1}.R.fit.typen{4}.X{1}.data.explanatory(:,i);
    rsd = myMdlFits{1}.R.fit.typen{4}.X{1}.MLERsd;
    plot(expl,rsd,'k.')
end

% Correlation matrix
figure(2); clf;
imagesc(corr([myMdlFits{1}.R.fit.typen{4}.X{1}.data.explanatory,myMdlFits{1}.R.fit.typen{4}.X{1}.MLERsd]))


