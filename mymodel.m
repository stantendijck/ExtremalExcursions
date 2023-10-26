function mdl = mymodel(data, maxLag)

mdl = struct();

%% storm peak
mdl.SP = fitModelPeak(data);

%% around peak:
mdl.AP.data = getHTdata(data, maxLag, 'aroundpeak');
mdl.AP.fit = fitModels_AP(mdl.AP.data);
mdl.AP.fit = fitRsdModels2(mdl.AP.fit, true);

%% Rising
mdl.R.data = getHTdata(data, maxLag, 'rising');
mdl.R.fit = fitModels_v3(mdl.R.data, false);
mdl.R.fit.rsd = fitRsdModels(mdl.R.fit, true, 'KD');

%% Falling
mdl.F.data = getHTdata(data, maxLag, 'falling');
mdl.F.fit = fitModels_v3(mdl.F.data, false);
mdl.F.fit.rsd = fitRsdModels(mdl.F.fit, true, 'KD');


end