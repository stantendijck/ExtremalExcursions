function data = readData()
%READDATA Summary of this function goes here
%   Detailed explanation goes here

load 'NORA10_6256N_0197E.mat' wam

data.X = wam.hs + rand(size(wam.hs))*0.1 - 0.05;
data.Y = wam.ws + rand(size(wam.hs))*0.1 - 0.05;

n = length(data.X);


marg_quant.X.upper = 0.8;
marg_quant.X.lower = 0.2;
marg_quant.Y.upper = 0.8;
marg_quant.Y.lower = 0.2;
data.unif = PIT(data, marg_quant);
data.exp.X = Exp_iCDF(data.unif.X);
data.exp.Y = Exp_iCDF(data.unif.Y);
data.lap.X = Laplace_iCDF(data.unif.X);
data.lap.Y = Laplace_iCDF(data.unif.Y);

% Plot
figure(1);clf;
plot( data.exp.X, 'k-');
hold on;
plot( data.exp.Y, 'r-');
xlim([0 100]);


q = 0.85;
data.thr = Laplace_iCDF(q);

% Initialisation
stormIn = data.lap.X(1) > data.thr | data.lap.Y(1) > data.thr;
if stormIn
    storm.data = [data.lap.X(1),data.lap.Y(1)];
    storm.index = 1;
else
    storm.data = nan(0,2);
    storm.index = [];
end

% Loop over the data
storms = {};
for i = 2:n
    if stormIn && (data.lap.X(i) > data.thr || data.lap.Y(i) > data.thr)
        storm.data(end+1,:) = [data.lap.X(i),data.lap.Y(i)];
        storm.index(end+1) = i;
    elseif stormIn
        storms{end+1} = storm;
        stormIn = false;
    elseif data.lap.X(i) > data.thr || data.lap.Y(i) > data.thr
        storm.data = [data.lap.X(i),data.lap.Y(i)];
        storm.index = i;
        stormIn = true;
    end
    if (data.lap.X(i) > data.thr || data.lap.Y(i) > data.thr) && i == n
        storms{end+1} = storm;
    end
end

figure(2); clf;
iCnt = 1;
iStorm = 0;
while iCnt <= 16
    iStorm = iStorm + 1;
    if length(storms{iStorm}.index)<10
        continue
    end
    subplot(4,4,iCnt);
    inds = storms{iStorm}.index;
    if storms{iStorm}.index(1) > 1
        inds = [inds(1)-1,inds];
    end
    if storms{iStorm}.index(end) < n
        inds = [inds,inds(end)+1];
    end
    plot(data.lap.X(inds));
    hold on;
    plot(data.lap.Y(inds));
    plot(xlim,data.thr*ones(2,1),'k--');
    iCnt = iCnt + 1;
end


data.storms = storms;

end

