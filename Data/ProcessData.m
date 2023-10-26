function data = ProcessData(data, q)
% Get storms

% Storm threshold
if nargin < 2
    q = 0.95;
end

% Defaults
n = length(data.lap.X);
data.thr = Laplace_iCDF(q);

% Initialisation
% stormIn = data.lap.X(1) > data.thr | data.lap.Y(1) > data.thr;
stormIn = data.lap.X(1) > data.thr;
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
    %if stormIn && (data.lap.X(i) > data.thr || data.lap.Y(i) > data.thr)
    if stormIn && data.lap.X(i) > data.thr
        storm.data(end+1,:) = [data.lap.X(i),data.lap.Y(i)];
        storm.index(end+1) = i;
    elseif stormIn
        storms{end+1,1} = storm;
        stormIn = false;
    %elseif data.lap.X(i) > data.thr || data.lap.Y(i) > data.thr
    elseif data.lap.X(i) > data.thr
        storm.data = [data.lap.X(i),data.lap.Y(i)];
        storm.index = i;
        stormIn = true;
    end
    %if (data.lap.X(i) > data.thr || data.lap.Y(i) > data.thr) && i == n
    if data.lap.X(i) > data.thr && i == n
        storms{end+1,1} = storm;
    end
end

figure(7); clf;
iCnt = 1;
iStorm = 0;
while iCnt <= 16 && iStorm + 1 <= length(storms)
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

nStorms = length(data.storms);
for iStorm = 1:nStorms
    [~,I] = max(data.storms{iStorm}.data(:,1));
    data.storms{iStorm}.peakInd = data.storms{iStorm}.index(I);
    data.storms{iStorm}.risingInd = data.storms{iStorm}.index(1:I);
    data.storms{iStorm}.fallingInd = data.storms{iStorm}.index(I:end);
end






end