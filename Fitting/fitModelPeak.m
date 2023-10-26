function mdl = fitModelPeak(data)


iPeak = cellfun(@(x)(x.peakInd),data.storms);
peakData = [data.lap.X(iPeak),data.lap.Y(iPeak)];

opts = optimset('MaxFunEvals',3000,'MaxIter',3000);
q =0.7:0.01:0.9;
phat = cell(length(q),1);
for iq = 1:length(q)
    m(iq) = quantile(peakData(:,1),q(iq));
    iExc = peakData(:,1)>m(iq);
    phat{iq} = fminsearch(@(p)(GPD_like(peakData(iExc,1),m(iq),p(1),p(2))),[1,0.1],opts);
end
mdl = struct();
mdl.X.q = q(21);
mdl.X.MLE = [m(21), phat{21}];
mdl.X.data = peakData(:,1);

dat.explanatory = peakData(:,1);
dat.logExplanatory = log(dat.explanatory);

% Fit Y(i) | X(1) > u                                   (i = 1, 2, ...)
mdl.Y = fitMdl1(dat, peakData(:,2), false, {1, 'Y', 0});





end