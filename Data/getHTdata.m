function htdata = getHTdata(data,maxLag,type)

if nargin < 3
    translate = @(inds,r)([false(r,1);inds(1:end-r)]);
    I = data.lap.X > data.thr;
    I(end-maxLag+1:end) = false;
    htdata.X.data = nan(sum(I),maxLag+1);
    htdata.Y.data = nan(sum(I),maxLag+1);
    for i = 0:maxLag
        htdata.X.data(:,i+1) = data.lap.X(translate(I,i));
        htdata.Y.data(:,i+1) = data.lap.Y(translate(I,i));
%         htdata.Y_unifDrc.data(:,i+1) = data.lap.Y_unifDrc(translate(I,i));
    end
    htdata.logX = log(data.lap.X(I));
else
    if strcmp(type,'aroundpeak')
        peakInds = cellfun(@(x)(x.peakInd),data.storms);
        m = maxLag - 1;
        htdata = cell(length(peakInds),1);
        for iPk = 1:length(peakInds)
            htdata{iPk}.X = data.lap.X(max(1,-m+peakInds(iPk)):min(end,m+peakInds(iPk)));
            htdata{iPk}.Y = data.lap.Y(max(1,-m+peakInds(iPk)):min(end,m+peakInds(iPk)));
            htdata{iPk}.logX = log(htdata{iPk}.X);
        end
    else
        translate = @(inds,r)(inds+r);
        factor = 1;
        switch type
            case 'falling'
                I = cell2mat(cellfun(@(x)(x.fallingInd'),data.storms,'UniformOutput',false));
%                 I = cell2mat(cellfun(@(x)(x.fallingInd),data.storms,'UniformOutput',false)');
            case 'rising'
                I = cell2mat(cellfun(@(x)((x.risingInd(1):(x.risingInd(end)))'),data.storms,'UniformOutput',false));
%                 I = cell2mat(cellfun(@(x)(x.risingInd(1):(x.risingInd(end))),data.storms,'UniformOutput',false)');
                factor = -1;
        end
        I(I>length(data.lap.X)-maxLag+1) = [];
        I(I<=maxLag) = [];
        
        inds = translate(I, factor * maxLag)>0 & translate(I, factor * maxLag)<=length(data.lap.X);
        I = I(inds);
        htdata.X.data = nan(length(I),maxLag+1);
        htdata.Y.data = nan(length(I),maxLag+1);
        for i = 0:maxLag
            htdata.X.data(:,i+1) = data.lap.X(translate(I, factor * i));
            htdata.Y.data(:,i+1) = data.lap.Y(translate(I, factor * i));
    %         htdata.Y_unifDrc.data(:,i+1) = data.lap.Y_unifDrc(translate(I, factor * i));
        end
        htdata.logX = log(data.lap.X(I));
    end
end
% translate = @(inds,r)(inds+r);
% I = cellfun(@(x)(x.index(1)),data.storms);
% I(I>length(data.lap.X)-maxLag+1) = [];
% htdata.X.data = nan(length(I),maxLag+1);
% htdata.Y.data = nan(length(I),maxLag+1);
% for i = 0:maxLag
%     htdata.X.data(:,i+1) = data.lap.X(translate(I,i));
%     htdata.Y.data(:,i+1) = data.lap.Y(translate(I,i));
%     htdata.Y_unifDrc.data(:,i+1) = data.lap.Y_unifDrc(translate(I,i));
% end
% htdata.logX = log(data.lap.X(I));
end
