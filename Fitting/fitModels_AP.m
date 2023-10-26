function mdl = fitModels_AP(htdata, printOn)

if nargin < 2
    printOn = true;
end
adpMCMCFlag = false;

maxLag = ceil(length(htdata{1}.X)/2);

varNms = {'X','Y'};

data = cell(maxLag,1);
for iLag = 1:maxLag
    if iLag == 4
        debug = true;
    end
    fprintf('%d\n',iLag');
    a = cellfun(@(x)([x.X(max(1,maxLag-iLag+1):min(end,maxLag+iLag-1));x.Y(max(1,maxLag-iLag+1):min(end,maxLag+iLag-1))]),htdata,'UniformOutput',false)';
    if length(a{end}) ~= length(a{end-1})
        data{iLag} = cell2mat(a(1:end-1))';
    else
        L = cellfun(@length,a);
        mL = median(L);
        data{iLag} = cell2mat(a(L==mL))';
    end
        
    
    d = struct();
    d.explanatory = data{iLag}(:,iLag);
    d.logExplanatory = log(d.explanatory);
    d.response = data{iLag}(:,[1:iLag-1,iLag+1:end]);
    phat = cell(size(d.response,2),1);
    for k = 1:size(d.response,2)
        keefProp = generateKeefProp(d.response(:,k),d.explanatory);
        like = @(d,p)(HT_like(d.response(:,k), d.explanatory, p, d.logExplanatory, keefProp));
        pDmn = 2;
        phat{k} = findMLE(@(p)(like(d, p)), pDmn);
    end
    mdl{iLag}.phat = phat;
    mdl{iLag}.data = d;
end

end