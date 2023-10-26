function mdl = fitMdl1(data, respData, adpMCMCFlag, printOpts)

if nargin < 4
    printOn = false;
else
    printOn = printOpts{1};
end

maxLag = size(respData,2);

mdl = cell(maxLag, 1);
for i = 1:maxLag
    
    if printOn
        fprintf('Fit %s(%d) | X(1) > u\n',printOpts{2},printOpts{3}+i);
    end
    
    data.response = respData(:,i);
    
    keefProp = generateKeefProp(data.response, data.explanatory);
    
    like = @(d,p)(HT_like(d.response, d.explanatory, p, d.logExplanatory, keefProp));
    
    if adpMCMCFlag
        mdl{i} = adptMCMC(data, like, [0.8;0], 1e5, 0.001);
    else
        pDmn = 2;
        phat = findMLE(@(p)(like(data,p)), pDmn);
        
        mdl{i} = saveMLE(data, phat, like, {});
    end
end

end