function mdl = fitMdln_v2(data, respData, mdlType1, p0, respNme, functionNme, adpMCMCFlag, printOpts)

if nargin < 4
    printOn = false;
else
    printOn = printOpts{1};
end

maxLag = size(respData,2);

mdl = cell(maxLag, 1);
for i = 1:maxLag
    if printOn
        str = '(X(1) > u, Y(1)';
        for itype = 3:functionNme
            str = strcat(str, sprintf(', X(%d), Y(%d)',itype-1,itype-1));
        end
        fprintf('Fit %s(%d) | %s)\n', printOpts{2}, printOpts{3}+i+1,str);
    end
    
    data.response = respData(:,i);
    
    keefProp = generateKeefProp(data.response, data.explanatory);
    
    ind = i + functionNme - 2;
    if strcmp(respNme,'Y')
        ind = ind + 1;
    end
    
    pDmn = 2 * functionNme - 1;
    
    MLEs = cell(2 * (functionNme - 1), 1);
    for j = 1:functionNme - 2
        MLEs{j} = mdlType1.X{j}.MLE;
    end
    for j = 1:functionNme - 1
        MLEs{functionNme-2+j} = mdlType1.Y{j}.MLE;
    end
    MLEs{end} = mdlType1.(respNme){ind}.MLE;
    
    like = @(d,p)(HTn_like(d.response, d.explanatory, functionNme, p, d.logExplanatory, MLEs, keefProp));

    phat = findMLE(@(p)(like(data,p)), pDmn, p0{i}.MLE);
    if adpMCMCFlag
        mdl{i} = adptMCMC(data, like, phat, 5e4, 0.0005);
    else
        mdl{i} = saveMLE(data, phat, like, MLEs);
    end
    
    mdl{i} = getAlphas(mdl{i}, MLEs);

end