function mdl = getAlphas(mdl, MLEs)

mdl.newMLE = nan(size(mdl.MLE));
mdl.newMLE([1,end]) = mdl.MLE([1,end]);

newfieldnames = cell(length(MLEs) - 1,1);
for i = 1:length(MLEs) - 1
    newfieldnames{i} = sprintf('alpha%d',i+1);
    
    if i == 1
        mdl.newMLE(i+1) = (MLEs{end}(1) - mdl.MLE(i))/MLEs{1}(1) + mdl.MLE(i+1);
    else
        %MLE1 = MLEs{1}(1);
        %MLE2 = MLEs{end}(1);
        
        %-mdl{i}.MLE(2) * mdlType1.X{1}.MLE(1) / mdlType1.Y{1}.MLE(1) + mdl{i}.MLE(3);
        mdl.newMLE(i+1) = -mdl.MLE(i) * MLEs{i-1}(1) / MLEs{i}(1) + mdl.MLE(i+1);
        %mdl.newMLE(i+1) = (MLE2 - mdl.MLE(1))/MLE1 + mdl.MLE(2);
    end
    

end


end