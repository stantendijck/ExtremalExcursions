function dataUnif = PIT(data, marg_quant)

%% PIT transform all data in data to uniform margins

nms = fieldnames(data);

% Remove time and directional components
remove_inds = [];
for iname = 1:length(nms)
    if strcmp(nms{iname},'t') || strcmp(nms{iname},'wd')
        remove_inds(end+1) = iname;
    end
end
nms(remove_inds)=[];

% Pre allocation
for iname = 1:length(nms)
    dataUnif.(nms{iname}) = zeros(length(data.(nms{iname})),1);
    
    dataUnif.prm.(nms{iname}).lower.quantile =  marg_quant.(nms{iname}).lower;
    dataUnif.prm.(nms{iname}).upper.quantile =  marg_quant.(nms{iname}).upper;
    
    % quantiles
    qExcLower = quantile(data.(nms{iname}), marg_quant.(nms{iname}).lower);
    qExcUpper = quantile(data.(nms{iname}), marg_quant.(nms{iname}).upper);
    
    % indices of exceedances
    iExcUpper = data.(nms{iname}) > qExcUpper;
    iExcLower = data.(nms{iname}) < qExcLower;
    iMiddle = ~iExcUpper & ~iExcLower;
    
    % fit GPD
    gpdprmUpper = fminsearch(@(p)(GPD_like(data.(nms{iname})(iExcUpper), qExcUpper, p(1), p(2))), [1, -0.001], optimset('Display','off'));
    gpdprmLower = fminsearch(@(p)(GPD_like(-data.(nms{iname})(iExcLower), -qExcLower, p(1), p(2))), [1, -0.001], optimset('Display','off'));
    gpdprm.(nms{iname}).upper = [qExcUpper, gpdprmUpper];
    gpdprm.(nms{iname}).lower = [qExcLower, gpdprmLower];
    
    dataUnif.prm.(nms{iname}).upper.xi = gpdprm.(nms{iname}).upper(3);
    dataUnif.prm.(nms{iname}).upper.mu = gpdprm.(nms{iname}).upper(1);
    dataUnif.prm.(nms{iname}).upper.sigma = gpdprm.(nms{iname}).upper(2);
    dataUnif.prm.(nms{iname}).lower.xi = gpdprm.(nms{iname}).lower(3);
    dataUnif.prm.(nms{iname}).lower.mu = gpdprm.(nms{iname}).lower(1);
    dataUnif.prm.(nms{iname}).lower.sigma = gpdprm.(nms{iname}).lower(2);

    
    % Transform bulk empirically
    qMiddleRange = marg_quant.(nms{iname}).upper - marg_quant.(nms{iname}).lower;
    dataUnif.(nms{iname})(iMiddle) = epit(data.(nms{iname})(iMiddle)) * qMiddleRange + marg_quant.(nms{iname}).lower;
    
    % Transform tails
    dataUnif.(nms{iname})(iExcUpper) = GPD_CDF(data.(nms{iname})(iExcUpper),gpdprm.(nms{iname}).upper(1),gpdprm.(nms{iname}).upper(2),gpdprm.(nms{iname}).upper(3)) * (1-marg_quant.(nms{iname}).upper) + marg_quant.(nms{iname}).upper;
    dataUnif.(nms{iname})(iExcLower) = (1 - GPD_CDF(-data.(nms{iname})(iExcLower),-gpdprm.(nms{iname}).lower(1),gpdprm.(nms{iname}).lower(2),gpdprm.(nms{iname}).lower(3))) * marg_quant.(nms{iname}).lower;
end

end
