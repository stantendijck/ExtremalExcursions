function rsd = getRsdHT(mdl)
% Get the residuals of the fitted HT model
%
% INPUT:
% mdl      struct, fit of HT model
%
% OUTPUT
% MLERsd    Residuals of the fitted HT model using the MLE
 
rsd = (mdl.data.response - mdl.data.explanatory * mdl.newMLE(1:end-1)) ./ exp(mdl.data.logExplanatory * mdl.newMLE(end));

end