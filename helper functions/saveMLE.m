function mdl = saveMLE(data, phat, like, MLES)

mdl.data = data;
mdl.MLE = phat;
if nargin >= 4
    mdl = getAlphas(mdl, MLES);
else
    mdl = getAlphas(mdl);
end
mdl.like = like(data, phat);
mdl.MLERsd = getRsdHT(mdl);

end