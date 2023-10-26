function out = transform2Org(val, drc, MM)

% Step 1: which bin?
A = BinAlc(drc, MM.Bn.Edg{1});

% Step 2: transform back to original margins
out = nan(size(val));
uVal = MM.CDF_Standard(val);

thr = MM.Thr(A,1);
nep = MM.NEP(1);

% above threshold
iExc = uVal > nep;
scl = MM.Scl(A,1);
shp = MM.Shp(1);
out(iExc) = GPD_iCDF((uVal(iExc) - nep)/(1 - nep), shp, thr(iExc), scl(iExc));

% below threshold
if any(~iExc)
    yA = MM.Y(MM.Bn.A==A);
    yA_below = yA(yA < thr);
    out(~iExc) = quickquantile(yA_below, uVal(~iExc) / nep);
end

% out = MM.INV(uVal,1,A);

end