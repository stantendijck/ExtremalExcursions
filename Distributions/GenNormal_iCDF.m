function F = GenNormal_iCDF(p, m, a, b)

if a < 0 || b < 0
    F = nan;
    return
end

F = sign(p-0.5) .* gaminv(2*abs(p-0.5),1/b,a^b).^(1/b) + m;

end
