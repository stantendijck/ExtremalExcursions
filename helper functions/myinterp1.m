function vq = myinterp1(x,v,xq)

[~,I] = unique(x);

vq = interp1(x(I),v(I),xq);

end