function p0 = findMLE(like, pDmn, p0)

if nargin < 3
    p0 = [0.8;zeros(pDmn-1,1)];
    p1 = zeros(pDmn,1);
else
    p1 = zeros(size(p0));
end
while sum(abs(p0 - p1)) > 0.01
    p1 = p0;
    p0 = fminsearch(like,p0,optimset('Display','off'));
end

end