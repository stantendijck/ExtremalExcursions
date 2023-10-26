function A=BinAlc(X,Edg)
%Allocate data to covariate bins
%Input:
%   X      n x 1 covariate directions (on 0,360)
%   Edg    nBin x 1  Bin edges [0,360]; bins will wrap around 0 degrees
% Output:
%   A         n x nDmn bin allocation
  


if any(X>360) || any(X<0)
    error('X expected on range [0,360]');
end
if any(Edg>360) || any(Edg<0)
    error('Edg expected on range [0,360]');
end
MdEdg=mod(Edg,360);
if numel(unique(MdEdg)) < numel(Edg)
    error('You have supplied 2 bin edges which correspond to same point on a circle. For a single bin, please provide 1D DrcEdg (any value will do). For 2 bins, please provide 2D DrcEdg containing values which differ on 360deg circle.');
end

%check Edg sorted
Edg=sort(Edg);
if size(Edg,2) > size(Edg,1)
    error('DrcEdg supplied is a row vector, please provide as col vector')
end
nBin=numel(Edg);

%find indices of directions X which fall in each directional bin
ExtEdg=unique([0;Edg;360]); %pad with 0, 360 if needed  %added..***
A=discretize(X,ExtEdg);
A(A==nBin+1)=1; %wrap bins around 0
