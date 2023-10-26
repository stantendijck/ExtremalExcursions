function data = getData(filenme)

load(filenme,'wam');

data.t = wam.t;

data.X = wam.hs + rand(size(wam.hs))*0.1 - 0.05;
data.Xdrc = mod(wam.wvd + rand(size(wam.wvd)),360);

data.Y = wam.ws + rand(size(wam.hs))*0.1 - 0.05;
data.Ydrc = mod(wam.wd + rand(size(wam.wd)),360);

end