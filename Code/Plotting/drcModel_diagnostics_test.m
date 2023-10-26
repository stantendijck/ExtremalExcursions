rl = linspace(0,30,1000);

a = -6; b = 6;
u = rand(1e5,2)*(b-a) + a;
s2 = nan(length(rl),1);
for ir = 1:length(rl)
    s2(ir) = std(atan(u(:,2)./(rl(ir)-u(:,1)))/(2*pi) * 360);
end

figure(24); clf;
plot(rl,s2)
hold on;
plot(xl(1:end-1),s,'r-')

%%
aY = min(data.Y([(1:end-1)',(2:end)']),[],2);
dY = data.Ydrc(1:end-1);
diffYdrc = mod(diff(data.Ydrc)+180,360)-180;
diffXdrc = mod(diff(data.Xdrc)+180,360)-180;
diffXY = mod(diffXdrc - diffYdrc+180,360)-180;
figure(25); clf;
scatter(aY,diffXY,5,dY,'filled');
colorbar;
figure(26); clf;
plot3(aY,diffXY,dY,'k.');
xlabel('WS'); ylabel('Drc(WS)-Drc(HS)'); zlabel('WS-Drc');
figure(27); clf;
for i = 1:36
    I = dY > (i-1)*10 & dY < i*10;
    subplot(6,6,i);
    plot(aY(I), diffXY(I),'k.');
    xlim([0 35]);
    ylim([-200 200])
    grid on;
end
% subplot(1,2,1);
% plot(aY,diffYdrc,'k.')
% subplot(1,2,2);
% plot(aY,diffXY,'r.')

% 

xl = linspace(min(aY),max(aY),100);
s = nan(length(xl)-1,1);
for i = 1:length(xl)-1
    I = aY >= xl(i) & aY < xl(i+1);
    s(i) = std(diffXY(I));
end

a = (max(aY)-min(aY))./(aY-min(aY));
inds = min(100,max(1,round(1./(1+a) * 100)));
% xl(inds)

% figure(26); clf;
% plot(aY,diffXY./s(inds),'k.');
%%
figure;
plot(xl(1:end-1),s,'r-')

%%
n = length(data.Y);
U = rand(n,1);
% E = -log(1-U);

% f = inverse of the variance
% f goes from 100 to 25 from 0 to 5; then from (5,25) to (22,13)
f = @(e)((100 - e * 75/5) .* (e >= 0 & e < 5) + (485/17 - 12/17 * e) .* (e >= 5));
% c = 0.05;
% f = @(e)(e./(e+1) * c); %e=0: variance = infinity; e = infinity: variance = 1



Z1 = trnd(7,n,1).*(f(data.Y));
Z2 = rand(n,1)*360 - 180;

Z = Z1;
I = rand(n,1) < 0.01;
Z(I) = Z2(I);

figure(27); clf;
plot(data.Y,mod(Z+180,360)-180,'k.');




