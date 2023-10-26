function [p, u, pconf] = chi(x, u)
%We need x to be on standard uniform margins, then this function 
%calculates chi(u) = P(X_1 > u | X_2 > u) for different values of u

if nargout == 3
    flag_conf = true;
else
    flag_conf = false;
end

n = size(x,1);
if nargin < 2
    max_u = 1 - 10/n;
    min_u = 0;
    u = linspace(min_u,max_u,1000);
end
p = zeros(size(u));
pconf = zeros(length(u),2);

% for i = 1:length(u)
% %     p(i) = mean(x(x(:,1)>u(i),2) > u(i));
%     p(i) = 2 - log(mean(x(:,1)<u(i) & x(:,2)<u(i)))/log(u(i));
%     if isnan(p(i))
%         p(i) = 0;
%     end
% end

for i = 1:length(u)
%     p(i) = mean(x(x(:,1)>u(i),2) > u(i));
    chi_i = mean(x(:,1)>u(i) & x(:,2)>u(i));
    p(i) = chi_i/(1-u(i));
    if isnan(p(i))
        p(i) = 0;
    elseif flag_conf
        if n*chi_i < 30
            binomconf = binoinv([0.025,0.975],n,chi_i);
        else
            if chi_i < 1 && chi_i > 0
                binomconf = norminv([0.025,0.975],n*chi_i,sqrt(n*chi_i*(1-chi_i)));
            else
                binomconf = [n*chi_i,n*chi_i];
            end
        end
        pconf(i,:) = (1/n*binomconf)/(1-u(i));
    end
end


end