function kcond = keefConstraints(X, Y, p, q, X_beta, z, keefProp)

if size(q,1) < size(q,2)
    q = q';
end

v = keefProp.X.max + 1;
nDmn = size(Y,2);

alpha = p(1,:); beta = p(2,:);

if any(abs(alpha) > 1) || any(beta > 1)
    kcond = zeros(size(q));
end

if q == 1
    z_plus_q = keefProp.plus.max;
    z_min_q = keefProp.min.max;
    z_q = max(z);
elseif q == 0
    z_plus_q = keefProp.plus.min;
    z_min_q = keefProp.min.min;
    z_q = min(z);
elseif all(q == [0;1])
    z_plus_q = [keefProp.plus.min;keefProp.plus.max];
    z_min_q = [keefProp.min.min;keefProp.min.max];
    z_q = [min(z);max(z)];
else
    z_plus_q = quantile(keefProp.plus.data,q);
    z_min_q = quantile(keefProp.min.data,q);
    z_q = quantile(z,q);
end

min_vec1 = ones(length(q),3);
min_vec2 = ones(length(q),3);

kcond = 1;
for id = 1:nDmn
    min_vec1(:,2) = 1 - beta(id)*v^(beta(id)-1)*z_q(:,id);
    min_vec1(:,3) = 1 - v^(beta(id)-1)*z_q(:,id) + z_plus_q(:,id)/v;
    
    a_beta = (1 - alpha(id))^(-beta(id)/(1-beta(id)));
    bz_beta = exp(log(beta(id)*z_q(:,id))/(1-beta(id)));
    pos_quant1 = (1 - 1/beta(id)) * a_beta * bz_beta + z_plus_q(:,id);
    
    cond1_flag1 = alpha(id) <= min(min_vec1,[],2);
    cond1_flag2 = min_vec1(:,2) < alpha(id) & pos_quant1 > 0;
    
    cond1 = cond1_flag1 | cond1_flag2;
    
    min_vec2(:,2) = 1 + beta(id)*v^(beta(id)-1)*z_q(:,id);
    min_vec2(:,3) = 1 + v^(beta(id)-1)*z_q(:,id) - z_min_q(:,id)/v;
    
    a_beta = (1 + alpha(id))^(-beta(id)/(1-beta(id)));
    bz_beta = exp(log(-beta(id)*z_q(:,id))/(1-beta(id)));
    pos_quant2 = (1 - 1/beta(id)) * a_beta * bz_beta - z_min_q(:,id);
    
    cond2_flag1 = -alpha(id) <= min(min_vec2,[],2);
    cond2_flag2 = min_vec2(:,2) < -alpha(id) & pos_quant2 > 0;
    
    cond2 = cond2_flag1 | cond2_flag2;
    
    kcond = kcond & (cond1 & cond2);
end


end % KeefConstraints