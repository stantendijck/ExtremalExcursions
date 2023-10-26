function NLOGL = negloglike(rsd, margins, copula, index)

%% Calculate MLEs of marginals
[n,d] = size(rsd);
switch margins
    case 'Gaussian'
        MLE.m = mean(rsd)';
        MLE.s = cov(rsd);
        
    case 'GeneralisedGaussian'
        p0 = [0,1,1];

        MLE.p = cell(d,1);
        for i = 1:d
            fun = @(p)(GenNormal_like(rsd(:,i), p(1), p(2), p(3)));
            MLE.p{i} = num2cell(fminsearch(fun, p0));
        end
        
    case 'TwoSidedGeneralisedGaussian'
        p00 = [0,1,1];

        MLE.p = cell(d,1);
        for i = 1:d
            fun0 = @(p)(GenNormal_like(rsd(:,i), p(1), p(2), p(3)));
            phat0 = fminsearch(fun0, p00);
            p0 = [phat0([1,2,2,3,3]),0.5];
            
            fun = @(p)(TwoSidedGenNormal_like(rsd(:,i), p(1), p(2:3), p(4:5), p(6)));
            opts = optimset('MaxIter',500*6,'MaxFunEvals',500*6);
            MLE.p{i} = num2cell(fminsearch(fun, p0, opts));
        end
end

switch margins
    case 'GeneralisedGaussian'
        pdf = @(x,p1,p2,p3)(GenNormal_PDF(x,p1,p2,p3));
        cdf = @(x,p1,p2,p3)(GenNormal_CDF(x,p1,p2,p3));
        
    case 'TwoSidedGeneralisedGaussian'
        pdf = @(x,p1,p2,p3,p4,p5,p6)(TwoSidedGenNormal_PDF(x,p1,[p2,p3],[p4,p5],p6));
        cdf = @(x,p1,p2,p3,p4,p5,p6)(TwoSidedGenNormal_CDF(x,p1,[p2,p3],[p4,p5],p6));
end

%% if nargin < 3: no copula
if nargin < 3
    switch margins
        case 'Gaussian'
            NLOGL = -sum(log(normpdf(rsd, MLE.m ,sqrt(MLE.s))));
%             NLOGL = std(rsd);
            
        case {'GeneralisedGaussian','TwoSidedGeneralisedGaussian'}
            NLOGL = -sum(log(pdf(rsd, MLE.p{1}{:})));
%             NLOGL = std(rsd);
    end
    return;
else
    negInds = (1:d)';
    negInds(index) = [];
    
    % Check if margins is Gaussian or other
    switch margins
        case 'Gaussian'
            % Compute NegLogLike for Independent and Gaussian copulas
            switch copula
                case 'independent'
                    NLOGL = -sum(log(normpdf(rsd(:,index), MLE.m(index), sqrt(MLE.s(index,index)))));
%                     NLOGL = std(rsd(:,index));
                    
                case 'Gaussian'
                    [condM,condS] = conditionalGaussian(MLE.m, MLE.s, index);
                    cM_vec = cellfun(condM,num2cell(rsd(:,negInds),2));
                    NLOGL = -sum(log(normpdf(rsd(:,index), cM_vec, sqrt(condS))));
%                     NLOGL = std(rsd(:,index)-cM_vec);
            end
            
        otherwise
            % Compute NegLogLike for independent and Gaussian copulas
            switch copula
                case 'independent'
                    NLOGL = -sum(log(pdf(rsd(:,index), MLE.p{index}{:})));
%                     NLOGL = std(rsd(:,index));
                    
                case 'Gaussian'
                    rsdNormal = nan(size(rsd));
                    for i = 1:d
                        rsdNormal(:,i) = norminv(cdf(rsd(:,i),MLE.p{i}{:}));
                    end

                    m = mean(rsdNormal)';
                    s = cov(rsdNormal);

                    [condM, condS] = conditionalGaussian(m, s, index);
                    cM_vec = cellfun(condM, num2cell(rsdNormal(:,negInds),2));
                    
                    % Method 1
                    %LL = zeros(5,1);
                    %LL(1) = - n / 2 * log(2*pi);
                    %LL(2) = - n / 2 * log(condS);
                    %LL(3) = - sum((rsdNormal(:,index) - cM_vec).^2 / (2*condS));
                    %LL(4) = - sum(log(normpdf(rsdNormal(:,index))));
                    %LL(5) = sum(log(pdf(rsd(:,index), MLE.p{index}{:})));
                    %NLOGL = -sum(LL);
                    
                    % Method 2
                    p00 = [0,1,1];
                    fun0 = @(p)(GenNormal_like(rsdNormal(:,index)-cM_vec, p(1), p(2), p(3)));
                    phat0 = fminsearch(fun0, p00);
                    p0 = [phat0([1,2,2,3,3]),0.5];

                    fun = @(p)(TwoSidedGenNormal_like(rsdNormal(:,index)-cM_vec, p(1), p(2:3), p(4:5), p(6)));
                    opts = optimset('MaxIter',500*6,'MaxFunEvals',500*6);
                    phat = fminsearch(fun, p0, opts);

                    LL = zeros(3,1);
                    LL(1) = fun(phat);
                    LL(2) = -sum(log(normpdf(rsdNormal(:,index))));
                    LL(3) = sum(log(pdf(rsd(:,index), MLE.p{index}{:})));
                    NLOGL = -sum(LL);
                    
%                     NLOGL = std(normcdf(rsdNormal(:,index)-cM_vec))
            end
    end
end

end