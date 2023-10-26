function MCMC = adptMCMC(data, like, p0, max_iter, init_stepsize)
%Compute max_iter iterations of an adaptive MCMC chain with likelihood 
%function like. Starting at p0.

if nargin < 5
    init_stepsize = 0.05;
end

%% Initialisation
% rng(12345);

iter = 1;

MCMC.burnIn = 1000;
MCMC.adpItr = 1000;

%% Pre-allocation
new_p = zeros(length(p0),0); % parameter estimates in matrix

%% Starting value
new_p(:,1) = p0;

%% Step-size
% h = [0.015*ones(1,n_cp+1);0.025*ones(1,n_cp+1);0.03*ones(1,n_cp+1);0.015*ones(1,n_cp+1)];
stepsize.h_in = init_stepsize * ones(size(p0));
stepsize.h_ad = init_stepsize * ones(size(p0));
% h(1,:) = 0.01;
% h = [0.05*ones(1,n_cp),0;0.05*ones(1,n_cp),0;0.05*ones(1,n_cp+1);0.05*ones(1,n_cp+1)];

%% MCMC
MCMC.accept_vec = zeros(max_iter,1);
MCMC.like = zeros(max_iter,1);
while iter <= max_iter
    %% Print progress
    if mod(iter, max_iter) == 0; fprintf('\n\n');
    elseif mod(iter, 1e5) == 0; fprintf('\n');
    elseif mod(iter, 1000) == 0; fprintf('+');
    end
    
    %% Get current parameter values
    curr_prm = new_p(:,iter);    
       
    %% Calculate/Update covariance matrix
    if iter == MCMC.burnIn + MCMC.adpItr + 1 % we don't need it before anyway
        stepsize.mu_upd = mean(new_p(:,MCMC.burnIn:end),2);
        stepsize.adp_h = cov(new_p(:,MCMC.burnIn:end)');
    elseif iter > MCMC.burnIn + MCMC.adpItr + 1
        if mod(iter,1000) == 0 % just making sure, the update does not suffer from numerical errors
            stepsize.mu_upd = mean(new_p(:,MCMC.burnIn:end),2);
            stepsize.adp_h = cov(new_p(:,MCMC.burnIn:end)');
        else
            stepsize.mu_upd = stepsize.mu_upd + 1/(iter-MCMC.burnIn+1) * (new_p(:,iter) - stepsize.mu_upd);
            stepsize.adp_h = (iter-MCMC.burnIn)/(iter-MCMC.burnIn+1)*stepsize.adp_h + (iter-MCMC.burnIn)/(iter-MCMC.burnIn+1)^2*(new_p(:,iter) - stepsize.mu_upd)*(new_p(:,iter)-stepsize.mu_upd)';
        end
    end
    
    %% Propose new parameters
    if iter <=  MCMC.burnIn +  MCMC.adpItr % fixed proposal for the initial iterations
        proposed = curr_prm + randn(size(p0)) .* stepsize.h_in;
    else % adaptive MCMC proposal
        if length(p0) == 2
            proposed = curr_prm+(1-0.05)*mvnrnd(zeros(size(p0)),stepsize.adp_h,1)'*2.38/sqrt(length(p0))+0.05*stepsize.h_ad/sqrt(length(p0)).*randn(length(p0),1);
        else
            proposed = curr_prm+(1-0.05)*mvnrnd(zeros(size(p0)),stepsize.adp_h,1)'*2.38/sqrt(length(p0))+0.05*stepsize.h_ad/sqrt(length(p0)).*randn(length(p0),1);
        end
    end
    if iter == 50000
        debug = 1;
    end
    
    %% Accept the new parameters
    if iter > 1
        NLOGL_old = MCMC.like(iter-1);
    else
        NLOGL_old = like(data, curr_prm);
    end
    NLOGL_new = like(data, proposed);
    accept = exp(NLOGL_old - NLOGL_new);
    
    if rand < accept
        curr_prm = proposed;
        MCMC.accept_vec(iter) = 1;
        MCMC.like(iter) = NLOGL_new;
    else 
        MCMC.accept_vec(iter) = 0;
        MCMC.like(iter) = NLOGL_old;
    end
    if iter == 1
        MCMC.accept_rate(iter) = MCMC.accept_vec(iter);
    elseif iter <= 100
        MCMC.accept_rate(iter) = ((iter - 1) * MCMC.accept_rate(iter-1) + MCMC.accept_vec(iter))/iter;
    else
        MCMC.accept_rate(iter) = MCMC.accept_rate(iter-1) + (MCMC.accept_vec(iter) - MCMC.accept_vec(iter-100))/100;
    end
    
    %% Save the new (resp. old) parameters if accepted (rejected)
    new_p(:,iter+1) = curr_prm;
    
    %% Next iteration
    iter = iter + 1;
end

MCMC.data = data;
MCMC.max_iter = max_iter;
MCMC.p = new_p;

[~,I] = min(MCMC.like);
MCMC.MLE = MCMC.p(:,I+1);
