function [sig_all,l2err_sig,st,obj,iter,lambda] = block_kaczmarz(Vmeas,sig_guess,sigTrue,J,lambda,numbit,p,e,t,Tol,body)
%%%%% Modified block Kaczmarz method for EIT Inverse Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% May. 2025, Sanwar Ahmad, sahmad@vsu.edu, VSU
%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% Vmeas - experimental or simulated data on voltages
% sig_guess - initial model for the electrical conductivity structure
% sigTrue - true E.C. model  
% lambda - model weighting parameter
% numbit - number of iterations
% p,e,t - Delaunay triangulation for mesh 
% Tol - tolerance
% body - stores the circular body properties
% OUTPUTS
% sig_all - All predicted electrical conductivity model
% l2err_sig - Relative Reconstruction error
% st - set of step lengths obtained by babacktracking
% obj - Residual error values
% iter - number of iteration taken
% lambda - model weighting parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kaczmarz-Method

Nt = size(t,2); % Nt: number of elements
sig_all = [sig_guess' zeros(Nt,numbit)]; % saving all runs of model 
obj = zeros(numbit,1); % obj = objective values at each iteration

sig_prior = sig_guess'; % model prior, assume equal weight for all parameters (homogeneous model)

% initializing some variables
iter = 1;
l2err_sig(iter) = norm(sig_all(:,iter)-sigTrue',2)/norm(sigTrue,2); % Reconstruction error

[Vpred,~,~] = fwd_solver_eit2D(p,e,t,sig_all(:,iter),body); % simulated data

obj(iter) = norm(Vpred-Vmeas,2)^2; % Residual error

nmeas = size(Vpred,1);
W1 = speye(nmeas,nmeas); % Weight matrix

lambda1 = lambda; c = 4;
st = [];
stmin = 10^-4;
% %%  Computing search direction and step length
    r = (Vmeas-Vpred);
    
    fv1 = J'*r; % gradient of cost function
    
    [step,pk] = optstep(r,J,lambda); % compute optimum step length
if step == 0
    return
end

while iter < numbit

    [step, J, Vpred] = alpha(Vmeas,Vpred,J,pk,sig_all(:,iter),sig_prior,lambda,p,e,t,obj(iter),step,body);
    st = [st; step];
    if (iter+1) > numbit
        fprintf('Maximum iteration exceeded');
        break;
    else
        iter = iter + 1; % updating iteration count         
    end
    
    sig_all(:,iter) = sig_all(:,iter-1) + step.*pk; % updating e.c. model
    
   % %%  Computing search direction
    lambda = lambda1*c/(c+iter);

    r = (Vmeas-Vpred);
    l2err_sig(iter) = norm(sig_all(:,iter)-sigTrue',2)/norm(sigTrue,2); % Relative reconstruction error
    obj(iter) = norm(r,2)^2; % Residual error

    pk = (J*J' + lambda.*W1)\r;
    pk = J'*pk;
    Jp = J*pk;
% Checking for stopping criterion     
    if obj(iter) < Tol
        fprintf('Discrepancy principle');
        break
    end
    if abs(obj(iter)/obj(iter-1)) > 1
            fprintf('Res error increasing')
            break
    end 
end
end
