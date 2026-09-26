function [sig_all,l2err_sig,st,obj,iter] = IRGN(Vmeas,sig_guess,sigTrue,J,lambda,numbit,p,e,t,Tol,body)
%%%%% Iteratively Regularized Gauss-Newton Function for EIT Inverse Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Mar. 2020, Sanwar Ahmad, suahmad@colostate.edu, CSU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS 
% Vmeas - experimental or simulated data on voltages
% sig_guess - initial model for the electrical conductivity structure
% sigTrue - true E.C. model 
% J - Jacobian
% lambda - model weighting parameter
% numbit - total number of iterations
% p,e,t - Delaunay triangulation for mesh 
% Tol - Tolerance
% body - structure 
% OUTPUTS
% sig_all - All predicted electrical conductivity model
% l2err_sig - Reconstruction error
% st - set of step lengths obtained by babacktracking
% obj - Residual error values
% iter - final iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt = size(t,2); % Nt: number of elements
sig_all = [sig_guess' zeros(Nt,numbit)]; % saving all runs of model 
obj = zeros(numbit,1); % obj = objective values at each iteration
sig_prior = sig_guess'; % model prior, assume equal weight for all parameters (homogeneous model)
W2 = speye(Nt,Nt); % model regularization matrix 
% initializing some variables
iter = 1;
l2err_sig(iter) = norm(sig_all(:,iter)-sigTrue',2)/norm(sigTrue,2); % Reconstruction error

[Vpred,~,~] = fwd_solver_eit2D(p,e,t,sig_all(:,iter),body);
obj(iter) = norm(Vpred-Vmeas,2)^2; % Residual error
W1 = speye(Nt,Nt); % Weight matrix
lambda1 = lambda; c = 4;
st = [];
stmin = 10^-4;
%%  Computing search direction and step length
    r = (Vpred-Vmeas);
    fv1 =  J'*r; % gradient J at sig1
    fv = -(fv1 + lambda*W1*(sig_all(:,iter)-sig_prior)); 
    pk = (J'*J + lambda.*W1)\fv;
    Jp = J*pk;
    stepmax = - (r'*Jp)/(Jp'*Jp)

while iter < numbit 
    [s, J,Vpred] = alpha(Vmeas,Vpred,J,pk,sig_all(:,iter),sig_prior,lambda,p,e,t,obj(iter),abs(stepmax),body);
    
    if (iter+1) > numbit
        fprintf('Maximum iteration exceeded');
        break;
    else
        iter = iter + 1; % updating iteration count         
    end
    
    sig_all(:,iter) = sig_all(:,iter-1) + s.*pk; % updating e.c. model
    % %%  Computing search direction and step length
    lambda = lambda1*c/(c+iter);    
    r = (Vpred-Vmeas);
    fv1 =  J'*r; % gradient J at sig1
    fv = -(fv1 + lambda*W1*(sig_all(:,iter)-sig_prior)); 
    pk = (J'*J + lambda.*W1)\fv;
    l2err_sig(iter) = norm(sig_all(:,iter)-sigTrue',2)/norm(sigTrue,2); % Reconstruction error    
    obj(iter) = norm(r,2)^2; % Residual error
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
