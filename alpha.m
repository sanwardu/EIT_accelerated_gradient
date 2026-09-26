function [s,Jn,V2] = alpha(Vmeas,V,Jc,pk,sig1,sig_prior,lambda,p,e,t,obj_err,in_st,body) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Last updated: May. 2025, Sanwar Ahmad, sahmad@vsu.edu, VSU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the step length at kth iteration
%  Input:
%  Vmeas: synthetic data; 
% V: approx data for sig_k
%  Jc: Jacobian evaluated at sig_k
% pk: search direction
%  sig1: sig_k
%  sig_prior: prior conductivity
%  lambda: Regularization parameter at kth iteration
%  in_st: Initial step length
%  obj_err: residula error for sig_k
% Output: 
% s: step length; pk: direction at kth step
% Jn: Jacobian at sig_(k+1)
% V2: approx data for sig_(k+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
it_st = 16; % maximum number of inner iteration
s = abs(in_st);
Nt = size(t,2);
it = 1;
fv1 = -Jc'*(V-Vmeas); % gradient J at sig1
stmin = 10^-4;

W2 = speye(Nt,Nt);
fv = -(fv1 + lambda*W2*(sig1-sig_prior)); 
q1 = 10^-4; q2 = 0.9;
dsig = 1e-4;
while it < it_st     
    sig2 = sig1 + s.*pk; % corrected    
    [V2,~,~] = fwd_solver_eit2D(p,e,t,sig2,body);
    new_obj_err = norm(V2-Vmeas,2)^2;
    % update jacobian
    Jn = cal_jacobian(sig2,dsig,p,e,t,body); % regular jacobian update at sig2
    % Jn = update(Jc,pk,s,body); % accelerated jacobian updatr at sig2 
      fv2 = - Jn'*(V2-Vmeas); % gradient J at sig2
    
    % Checking strong Wolf's conditions
    wolf2nd_2 = fv2'*pk; %corrected
    wolf2nd_1 = q2.*(fv1')*pk; %corrected
    if new_obj_err <= (obj_err + q1*s*(-fv')*pk) && abs(wolf2nd_2) <= abs(wolf2nd_1) %corrected
        break
    end
	it = it + 1; % updating iteration count         
    if s/2 > stmin %updating step length
        s = s/2;
    else
        break
    end
end
