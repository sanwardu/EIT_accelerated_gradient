function [J] = calc_jacobian(sig,dsig,p,e,t,body)
%%%%% Computing Jacobian using Complete Electrode Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Sanwar Ahmad, suahmad@colostate.edu, Colorado State University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% sig - electrical conductivity structure
% p,t - Delaunay triangulation for mesh 
% e: node indices on each electrodes
% dsig: purturbation
% OUTPUTS
% J - sensitivity matrix or jacobian matrix from complete electrode model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% perturbing electrical conductivity model
NumSrc = body.NumSrc;
Ic = body.current;
npat = size(Ic,2);
nodes = size(p,2);
% Apply perturbation dsig to the nodes corresponding to each element (triangle)
 Nt = size(t,2); % Nt: # of elements
 ds = sparse(nodes,Nt);
 for i = 1:Nt
     ds(t(1,i),i) = dsig;
     ds(t(2,i),i) = dsig;
     ds(t(3,i),i) = dsig;
 end
ds(e(1,:),1:Nt) = 0;
[~,alphaT,Kn] = fwd_solver_eit2D(p,e,t,sig,body); % running CEM
% Computing Jacobian
for j = 1:Nt % change nodes to nt (# of elements)
    dU = []; 
    
    [fj,~,~]=assema(p,t,pdeintrp(p,t,ds(:,j)),0,0);
    rhs = [-fj*alphaT; sparse(NumSrc,npat)]; 
    thetaj = Kn\rhs; % calculating electric potential coefficients
    betaTj = thetaj(nodes+1:end,:); % potential on the boundary at electrodes
    for i = 1:npat
        dU = [dU; betaTj(:,i)]; % computing change in measured potential at electrodes
    end
    J(:,j) = dU; % storing each set of dUs as we perturb each node
end
end
