function [Vpred,alphaT,Kn] = fwd_solver_eit2D(p,e,t,sig,body)

% Computes the solution of the 2D forward EIT problem
% Sanwar Ahmad, suahmad@colostate.edu
% INPUT
% p,t: mesh information, nodes and elements
% e: node on the boudary
% sig: conductivity distribution at each nodes
% Ic: Input current
% body: structure, it contains information about source, CP, current vector, contact impedance, electrode length etc.
% OUTPUT
% alphaT: potential distribution at p
% Vpred: simulated/computed voltages
% Kn: stiffness matrix
% % % % % % % % % % % % % % % % % % % % % % 

np = length(p(1,:));
nt = size(t,2);
NumSrc = body.NumSrc; % number of electrodes
Ic = body.current;
npat = size(Ic,2);
elecInd = body.elec;

%% Global Stiffness Matrix Assembly 
% refer: Ph.D. thesis, "ITERATIVE IMAGE RECONSTRUCTION FOR ELECTRICAL IMPEDANCE TOMOGRAPHY USING ADAPTIVE TECHNIQUES"
% by Taoran Li (2014) for more information

% Computing K
if length(sig) == np
    sig = pdeintrp(p,t,sig);
end

K1 = zeros(np,np);
Kl = localstiff(p,t);

for m = 1:nt
    ind = t(1:3,m);
    for i = 1:3
        for j = 1:3
            K1(ind(i),ind(j)) = K1(ind(i),ind(j)) + sig(m)*Kl(m).stiff(i,j);
        end
    end
end

[B,C,D] = fem_matrices(p,elecInd,body); % other two global stiffness components

Ac = K1 + B; % global stiffness component
 
% finally, assemble global stiffness matrix with shunt and contact impedance models 
Kn = [Ac C; C' D];
KK = [Kn; zeros(1,np) ones(1,NumSrc)]; % Block matrix
%% Generating data and taking inversion for linear basis coefficients

% transforming data
Fn=[sparse(np,npat); Ic]; % rhs 
Fn = [Fn; zeros(1,npat)];
% solving for linear basis coefficients that will be used to
% reconstruct the predicted interior and measured electric potentials 
theta = KK\Fn;
% for internal potential 
alphaT = theta(1:np,:); 
% for potential on the boundary at the electrodes
betaT = theta(np+1:end,:); 
 
%% Computing electric potential and then displaying results

U = betaT;
adj = sum(U)/NumSrc;
adjust = zeros(NumSrc,npat);
adjust(1:npat,:) = meshgrid(adj);
adjust(NumSrc,:) = adj;
Vpred = U - adjust;

Vpred = reshape(Vpred,NumSrc*npat,1);
end
