function J = cem_jac(p,t,u,K,body)
% Input:
% p: nodes
% t: triangular elements
% sig: conductivity vector at each element
% u: FEM solution for sig
% K: inverse of the FEM matrix of size (N+L X N+L)
% body: cylindrical geometry
% Output:
% J: Jacobian
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Author: Sanwar Ahmad, suahmad@colostate.edu
% May 2020
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
L = body.NumSrc; % number of electrodes
npat = body.npat;
B = localstiff(p,t); % local stiffness matrices
np = length(p(1,:));
nt = length(t(1,:));
J = zeros(L*npat,nt); 

Ct = [zeros(L,np) eye(L,L)]; % (LX N+L) % not using any basis
M = Ct*K; 
G = M(:,1:np);

% Find the stiffness matrix
for m = 1: nt
    dB = B(m).stiff; % dB/ds for mth element
    dK = zeros(np,np);
    ind = t(1:3,m);
    for i = 1:3
        for j = 1:3
            dK(ind(i),ind(j)) = dB(i,j); % dK + B
        end
    end

    dU = - G*dK*u;
    J(:,m) = reshape(dU,L*npat,1);
end
end
