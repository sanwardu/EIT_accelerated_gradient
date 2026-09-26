function [A,G,D] = fem_matrices(p,e,body)
% % % % % % % % % % % % % % % % % % % % % % 
% Input:
% p: mesh points
% e: node indices on each electrodes
% Output:
% A: sum_l (1/zl) \int_el phi_i phi_j ds
% G: -1/zl \int_el phi_i ds 
% D: diag(1/zl (\int_el ds) = |el|/zl
% Reference: https://www.mathworks.com/help/releases/R2015a/pde/ug/elliptic-pdes.html
% Sanwar Ahmad, suahmad@colostate.edu
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
zl = body.zl; % contact impedance
el = body.el; % length of electrodes
nelec = body.NumSrc; % number of electrodes
np = length(p(1,:));
A = zeros(np,np); 
G = zeros(np,nelec);
for k = 1:nelec   
    ind = e(k).ind; % nodes on the kth electrode
    nl = size(ind,2); % number of nodes on each electrode
    for i = 1:(nl-1) % (nl-1): numer of elements on kth electrode        
        P1 = [p(1,ind(i)) p(2,ind(i))];
        P2 = [p(1,ind(i+1)) p(2,ind(i+1))];
        B = [norm(P1-P2,2)/3 norm(P1-P2,2)/6; norm(P1-P2,2)/6 norm(P1-P2,2)/3];      
        A(ind(i),ind(i)) = A(ind(i),ind(i))+ B(1,1);
        A(ind(i),ind(i+1)) = A(ind(i),ind(i+1))+ B(1,2);
        A(ind(i+1),ind(i)) = A(ind(i+1),ind(i))+ B(2,1);
        A(ind(i+1),ind(i+1)) = A(ind(i+1),ind(i+1))+ B(2,2);        
        
        G(ind(i),k) = G(ind(i),k) + norm(P1-P2,2)/(2);
        G(ind(i+1),k) = G(ind(i+1),k) + norm(P1-P2,2)/(2);  
    end
end
A = A./zl;    
G = -G./zl;

D = el/zl.*eye(nelec,nelec);
end
