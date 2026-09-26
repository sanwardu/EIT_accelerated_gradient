function B = localstiff(p,t)
% Computing the integral int grad phi_i . grad phi_j dx over a triangle
% for linear basis functions
% % % % % % % % % % % % % % % % 
% Author: Sanwar Ahmad, suahmad@colostate.edu
% May 2020
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% OUTPUT:
% B: Structure that stores the local stiffness matrices
% % % % % % % % % % % % % % % % % % % % % % % % % % % % 

nt = size(t,2);
[ar,~,~,~]= pdetrg(p,t); % ar: area of each triangular element
for m = 1:nt
    ind = t(1:3,m); % indices of the m-th element
    y1 = p(2,ind(1)); x1 = p(1,ind(1));
    y2 = p(2,ind(2)); x2 = p(1,ind(2));
    y3 = p(2,ind(3)); x3 = p(1,ind(3));
    
    G = [y2-y3 y3-y1 y1-y2; x3-x2 x1-x3 x2-x1]; % [grad phi_1 grad phi_2 grad phi_3]   
    for i = 1:3
        gri = G(:,i); % grad phi_i
        for j = 1:3
            grj = G(:,j); % grad phi_j
            B(m).stiff(i,j) = gri'*grj/(4*ar(m)); % local stiffness matrix
        end
    end
end
