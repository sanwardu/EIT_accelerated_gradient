function elec = elec_ind(body,p,e)

% Find the boundary nodes on each electrode
% Sanwar Ahmad, suahmad@colostate.edu, sahmad@vsu.edu
% INPUT
% p: nodes
% e: node on the boudary
% OUTPUT
% elec: structure, contains arrays of boundary nodes on each electrode


% Sorting the boundary element indices from [0,2pi)
boundIndPos = [e(1,:)' p(:,e(1,:))'];

% finding x,y positions in first and second quadrant
[ind12,~] = find(boundIndPos(:,3)>=0);
IndPos12 = boundIndPos(ind12,:);
[~,ind12] = sort(IndPos12(:,2),'descend');
IndPos12 = IndPos12(ind12,:);

% in third and fourth quadrant
[ind34,~] = find(boundIndPos(:,3)<0);
IndPos34 = boundIndPos(ind34,:);
[~,ind34] = sort(IndPos34(:,2),'ascend');
IndPos34 = IndPos34(ind34,:);

% combining all x,y positions 
boundIndPos = [IndPos34; IndPos12]; % to sort in [-pi, pi) reverse the order
% boundIndPos = [IndPos12; IndPos34]; % to sort in [0, 2pi) reverse the order
% Find the indices covering electrodes
el = body.el; % length of the electrodes
nelec = body.NumSrc; % number of electrodes
nbd = size(e,2);
ncover = nbd/nelec; % ncover: least number of nodes on each electrode, nelec must divide nbd
leftind = boundIndPos(1:ncover:end,1,1);
for i = 1:nelec
    k = 1;
    d = 0;
    ind1 = leftind(i); 
    elec(i).ind(k) = ind1;
    pos = find(ind1 == e(1,:)); % index of the left end
    ind2 = e(2,pos); % adjacent node to the left end
    P1 = [p(1,ind1) p(2,ind1)];
    P2 = [p(1,ind2) p(2,ind2)];
    d = d + norm(P1-P2,2);
    while d <= el
        k = k + 1;
        elec(i).ind(k) = ind2;
        ind1 = ind2; % update left end
        P1 = P2;
        pos = find(ind1 == e(1,:)); % index of the left end
        ind2 = e(2,pos); % adjacent node to the left end
        P2 = [p(1,ind2) p(2,ind2)];
        d = d + norm(P1-P2,2);  
    end
end 
