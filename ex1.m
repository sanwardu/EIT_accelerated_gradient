function [sigTrue,sigTrue1] = ex1(p,p1,inArea,outArea,body)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Computes True conductivity distribution inside the domain
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
nodes1=length(p);
nodes=length(p1);
rc = body.rc;
el = body.el;
center1 = [-0.0508/.6 0];
radius = el/.95;
for i=1:nodes
    if norm(center1'-p1(:,i))<=radius
        sigTrue(i)=inArea;
        
    else
        sigTrue(i)=outArea; 
    end
end

for i=1:nodes1
    if norm(center1'-p(:,i))<=radius
        sigTrue1(i)=inArea;
    else
        sigTrue1(i)=outArea; 
    end
end
