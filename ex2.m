function [sigTrue,sigTrue1] = geom5(p,p1,inArea,outArea,body)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Computes True conductivity distribution inside the domain
% One conductive and one resistive target
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
nodes1=length(p);
nodes=length(p1);
rc = body.rc;
el = body.el;
center1 = [(rc-2.5*el)*cos(pi/3) (rc-2.5*el)*sin(pi/3)];
center2 = [(rc-2.5*el)*cos(pi+pi/4) (rc-2.5*el)*sin(pi+pi/4)];
radius = el/.95;
for i=1:nodes
    if norm(center1'-p1(:,i))<=radius
        sigTrue(i)=inArea;
        
    elseif norm(center2'-p1(:,i))<=radius
        sigTrue(i)= 0.02;        
    else
        sigTrue(i)=outArea; 
    end
end

for i=1:nodes1
    if norm(center1'-p(:,i))<=radius
        sigTrue1(i)=inArea;        
    elseif norm(center2'-p(:,i))<=radius
        sigTrue1(i)= 0.02;         
    else
        sigTrue1(i)=outArea; 
    end
end
