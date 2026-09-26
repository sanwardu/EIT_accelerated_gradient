function [sigTrue,sigTrue1] = ex3(p,p1,inArea,outArea,body)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Computes True conductivity distribution inside the domain
% one circular, one rectangular and one triangular target
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
nodes1=length(p);
nodes=length(p1);
rc = body.rc;
el = body.el;
center1 = [(rc-2.5*el)*cos(pi/2) (rc-2.5*el)*sin(pi/2)];
radius = el/.95;

for i=1:nodes
    x = p1(1,i); y = p1(2,i);
    if norm(center1'-p1(:,i))<=radius
        sigTrue(i)=inArea;
    
    elseif (-0.1 <= x) && (x <= -0.08) && (-0.07 <= y) &&(y <= -0.01 )
        sigTrue(i)= inArea;
    elseif (0.1 + y <= x) && (x <= 0.04 - y) && (-0.06 <= y) %&&(y <= -0.01 )
        sigTrue(i)= inArea;        
    else
        sigTrue(i)=outArea; 
    end
end

for i=1:nodes1
    x = p1(1,i); y = p1(2,i);
    if norm(center1'-p(:,i))<=radius
        sigTrue1(i)=inArea;

    elseif (-0.1 <= x) && (x <= -0.08) && (-0.07 <= y) &&(y <= -0.01 )
        sigTrue1(i)= inArea;   
    elseif (0.1 + y <= x) && (x <= 0.04 - y) && (-0.06 <= y)
        sigTrue1(i)= inArea;
        
    else
        sigTrue1(i)=outArea; 
    end
end
