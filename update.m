function J = update(Jc,pk,st,body)
% % % % % % % % % % % % % % 
% Accelerated update of the jacobian
% Jc: current jacobiam
% pk: search direction
% st: step length
% J: updated jacobian
% % % % % % % % % % % % % % 

NumSrc = body.NumSrc;
Ic = body.current;
npat = length(Ic(1,:));
nt = length(Jc(1,:));
J = sparse(NumSrc*npat,nt);
Jk = Jc' + st.*Jc'*(Jc*pk);
J = Jk';
end
