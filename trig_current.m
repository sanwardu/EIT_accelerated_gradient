function Ic = trig_current(body,c)
% c: Amplitude of the current
% Ic: input current 

nl = body.NumSrc; % number of source
npat = nl-1; % number of pattern
Ic = zeros(nl,npat);
theta = 2*pi/nl.*(1:nl); % electrode positions
for i = 1:nl/2
    Ic(:,i) = c*cos(i*theta);
end

for i = (nl/2+1):npat
    Ic(:,i) = c*sin((i-nl/2)*theta);
end
