function zmap = p2z(pmap)
% converts p map to z map

% note that this is for 1-tailed analyses only!
% one tail
z_1tail = @(p) -sqrt(2) * erfcinv(p*2); % this anonymous function converts a 1 tailed p value to a z value.

zmap = z_1tail(pmap);

% two-tail
%z_2tail = @(p) abs(norminv(p/2)) 
