function zmap = p2z(pmap)
% converts p map to z map

if any(pmap(:)>1) 
    error('Input data must be p values, which cannot exceed 1.')
elseif any(pmap(:) <= 0)
    min(pmap(:))
    error('Input data must be p values, which cannot be less than or equal to zero.')    
end

% note that this is for 1-tailed analyses only!
% one tail
z_1tail = @(p) -sqrt(2) * erfcinv(p*2); % this anonymous function converts a 1 tailed p value to a z value.

zmap = z_1tail(pmap);

% two-tail
%z_2tail = @(p) abs(norminv(p/2)) 
