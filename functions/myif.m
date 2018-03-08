function val = myif(teststatement,return_true,return_false)
if nargin < 3
    error('Insufficient arguments.')
end
% convenience function
if teststatement 
    val = return_true;
else
    val = return_false;
end
