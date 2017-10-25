function success = CreateDirectory(path)
% Just a small wrapper to conditionally create directory if it doesn't exist.
success = 1;
try, if ~exist(path,'dir'), mkdir(path); end %#ok<NOCOM>
catch, success = 0;
end