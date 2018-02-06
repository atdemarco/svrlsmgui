function yesno = isint(val)
    % is the input value 'val' an integer or not...
    % returns boolean true false.
    yesno = ~(val-ceil(val));