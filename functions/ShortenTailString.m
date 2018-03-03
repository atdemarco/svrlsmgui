function parameters = ShortenTailString(parameters)
    % for convenience....
    switch parameters.tails % get a shortened tail name for convenience
        case  {'One-tailed (negative)','High scores are good'}
            parameters.tailshort = 'neg';
        case {'One-tailed (positive)','High scores are bad'}
            parameters.tailshort = 'pos';
        otherwise
            parameters.tailshort = 'two';
    end
