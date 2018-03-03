function p_vec = betas2pvals(betavec,tail) 
    % the input variable betavec should NOT be sorted already
    switch tail
        case 'pos' % highest (beta) value have significant p values
            [~,p] = sort(betavec,'descend');
        case 'neg' % smallest (beta) values has significant p values
            [~,p] = sort(betavec,'ascend');
        case 'two'
            error('not supported atm')
    end   
   
    r = 1:length(betavec);
    r(p) = r;
    p_vec = r/length(betavec);
    
    
%% Visual demo of the high scores good:
%     [~,p] = sort(betavec,'descend');
%     r = 1:length(betavec);
%     r(p) = r;
%     p_vec = r/length(betavec);
%     figure;
%     scatter(betavec,testpvals);
%     xlabel('beta val');
%     ylabel('p value')
