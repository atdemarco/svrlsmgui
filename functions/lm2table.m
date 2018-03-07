function html = lm2table(mdl,caption)

%% Pretty css ...
% table { 
%   color: #333;
%   font-family: Helvetica, Arial, sans-serif;
%   width: 640px;
%   /* Table reset stuff */
%   border-collapse: collapse; border-spacing: 0; 
% }
% 		
% td, th {  border: 0 none; height: 30px; }
% 			
% th {
%   /* Gradient Background */
% 	background: linear-gradient(#333 0%,#444 100%);
% 	color: #FFF; font-weight: bold;
% 	height: 40px;
% }
% 		
% td { background: #FAFAFA; text-align: center; }
% 
% /* Zebra Stripe Rows */
% 		
% tr:nth-child(even) td { background: #EEE; } 
% tr:nth-child(odd) td { background: #FDFDFD; }
% 
% /* First-child blank cells! */
% tr td:first-child, tr th:first-child { 
% 	background: none; 
% 	font-style: italic;
% 	font-weight: bold;
% 	font-size: 14px;
% 	text-align: right;
% 	padding-right: 10px;
% 	width: 80px; 
% }
% 
% /* Add border-radius to specific cells! */
% tr:first-child th:nth-child(2) { 
%   border-radius: 5px 0 0 0; 
% } 
% 
% tr:first-child th:last-child { 
%   border-radius: 0 5px 0 0; 
% }
% 
% <style>
% table, th, td {
%     border: 1px solid black;
% }
% </style>

% f=figure;a=axes(f);
% mdl.plot('parent',a)

% can we pass a model title in the model object?
% the goal here is to imitate something like Mdl.display in html format...

html = createheaders({'Term','Estimate','SE','T Stat','P Value'});
html = add_coefficient_data(html,mdl);
html = add_gen_lm_info(html,mdl);

html = wrap_w_table_html(html,caption);% finish this off...

function rout = nice_r(rho) 
    rout = strrep(sprintf('r = %.2f', rho),'0.','.'); % enforce no leading unnecessary 0 
    
function rout = nice_r_nolet(rho) 
    rout = strrep(sprintf('%.2f', rho),'0.','.'); % enforce no leading unnecessary 0 

function pout= nice_p(pval) 
    pout= myif(pval < .001,'P < .001',strrep(sprintf('P = %.2f', pval),' 0.',' .')); % enforce no leading unnecessary 0 

function pout= nice_p_nolet(pval) % ie no p = or p < ...
    pout= strrep(myif(pval < .001,'<.001',sprintf('%.3f', pval)),'0.','.'); % enforce no leading unnecessary 0 

function html = add_gen_lm_info(html,mdl)
    % Add in the num obs, err df, rmse, rsquared, adj rsquared, f--stat ...

    [fstat_vs_constant_pval,fstat_vs_constant_fval] = coefTest(mdl);
    [dw_P,dw_DW] = mdl.dwtest;
    
    % new rows to add in this cell array...
    rows = {add_col(5,'',sprintf('Number of observation: %d, Error degrees freedom: %d',mdl.NumObservations,mdl.DFE)), ...
            add_col(5,'',sprintf('Root Mean Square Error: %0.2f',mdl.RMSE)), ...
            add_col(5,'',sprintf('R&#178;: %s, Adjusted R&#178;: %s',nice_r_nolet(mdl.Rsquared.Ordinary),nice_r_nolet(mdl.Rsquared.Adjusted))), ...
            add_col(5,'',sprintf('F-statistic vs. constant model: %0.2f, %s',fstat_vs_constant_fval,nice_p(fstat_vs_constant_pval))), ...
            add_col(5,'',sprintf('Durbin-Watson test for corr. resids, DW = %0.2f, %s',dw_DW,nice_p(dw_P)))};
    
    for r = 1 : numel(rows) % go through and append each...
        newrow = rows{r};
        %newrow = ['<font size="9">' newrow '</font>'];
        html = [html wrap_w_row_html(newrow)]; % add the row wrapper.
    end
    
    %     Add in the durbin watson test - info from http://www.statisticshowto.com/durbin-watson-test-coefficient/
    %     The Hypotheses for the Durbin Watson test are:
    %     H0 = no first order autocorrelation.
    %     H1 = first order correlation exists.
    %     (For a first order correlation, the lag is one time unit).
    %
    %     Assumptions are:
    %     That the errors are normally distributed with a mean of 0.
    %     The errors are stationary.
    % 
    % The Durbin Watson test reports a test statistic, with a value from 0 to 4, where:
    % 
    % 2 is no autocorrelation.
    % 0 to <2 is positive autocorrelation (common in time series data).
    % >2 to 4 is negative autocorrelation (less common in time series data).


function html = add_coefficient_data(html,mdl)
    ncoeffs = numel(mdl.CoefficientNames);
    for c = 1 : ncoeffs
        newrow = [];
        newrow = add_col(1,newrow,mdl.CoefficientNames{c}); % Coefficient name
        newrow = add_col(1,newrow,mdl.Coefficients.Estimate(c)); % Estimate
        newrow = add_col(1,newrow,mdl.Coefficients.SE(c)); % SE
        newrow = add_col(1,newrow,mdl.Coefficients.tStat(c)); % tStat
        newrow = add_col(1,newrow,nice_p_nolet(mdl.Coefficients.pValue(c))); % pValue
        html = [html wrap_w_row_html(newrow)]; % add the row wrapper.
    end

function rowdat = add_col(span,rowdat,newval)
    if span > 1 % e.g., %colspan="2"
        spanstring = [' colspan="' num2str(span) '"'];
    else
        spanstring ='';
    end
    
    switch class(newval)
        case 'double'
            htmladd = sprintf('<td%s>%0.3f</td>',spanstring,newval); % insert span here as necessary
        case 'char'
            htmladd = sprintf('<td%s>%s</td>',spanstring,newval); % insert span here as necessary
        otherwise
            error('unknown data class')
    end
    rowdat = [rowdat htmladd]; % concat and return

%    {'(Intercept)'}    {'x1'}    {'x2'}    {'x3'}    {'x4'}    {'x5'}
%mdl.anova
%             SumSq     DF    MeanSq      F         pValue  
%              ______    __    ______    ______    __________
% 
%     x1       116.77     1    116.77    108.57    2.3674e-17
%     x2       490.59     1    490.59    456.14    7.7393e-38
%     x3       728.86     1    728.86    677.67    9.3215e-45
%     x4         1189     1      1189    1105.5    9.0239e-54
%     x5       1832.2     1    1832.2    1703.5    4.9317e-62
%     Error     101.1    94    1.0755                        

% addTerms               dwtest                 plotEffects            random                 
% anova                  feval                  plotInteraction        removeTerms            
% coefCI                 plot                   plotPartialDependence  step                   
% coefTest               plotAdded              plotResiduals          
% compact                plotAdjustedResponse   plotSlice              
% disp                   plotDiagnostics        predict                

function html = createheaders(headers)
    html = [];
    for h = 1 : numel(headers)
        curheader = headers{h}; % {'Coefficient','Estimate','SE','tStat','pValue'};
        html = [html '<th>' curheader '</th>'];
    end
    html = wrap_w_row_html(html);

function wrapped = wrap_w_row_html(html)
    wrapped = [newline '<tr>' newline html newline '</tr>'];

function wrapped = wrap_w_table_html(html,caption)
    cap = [newline '<caption>' caption '</caption>' newline];
    wrapped = ['<center><table style="width:80%">' cap html '</table></center>'];