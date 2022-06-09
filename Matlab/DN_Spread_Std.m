function out = DN_Spread_Std(y)

    % no combination of single functions
    coder.inline('never');
    
    out = DN_Spread_Std(y, 5);