function out = decimate_(y, q)

    nfilt = 4;
    
    nd = length(y);

    [b,a] = myButter(nfilt, 0.8/q);

    % be sure to filter in both directions to make sure the filtered data has zero phase
    % make a data vector properly pre- and ap- pended to filter forwards and back
    % so end effects can be obliterated.
    out_filt = filtfilt_(b,a,y);
    nbeg = 1;
    out = out_filt((nbeg:q:nd)');
    