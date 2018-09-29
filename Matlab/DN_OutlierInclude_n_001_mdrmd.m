function out = DN_OutlierInclude_n_001_mdrmd(y)

    % no combination of single functions
    coder.inline('never');

    out = DN_OutlierInclude_001_mdrmd(y, 'n');
