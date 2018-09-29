function out = SP_Summaries_welch_rect_centroid(y)

    % no combination of single functions
    coder.inline('never');

    outStruct = SP_Summaries_welch_rect(y);
    out = outStruct.centroid;
