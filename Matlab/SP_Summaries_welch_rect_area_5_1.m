function out = SP_Summaries_welch_rect_area_5_1(y)

    % no combination of single functions
    coder.inline('never');

    outStruct = SP_Summaries_welch_rect(y);
    out = outStruct.area_5_1;