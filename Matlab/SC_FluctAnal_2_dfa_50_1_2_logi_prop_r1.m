function out = SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(y)

    % no combination of single functions
    coder.inline('never');

    out = SC_FluctAnal_q2_taustep50_k1_logi_prop_r1(y, 'dfa', 2);