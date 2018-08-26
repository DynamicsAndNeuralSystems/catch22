function out = catch22(y)

out = zeros(22,1);

out(1) = DN_HistogramMode_5(y);
out(2) = DN_HistogramMode_10(y);
out(3) = CO_f1ecac(y);
out(4) = CO_FirstMin_ac(y);
out(5) = CO_HistogramAMI_even_2_5(y);
out(6) = CO_trev_1_num(y);
out(7) = SB_TransitionMatrix_41_ondiag(y);
out(8) = SB_TransitionMatrix_2ac_sumdiagcov(y);
out(9) = SB_TransitionMatrix_3ac_maxeigcov(y);
out(10) = PD_PeriodicityWang_th0_01(y);
out(11) = CO_Embed2_Dist_tau_d_expfit_meandiff(y);
out(12) = PH_Walker_prop_01_sw_propcross(y);
out(13) = FC_LocalSimple_mean1_tauresrat(y);
out(14) = DN_OutlierInclude_p_001_mdrmd(y);
out(15) = DN_OutlierInclude_n_001_mdrmd(y);
out(16) = EX_MovingThreshold_01_01_meankickf(y);
out(17) = SB_BinaryStats_diff_longstretch0(y);
out(18) = SB_MotifThree_quantile_hh(y);
out(19) = SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1(y);
out(20) = SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1(y);
out(21) = SP_Summaries_welch_rect_centroid(y);
out(22) = CO_trev_2_denom(y); % replacement for  MF_arfit_1_8_sbc_meanA