catch22_all <- function(data, catch24 = FALSE){

	names <- c('DN_OutlierInclude_n_001_mdrmd',
			   'DN_OutlierInclude_p_001_mdrmd',
			   'DN_HistogramMode_5',
				'DN_HistogramMode_10',
				'SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1',
				'SB_TransitionMatrix_3ac_sumdiagcov',
				'FC_LocalSimple_mean1_tauresrat',
				'SB_MotifThree_quantile_hh',
				'CO_HistogramAMI_even_2_5',
				'CO_Embed2_Dist_tau_d_expfit_meandiff',
				'SB_BinaryStats_diff_longstretch0',
				'MD_hrv_classic_pnn40',
				'SB_BinaryStats_mean_longstretch1',
				'FC_LocalSimple_mean3_stderr',
				'SP_Summaries_welch_rect_area_5_1',
				'SP_Summaries_welch_rect_centroid',
				'CO_f1ecac',
				'CO_FirstMin_ac',
				'IN_AutoMutualInfoStats_40_gaussian_fmmi',
				'PD_PeriodicityWang_th0_01',
				'SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1',
				'CO_trev_1_num')

	if(catch24){
      names <- append(names, c("DN_Mean", "DN_Spread_Std"))
	}

	values = c()

	for (feature in names){
		fh = get(feature)
		values = append(values, fh(data))
	}

	outData = data.frame(names = names, values = values)
	return(outData)
}