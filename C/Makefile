
OBJ:= ./bin

run_features: $(OBJ)/main.o $(OBJ)/butterworth.o $(OBJ)/CO_AutoCorr.o $(OBJ)/DN_HistogramMode_5.o $(OBJ)/DN_HistogramMode_10.o $(OBJ)/DN_OutlierInclude.o $(OBJ)/FC_LocalSimple.o $(OBJ)/fft.o $(OBJ)/helper_functions.o $(OBJ)/histcounts.o $(OBJ)/IN_AutoMutualInfoStats.o $(OBJ)/MD_hrv.o $(OBJ)/PD_PeriodicityWang.o $(OBJ)/SB_BinaryStats.o $(OBJ)/SB_CoarseGrain.o $(OBJ)/SB_MotifThree.o $(OBJ)/SB_TransitionMatrix.o $(OBJ)/SC_FluctAnal.o $(OBJ)/SP_Summaries.o $(OBJ)/splinefit.o $(OBJ)/stats.o
	gcc $(OBJ)/main.o $(OBJ)/butterworth.o $(OBJ)/CO_AutoCorr.o $(OBJ)/DN_HistogramMode_5.o $(OBJ)/DN_HistogramMode_10.o $(OBJ)/DN_OutlierInclude.o $(OBJ)/FC_LocalSimple.o $(OBJ)/fft.o $(OBJ)/helper_functions.o $(OBJ)/histcounts.o $(OBJ)/IN_AutoMutualInfoStats.o $(OBJ)/MD_hrv.o $(OBJ)/PD_PeriodicityWang.o $(OBJ)/SB_BinaryStats.o $(OBJ)/SB_CoarseGrain.o $(OBJ)/SB_MotifThree.o $(OBJ)/SB_TransitionMatrix.o $(OBJ)/SC_FluctAnal.o $(OBJ)/SP_Summaries.o $(OBJ)/splinefit.o $(OBJ)/stats.o -o $@

$(OBJ)/main.o: main.c main.h 
	gcc -c $< -o $@

$(OBJ)/butterworth.o: butterworth.c butterworth.h
	gcc -c $< -o $@

$(OBJ)/CO_AutoCorr.o: CO_AutoCorr.c CO_AutoCorr.h
	gcc -c $< -o $@

$(OBJ)/DN_HistogramMode_5.o: DN_HistogramMode_5.c DN_HistogramMode_5.h
	gcc -c $< -o $@

$(OBJ)/DN_HistogramMode_10.o: DN_HistogramMode_10.c DN_HistogramMode_10.h
	gcc -c $< -o $@

$(OBJ)/DN_OutlierInclude.o: DN_OutlierInclude.c DN_OutlierInclude.h
	gcc -c $< -o $@

$(OBJ)/FC_LocalSimple.o: FC_LocalSimple.c FC_LocalSimple.h
	gcc -c $< -o $@

$(OBJ)/fft.o: fft.c fft.h
	gcc -c $< -o $@

$(OBJ)/helper_functions.o: helper_functions.c helper_functions.h
	gcc -c $< -o $@

$(OBJ)/histcounts.o: histcounts.c histcounts.h
	gcc -c $< -o $@

$(OBJ)/IN_AutoMutualInfoStats.o: IN_AutoMutualInfoStats.c IN_AutoMutualInfoStats.h
	gcc -c $< -o $@ 

$(OBJ)/MD_hrv.o: MD_hrv.c MD_hrv.h 
	gcc -c $< -o $@

$(OBJ)/PD_PeriodicityWang.o: PD_PeriodicityWang.c PD_PeriodicityWang.h 
	gcc -c $< -o $@

$(OBJ)/SB_BinaryStats.o: SB_BinaryStats.c SB_BinaryStats.h 
	gcc -c $< -o $@

$(OBJ)/SB_CoarseGrain.o: SB_CoarseGrain.c SB_CoarseGrain.h 
	gcc -c $< -o $@

$(OBJ)/SB_MotifThree.o: SB_MotifThree.c SB_MotifThree.h 
	gcc -c $< -o $@

$(OBJ)/SB_TransitionMatrix.o: SB_TransitionMatrix.c SB_TransitionMatrix.h 
	gcc -c $< -o $@

$(OBJ)/SC_FluctAnal.o: SC_FluctAnal.c SC_FluctAnal.h 
	gcc -c $< -o $@

$(OBJ)/SP_Summaries.o: SP_Summaries.c SP_Summaries.h
	gcc -c $< -o $@

$(OBJ)/splinefit.o: splinefit.c splinefit.h 
	gcc -c $< -o $@

$(OBJ)/stats.o: stats.c stats.h 
	gcc -c $< -o $@

clean:
	rm $(OBJ)/* run_features