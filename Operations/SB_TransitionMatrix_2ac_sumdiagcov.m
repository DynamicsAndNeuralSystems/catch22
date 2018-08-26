function out = SB_TransitionMatrix_2ac_sumdiagcov(y)

outStruct = SB_TransitionMatrix(y, 2, 'ac');

out = outStruct.sumdiagcov;