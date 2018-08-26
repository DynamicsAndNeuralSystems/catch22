function out = SB_TransitionMatrix_3ac_maxeigcov(y)

outStruct = SB_TransitionMatrix(y, 3, 'ac');

out = outStruct.maxeigcov;