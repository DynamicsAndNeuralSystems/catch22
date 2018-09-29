function out = SB_TransitionMatrix_3ac_sumdiagcov(y)

% no combination of single functions
coder.inline('never');

outStruct = SB_TransitionMatrix(y, 3, 'ac');

out = outStruct.sumdiagcov;