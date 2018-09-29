function out = CO_trev_1_num(y)

% no combination of single functions
coder.inline('never');

outStruct = CO_trev(y, 1);

out = outStruct.num;