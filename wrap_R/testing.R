library(catch22);
library(foreign);

data = read.table("../testData/test.txt");
dataVec = c(t(data));

catch22_out = catch22_all(dataVec);

print(catch22_out)