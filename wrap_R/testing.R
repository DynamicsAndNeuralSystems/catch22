library(catch22)
library(foreign)

data = read.table("../testData/test.txt")
dataVec = c(t(data))

catch22_out = catch22_all(dataVec)
catch24_out = catch22_all(dataVec, TRUE)

print(catch22_out)
print(catch24_out)