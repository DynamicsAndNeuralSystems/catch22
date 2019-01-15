library(catch22);
library(foreign);

data = read.table("../testData/test.txt");
dataVec = c(t(data));

features = read.table("../featureList.txt");
featureVec = c(t(features));

for (feature in featureVec){
    cat(paste(feature, ': '))
    fh = get(feature);
    cat(paste(fh(dataVec), '\n'))
}
