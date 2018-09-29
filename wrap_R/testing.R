library(catch22);
library(foreign);

data = read.table("../testData/test.txt");
dataVec = c(t(data));

features = read.table("../featureList.txt");
featureVec = c(t(features));

for (feature in featureVec){
    print(feature)
    fh = get(feature);
    print(fh(dataVec))
}
