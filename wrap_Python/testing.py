import catch22

# get feature names from package content
# features = dir(catch22)
# features = [item for item in features if not '__' in item]

# alternatively, get the feature names from the same source as the other wrappers
with open('../featureList.txt') as f:
    features = f.readlines()
features = [x.strip() for x in features] 

for dataFile in ['../testData/test.txt', '../testData/test2.txt']:

    print '\n',  dataFile

    data = [line.rstrip().split(' ') for line in open(dataFile)]
    flat_data = [float(item) for sublist in data for item in sublist]

    for testFun in features:
        featureFun = getattr(catch22, testFun)
        print '%s : %1.6f' % (testFun, featureFun(flat_data))