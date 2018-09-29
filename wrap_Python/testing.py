import catch22

features = dir(catch22)
features = [item for item in features if not '__' in item];

for dataFile in ['../testData/test.txt', '../testData/test2.txt']:

    print '\n',  dataFile

    data = [line.rstrip().split(' ') for line in open(dataFile)]
    flat_data = [float(item) for sublist in data for item in sublist]

    for testFun in features:
        featureFun = getattr(catch22, testFun)
        print '%s : %1.6f' % (testFun, featureFun(flat_data))