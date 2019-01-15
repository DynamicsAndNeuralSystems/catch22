import catch22

for dataFile in ['../testData/test.txt', '../testData/test2.txt']:

    print '\n',  dataFile

    data = [line.rstrip().split(' ') for line in open(dataFile)]
    flat_data = [float(item) for sublist in data for item in sublist]

    catchOut = catch22.catch22_all(flat_data)

    featureNames = catchOut['names']
    featureValues = catchOut['values']

    for featureName, featureValue in zip(featureNames, featureValues):
        print '%s : %1.6f' % (featureName, featureValue)