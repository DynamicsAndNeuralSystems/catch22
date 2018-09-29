% testFilter

fileID = fopen('/Users/carl/PycharmProjects/catch22/C/timeSeries/tsid0133.txt','r');
tsData = fscanf(fileID,'%f');
fclose(fileID);

[b, a] = myButter(4, 0.8/24);

tsDatafilt = filter(b,a,tsData)