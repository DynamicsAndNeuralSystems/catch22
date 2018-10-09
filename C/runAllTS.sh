mkdir -p "featureOutputs"

for entry in ./timeSeries/*
do
	echo "$entry"
	
	filename=$(basename "$entry")
	extension="${filename##*.}"
	filename="${filename%.*}"

	./run_features $entry #"./featureOutputs/"$filename"output.txt"

done
