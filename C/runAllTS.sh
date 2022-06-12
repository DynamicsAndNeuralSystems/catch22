#!/bin/bash

help()
{
   echo ""
   echo "Usage: $0 -i indir -o outdir -a append_string -s"
   echo -e "\t-h Show this help message"
   echo -e "\t-i Path to a directory containing input time series files (.txt with one time series value per line). Default: './timeSeries'"
   echo -e "\t-o Path to an existing directory in which to save output feature values. Default: './featureOutput'"
   echo -e "\t-a A string (minus ext.) appended to the input file names to create the output file names. Default: 'output'"
   echo -e "\t-s Evaluate catch22 (0) or catch24 (1). Default: 0"
   exit 1
}

while getopts "i:o:a:s:h" opt
do
   case "$opt" in
      i) indir="$OPTARG" ;;
      o) outdir="$OPTARG" ;;
      a) append="$OPTARG" ;;
      s) catch24="$OPTARG" ;;
      h) help ;;
   esac
done

srcdir=$(dirname "$0}")

if [ -z "$indir" ]
then
   indir="./timeSeries"
fi

if [ -z "$outdir" ]
then
   outdir="./featureOutput"
fi

if [ -z "$append" ]
then
   append="output"
fi

if [ -z "$catch24" ]
then
   catch24=0
fi

indir="$(dirname $indir)/$(basename $indir)"
outdir="$(dirname $outdir)/$(basename $outdir)"
mkdir -p $outdir

# Loop through each file in indir and save the feature outputs
for entry in "${indir}"/*.txt
do
    filename=$(basename "$entry")
    extension="${filename##*.}"
    filename="${filename%.*}"
    fullfile="${outdir}/${filename}${append}.${extension}"
    if [ "${filename: -${#append}}" != "${append}" ]
    then
        yes $catch24 | "${srcdir}/run_features" $entry $fullfile > /dev/null

        if [ -s $fullfile ] # Remove file if catch22 errors
        then
            echo "Output written to ${fullfile}"
        else
            rm $fullfile
        fi

    fi
done
