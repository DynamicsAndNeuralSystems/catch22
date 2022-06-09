#!/bin/bash

help()
{
   echo ""
   echo "Usage: $0 -i indir -o outdir -a append_string -s"
   echo -e "\t-i Path to a directory containing input time series files (.txt with one time series value per line). Default: './'"
   echo -e "\t-o Path to an existing directory in which to save output feature values. Default: './'"
   echo -e "\t-a A string appended to the input file names to create the output file names. Default: '_output'"
   echo -e "\t-s Evaluate catch22 (0) or catch24 (1). Default: 0"
   exit 1
}

while getopts "i:o:a:s:" opt
do
   case "$opt" in
      i) indir="$OPTARG" ;;
      o) outdir="$OPTARG" ;;
      a) append="$OPTARG" ;;
      s) catch24="$OPTARG" ;;
      ?) help ;; # Print help in case parameter is non-existent
   esac
done

srcdir=$(dirname "$0}")

# Print help in case parameters are empty
if [ -z "$indir" ]
then
   indir="./"
fi

if [ -z "$outdir" ]
then
   outdir=$indir
fi

if [ -z "$append" ]
then
   append="_output"
fi

if [ -z "$catch24" ]
then
   catch24=0
fi

indir=$(basename "$indir")
outdir=$(basename "$outdir")

# Loop through each file in indir and save the feature outputs
for entry in "${indir}"/*
do
    filename=$(basename "$entry")
    extension="${filename##*.}"
    filename="${filename%.*}"
    fullfile="${outdir}/${filename}${append}.${extension}"
    if [ "${filename: -${#append}}" != "${append}" ]
    then
        yes $catch24 | "${srcdir}/run_features" $entry $fullfile

        if [ -s $fullfile ]
        then
            echo $fullfile
        else
            rm $fullfile
        fi

    fi
done
