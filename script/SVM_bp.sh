#!/bin/bash
input_fa=$1
mountme=$( cd "$( dirname "${input_fa}" )" && pwd )
file=`basename $input_fa`
cmd="docker run --rm -it -v $mountme:/svmbp_stage mlovci/svm_bpfinder /usr/bin/SVM_BP/svm_bpfinder.py $file $2 $3"
eval $cmd | tail -n+2
