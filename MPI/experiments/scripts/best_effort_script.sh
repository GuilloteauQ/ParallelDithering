#!/bin/bash

for i in {1..20..10}
do
    lower_bound=${i}
    upper_bound=$((i+10))
    oarsub -t besteffort -l node=1 "./run_experiments.R ${lower_bound} ${upper_bound} > ../data/output_${i}.csv"
done
