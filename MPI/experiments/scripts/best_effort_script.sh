#!/bin/sh

for ((i = 0; i < 100; i += 10))
do
    start= ${i}
    final= ${i} + 9
    ./run_experiments ${start} ${final} > ../data/output_${i}.csv
end

