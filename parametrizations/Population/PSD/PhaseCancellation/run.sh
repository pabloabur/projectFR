#!/bin/sh

./make.sh

for trial in {1..10}; do
    echo "------trial $trial"
    for mod in rc dc; do
        echo "-------mod $mod"
        ./FrequencyAnalysis $mod $trial
    done
done
