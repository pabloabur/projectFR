#!/bin/sh

./make.sh

for trial in {1..5}; do
    echo "------trial $trial"
    for mod in s o; do
        echo "-------mod $mod"
        ./FrequencyAnalysis $mod $trial
    done
done
