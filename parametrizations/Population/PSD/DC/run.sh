#!/bin/sh

./make.sh

for trial in {1..10}; do
    echo "------trial $trial"
    for mod in d s h o; do
        echo "-------mod $mod"
        for mvc in 05 70; do
            echo "--------mvc $mvc"
            ./FrequencyAnalysis $mod $mvc $trial
        done
    done
done
