#!/bin/sh

./make.sh

for simType in low high; do # should be simType in max low high for complete simulation
    for trial in {6..10}; do
        echo "------trial $trial"
        for mod in d s h o; do
            ./modulatingRC $simType $mod $trial
        done
    done
done
