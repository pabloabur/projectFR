#!/bin/sh

./make.sh

for simType in max; do # should be simType in max low high for complete simulation
    for trial in {11..20}; do
        echo "------trial $trial"
        for mod in s o; do
            ./modulatingRC $simType $mod $trial
        done
    done
done
