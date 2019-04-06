#!/bin/sh

./make.sh

for simType in max low high; do
    for trial in {1..5}; do
        echo "------trial $trial"
        for mod in d s h o; do
            ./modulatingRC $simType $mod $trial
        done
    done
done
