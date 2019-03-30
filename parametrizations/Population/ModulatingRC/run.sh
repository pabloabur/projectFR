#!/bin/sh

for trial in {1..5}; do
    echo "------trial $trial"
    for mod in d s h o; do
        ./modulatingRC max $mod $trial
    done
done
