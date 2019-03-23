#!/bin/sh

for trial in "1" "2" "3" "4" "5"; do
    echo "------trial $trial"
    for mod in "o" "h" "s" "d"; do
        ./modulatingRC high $mod $trial
    done
done
