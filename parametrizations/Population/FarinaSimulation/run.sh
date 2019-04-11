#!/bin/sh

./make.sh

for trial in {1..1}; do
    echo "------trial $trial"
    for mod in s o; do
        echo "-------mod $mod"
        ./Farina $mod $trial
    done
done
