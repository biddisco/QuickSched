#!/bin/bash

# Executable and parameters
prog=test_qr
params="-m 32 -n 32 -r 10"

# Main loop
for n in {1..64}
do

    if [ ! -e ${prog}_${n}.dump ]
    then
        ./${prog} ${params} -t ${n} > ${prog}_${n}.dump
    fi

done
