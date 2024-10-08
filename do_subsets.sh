#!/bin/bash

# generate names based on gyrA status
#R --vanilla -q -s <zygCMD.R

# convert to genome names
for i in 0050 0100 0250 0500 1000 ; do
    echo $i
    mkdir -p $i.genomes
    for j in $i/* ; do 
        cat $j | sed -e 's/.gyrA.*.faa//g' > $i.genomes/`basename $j`
    done
done
