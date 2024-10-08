#!/bin/bash
#
# run in ../sra-contigs
#
SEQUENTIAL="YES"


for gene in gyrA ; do
    output=blast.$gene
    mkdir $output

    reference=../gyrA/gyrA.fasta

    if ! makeblastdb -in $reference -dbtype nucl 
    then
        exit    
    fi
    

    for seq in fasta/*.fasta ; do
	echo $seq
	name=`basename $seq .fasta`
        if [ -s $output/$name.b7 ] ; then continue ; fi
	blastn -query $seq -db $reference -out $output/$name.b7 -outfmt 7
    done

done

exit
#------------------------------
