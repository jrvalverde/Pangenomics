if [ ! -e metadata.full.tsv ] ; then
    echo $(head -n 1 ../GyrA/zyg_metadata.tsv)" reference outgroup" \
    | tr ' ' '	' \
    > metadata.full.tsv
    tail -n +2 ../GyrA/zyg_metadata.tsv \
        | while read id f s83 d87 s83l d87n s83ld87n res ; 
        do
            id=${f%.gyrA*}
            echo "$id	$id.fasta	$s83	$d87	$s83l	$d87n	$s83ld87n	$res	N	N"
        done \
        >> metadata.full.tsv

    echo "Efergusonii,Efergusonii.fasta,N,N,N,N,N,N,N,Y" | tr ',' '	' >> metadata.full.tsv
    echo "Reference,EcoliK12.fasta,N,N,N,N,N,N,Y,N" | tr ',' '	' >> metadata.full.tsv
fi

for i in 0*/* ; do
    echo $i
    head -n 1 metadata.full.tsv > $i/zyg_metadata.tsv
    for fa in $i/*.fasta ; do
        #echo $fa
        name=`basename $fa`
        #grep -m 1 $name metadata.full.tsv 
        grep -m 1 $name metadata.full.tsv >> $i/zyg_metadata.tsv
    done
    cat $i/zyg_metadata.tsv | tr '	' ',' > $i/zyg_metadata.csv
done
