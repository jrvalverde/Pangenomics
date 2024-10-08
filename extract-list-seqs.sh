find ../genomes -name '*_genomic.fna.gz' -print \
| while read fasta ; do
    if [[ "$fasta" == *"cds_from"* ]] ; then continue ; fi
    name=`basename "$fasta" _genomic.fna.gz`
    echo "Uncompressing $fasta to $name.fasta"
    zcat "$fasta" > ./"$name".fasta
done

find . -name '*.fasta' -exec grep '^>' {} \; > zygENTRIES

cat zygENTRIES | sort | uniq > zygENTRIES.s.u

cat zygENTRIES | sort | uniq -c > zygENTRIES.s.uc
