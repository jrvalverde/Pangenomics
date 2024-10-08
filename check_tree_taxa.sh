
tree_file=zyg_gubbins.final_tree.tre
metadata=zyg_metadata.tsv

cat $tree_file \
| sed -e 's/,/\n/g' -e 's/(//g' \
| sed -e 's/:.*//g' \
| while read taxon ; do
    #echo $taxon
    if ! grep -q $taxon $metadata ; then
        echo "$taxon missing in $metadata"
    fi
done
