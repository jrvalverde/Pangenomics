tail -n +3 zyg-ncbi-e_coli-assembly_summary.txt \
|    cut -d'	' -f20 \
|   while read url ; do
        if [[ "$url" == "#"* ]] ; then continue ; fi
        echo $url
        name=${url##*/}
        echo $name
        #continue
        if [ ! -d "$name" ] ; then
            mkdir "$name"
        fi
        cd "$name"
        if [ ! -s "${name}_genomic.fna.gz" ] ; then
            wget "$url/${name}_genomic.fna.gz" #-O "${name}_genomic.fna.gz"
        fi
        if [ ! -s "${name}_genomic.gff.gz" ] ; then
            wget "$url/${name}_genomic.gff.gz"
        fi
        if [ ! -s "${name}_genomic.gtf.gz" ] ; then
            wget "$url/${name}_genomic.gtf.gz"
        fi
        if [ ! -s "${name}_cds_from_genomic.fna.gz" ] ; then
            wget "$url/${name}_cds_from_genomic.fna.gz"
        fi
        if [ ! -s "${name}_protein.faa.gz" ] ; then
            wget "$url/${name}_protein.faa.gz"
        fi
        cd -
        sleep 1
    done \
|& tee zygGET.log
