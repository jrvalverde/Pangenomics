hostname
date

# link in all target genomes
ln -s ../fasta-named/* .

# get list of genome files
ls -1 > zyg_ls_-1.log

# this we do afterwards to avoid including it in the list (this is the reference)
ln -s ../zygREF-Ec-ST131-no-plasmids.fasta .

# use the list of genomes to build a guide TAB file
grep fasta zyg_ls_-1.log \
| while read line ; do
    echo "${line%_*}    $line"
done > zyg_input.tab

export PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/sbin
export PATH=$PATH:~/work/jblazquez/src/snippy/snippy-master/bin:~/bin:~/.local/bin

# use snippy-multi to generate a script to run snippy
snippy-multi zyg_input.tab \
        --ref zygREF-Ec-ST131-no-plasmids.fasta
        --cpus 50 > zyg_snippy.sh 

# run snippy and snippy-core
bash zyg_snippy.sh |& tee zyg_snippy.log

for i in core* ; do cp $i zyg_$i ; done

# remove non-sequence characters changing them by N
snippy-clean_full_aln zyg_core.full.aln > zyg_core.full.clean.aln

export PATH=$PATH:~/contrib/miniconda3/bin

# run Gubbins with 50 threads
echo "Running Gubbins"
run_gubbins.py --verbose \
        --c 50 \
        --prefix zyg_gubbins \
        zyg_core.full.clean.aln \
        |& tee zyg_gubbins.log


echo "Computing SNPs and cleaning alignment"
snp-sites -c zyg_gubbins.filtered_polymorphic_sites.fasta > zyg_core.full.clean.gubbins.aln

echo "Computing tree with FastTreeMP"
FastTree=fasttreeMP     # use OpenMP version
$FastTree -gtr -nt zyg_core.full.clean.gubbins.aln > zyg_core.full.clean.gubbins.tree

date

exit

