#!/bin/bash

UNIQUE_TAG="%TAG%"		# this will be substitued on creation

# set some default values
ref=${ref:-EcoliK12}
og=${outgroup:-Efergusonii}
ALN=${ALN:-SKA}
TREEBUILDER=${TREEBUILDER:-fasttree}
NPROC=$(( `nproc` / 2 ))

#snp_sites=$HOME/bin/snp-sites

export PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/sbin
#export PATH=$PATH:~/work/jblazquez/src/snippy/snippy-master/bin
export PATH=:$PATH:~/.local/bin			# pip installed
export PATH=$PATH:~/contrib/miniconda3/bin	# conda installed
export LANG=C
source ~/contrib/anaconda3/etc/profile.d/conda.sh
conda activate

if [ -e zRUNNING ] ; then
    echo "A computing job is already running"
    echo "    "`pwd`
    echo "    ...skipping"
    echo "    (kill it and remove file zRUNNING if you want to restart it)"
    exit
elif [ -e zDONE ] ; then
    echo "This jobs seems to be already completed"
    echo "    "`pwd`
    echo "    ...skipping"
    echo "    (remove file zDONE is you want to recalculate it)"
    exit
#elif [ -e zFAIL* ] ; then we will try again to finish the calculation
#else we have to carry out the calculation
fi

echo ">>> host: " $(hostname) |& tee zRUNNING
echo ">>> start:" $(date) |& tee -a zRUNNING

echo ">>> Computing core alignmente using $ALN"
date
if [ "$ALN" == "SNIPPY" ] ; then
    # Compute core phylogeny using SNIPPY (snippy-multi) or SKA
    #export PATH=$PATH:~/work/jblazquez/src/snippy/snippy-master/bin:~/bin
    alignment=zyg_core.full.snippy.aln

    if [ ! -s "$alignment" ] ; then
        # use snippy-multi to generate a script to run snippy
        snippy-multi zyg_input.tab \
                --ref $ref.fasta
                --cpus $NPROC > zyg_snippy.sh 

        # run snippy and snippy-core
        bash zyg_snippy.sh |& tee zyg_snippy.log

        for i in core* ; do ln $i zyg_$i ; done

        # remove non-sequence characters changing them by N
        snippy-clean_full_aln zyg_core.full.aln > $alignment

        echo ">>> snippy:" $(date) |& tee -a zRUNNING
    fi
elif [ "$ALN" == "SKA" ] ; then
    # compute core alignment using SKA
    alignment=zyg_core.full.ska.aln
    
    if [ ! -s "$alignment" ] ; then
    
        # first try to save time (these alignments will hopefully not be
        # calculated)
        for i in *.fasta ; do
            # we use the same reference in all cases, so the
            # alignment should be the same as well
            g=${i%.fasta}
            if [ ! -s $g.skf -a -s ../../../ska/$g.skf ] ; then
                # force overwrite to avoid errors (and to use latest alignment)
                ln -f ../../../ska/$g.skf .
            else
                rm -f $g.skf
            fi
            if [ ! -s $g.map.aln -a  -s ../../../ska/$g.map.aln ] ; then
                ln -f ../../../ska/$g.map.aln .
            else
                rm -f $g.map.aln
            fi
        done

	generate_ska_alignment=./generate_ska_alignment.py
        python3 $generate_ska_alignment \
            --reference $ref.fasta \
            --fasta zyg_input.tab \
            --threads $NPROC \
            --out $alignment
        ok=$?

       for i in *.skf ; do
           g=${i%.skf}
           if [ ! -s ../../../ska/$g.skf ] ; then
               ln -f $g.skf ../../../ska/$g.skf
           fi
           if [ ! -s ../../../ska/$g.map.aln ] ; then
               ln -f $g.map.aln ../../../ska/$g.map.aln
           fi
       done
        
        echo ">>> ska:" $(date) |& tee -a zRUNNING

	if [ $ok -ne 0 ] ; then
            mv zRUNNING zFAIL.ska
            exit
        fi
    fi
else
    echo "Unknown alignment method"
    mv zRUNNING zFAIL.ALN
    exit
fi

# run Gubbins with $nproc threads

if [ ! -s zyg_gubbins.filtered_polymorphic_sites.fasta ] ; then
    echo "Running Gubbins"
    
    treebuilder="$TREEBUILDER"
    #treebuilder=fasttree
    #treebuilder=raxml
    #treebuilder=raxmlng
    
    # this will create (among others) a $prefix.log (zyg_gubbins.log)
    run_gubbins.py --verbose \
            --threads $NPROC \
            --prefix zyg_gubbins \
            ${comment# --no-cleanup} \
            --tree_builder $treebuilder ${comment# version 2.4.1 in the cluster } \
            --filter_percentage 35 \
            ${comment# --tree-builder $treebuilder }  ${comment# version 3.0.0 } \
            ${comment# --first-tree-builder $treebuilder } \
            ${comment# --filter-percentage 35 } \
            ${comment# --outgroup $og } \
            $alignment 
        ok=$?

        echo ">>> gubbins:" $(date) |& tee -a zRUNNING

	if [ $ok -ne 0 ] ; then
            mv zRUNNING zFAIL.gubbins
            exit
        fi
fi

echo ">>> end:" $(date) |& tee -a zRUNNING
mv zRUNNING zDONE
exit


# the next steps are only needed is we want to build an additional tree
# which we do not need tright now.

if [ ! -s zyg_core.full.gubbins.aln ] ; then
    echo "Computing SNPs and cleaning alignment"
    # -c means 'only output columns containing only ATGC' 
    # works i 2.3.3 and above
    # 2.3.3 is installed in NGS, 2.1 in XISTRAL
    # for 2.3+
    snp-sites -c zyg_gubbins.filtered_polymorphic_sites.fasta > zyg_core.full.gubbins.aln
    # for 2.1
    #snp-sites zyg_gubbins.filtered_polymorphic_sites.fasta > zyg_core.full.gubbins.aln
    ok=$?
    
    echo ">>> snp:" $(date) |& tee -a zRUNNING
    if [ $ok -ne 0 -o ! -s zyg_core.full.gubbins.aln ] ; then
        mv zRUNNING zFAIL.snp
        exit
    fi
fi 

if [ ! -s zyg_core.full.gubbins.tree ] ; then
    echo "Computing tree with FastTreeMP"
    #alias FastTree=fasttreeMP     # use OpenMP version
    alias FastTree="VeryFastTree -threads $NPROC -fastexp 2"
    FastTree -gtr -nt zyg_core.full.gubbins.aln > zyg_core.full.gubbins.tree
    ok=$?
    
    echo ">>> tree:" $(date) |& tee -a zRUNNING
    if [ $ok -ne 0 -o ! -s zyg_core.full.gubbins.tree ] ; then
        mv zRUNNING zFAIL.tree
        exit
    fi
fi

echo ">>> end:" $(date) |& tee -a zRUNNING
mv zRUNNING zDONE
exit

