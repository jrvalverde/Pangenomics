#!/bin/bash

num='0500'		# 0050 0250 500
bootsize=100
# we will loop over files in ../subsets/xxx.genomes

# define our preferences
ref=EcoliK12		# reference sequence
outgroup=Efergusonii	# outgroup
ALN=SKA			# SNIPPY SKA
TREEBUILDER=fasttree	# fasttree raxml raxmlng

export ref outgroup ALN TREEBUILDER

# first we will do only one of each
#for i in 0050 0100 0500 1000 ; do
#for i in 0050 0500 ; do
for i in "$num" ; do
    origdir=../subsets/$i.genomes
    
    for j in {1..$bootsize} ; do
        suffix=`printf "%04d" $j`
        workdir="$i/$suffix"
        
        if [ -d "$workdir" ] ; then continue ; fi
        if [ ! -e "$origdir/sensible.$suffix" ] ; then break ; fi
        
        echo "Creating $workdir job"
        mkdir -p $workdir
        mkdir -p $workdir/log
        
        sennames="$origdir/sensible.$suffix"
        resnames="$origdir/s83l.$suffix"
            
        # populate $workdir and generate auxiliary files
        echo "ID,file,sensible,S83L,REF,OUTGRP" > $workdir/zyg_metadata.csv
	truncate -s 0 $workdir/zyg_input.tab

        # add sensible genomes
        cat $sennames \
        | while read genome ; do
            fasta=`realpath ../fasta/$genome.fasta`
            # add genome to $workdir
            ln -sf $fasta $workdir/
            # add to metadata file
            echo "$genome,$genome.fasta,Y,N,N,N" >> $workdir/zyg_metadata.csv
            # add to input file
            echo "$genome	$genome.fasta" >> $workdir/zyg_input.tab
        done
        
        # repeat for resistant genomes
        cat $resnames \
        | while read genome ; do
            fasta=`realpath ../fasta/$genome.fasta`
            # add genome to $workdir
            ln -sf $fasta $workdir/
            # add to metadata file
            echo "$genome,$genome.fasta,N,Y,N,N" >> $workdir/zyg_metadata.short.csv
            # add to input file
            echo "$genome	$genome.fasta" >> $workdir/zyg_input.tab
        done

	# add the outgroup
        og='Efergusonii'
        ln -s `realpath ../outgroup/$og.fasta $workdir`/$og.fasta
        echo "$og,$og.fasta,Y,N,N,Y" >> $workdir/zyg_metadata.csv
        echo "$og	$og.fasta" >> $workdir/zyg_input.tab
        
        # now add the reference
        ref="EcoliK12"
        ln -s `realpath ../ref/$ref.fasta $workdir`/$ref.fasta
        echo "Reference,$ref.fasta,Y,N,Y,N" >> $workdir/zyg_metadata.short.csv
        echo "$ref	$ref.fasta" >> $workdir/zyg_input.tab

	# generate tab metadata
        cat $workdir/zyg_metadata.short.csv \
        | tr ',' '	' \
        > $workdir/zyg_metadata.short.tsv
        
        # generate full metadata
        head -n 1 metadata.full.tsv > $i/zyg_metadata.tsv
        for fa in $workdir/*.fasta ; do
            #echo $fa
            name=`basename $fa`
            #grep -m 1 $name metadata.full.tsv 
            grep -m 1 $name metadata.full.tsv >> $workdir/zyg_metadata.tsv
        done
        cat $workdir/zyg_metadata.tsv | tr '	' ',' > $workdir/zyg_metadata.csv
        
        
        # copy scripts
        #cp do_gubbins.sh $workdir/zyg_gubbins.sh
        # make a unique copy of the script in each folder (to ensure
        # we have a unique copy of the script used to calculate this).
        cat do_gubbins.sh \
        | sed -e "s|%TAG%|$i/$j|g" \
        > $workdir/zygGUBBINS.sh	# seems fastest
        #cp do_gubbins_snippy.sh $workdir/zyg_gubbins.sh
        
        cp ./generate_ska_alignment.py $workdir/
    done
done



# We want some fine control of the next steps for this is likely to be 
# a very computationally heavy bootstrap

start=${start:-3}		# for running in the ILLUMINA cluster (1..$nnodes)
nnodes=${nnodes:-4}		# ditto (4 nodes)
njobs=${njobs:-50}
echo "Processing from $start to $njobs by $nnodes"
for n in $( seq $start $nnodes $njobs ) ; do
    i="$num"/`printf "%04d" $n`
#for i in [01]*/* ; do
    echo $i
    #continue
    cd $i
    # defaults will be taken from the exported variables
    bash zygGUBBINS.sh |& sed -e "s|^|${i}: |g" | tee log/zygGUBBINS.log
    cd -
done

nnodes=10
njobs=50
for st in $( seq 1 $nnodes ) ; do
  start=$st
  echo "Processing jobs from $start to $njobs using $nnodes nodes"
  (
    for n in $( seq $start $nnodes $njobs ) ; do
        i="$num"/`printf "%04d" $n`
    #for i in 00*/* ; do
    #for i in 02*/* ; do
    #for i in 05*/* ; do
        cd $i
        pwd

        if [ -e DONE ] ; then cd - ; continue ; fi
        if [ ! -e zyg_gubbins.final_tree.tre.bck ] ; then
            cp zyg_gubbins.final_tree.tre zyg_gubbins.final_tree.tre.bck
        fi
        if [ ! -e zyg_gubbins.final_tree.tre.sav ] ; then
            cp zyg_gubbins.final_tree.tre zyg_gubbins.final_tree.tre.sav
	    cat zyg_gubbins.final_tree.tre.sav \
            | sed -e 's/\bEcoliK12\b/Reference/' \
            > zyg_gubbins.final_tree.tre
        fi
        cp ../../zygRCandy.R .
        cp  ../../../ref/Escherichia_coli_str._K-12_substr._MG1655,_complete_genome.gff3 \
            EcoliK12.gff3
        mkdir -p zyg_R
        Rscript zygRCandy.R |& tee zygRCandy.log    
        touch DONE

        cd -
    done
  ) > log/R_"$num".$st.log 2>&1 &
done
