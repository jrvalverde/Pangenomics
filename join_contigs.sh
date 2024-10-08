#!/bin/bash
#
#	Given a multi-fasta file containing many contigs as separate
# fasta entries, e.g.
# >contig1
# ATGC.......
# >contig2
# ATGC....
# >etc
#
# join all the contigs in the same order as a single sequence in a
# new file stored in a parallel directory ../joined_contigs with the
# same name and a single fasta entry named like the file (without the
# .fasta extension
#

file=$1
name=`basename $1 .fasta`

echo ">$name" > ../joined_contigs/$file
grep -v '>' $file >> ../joined_contigs/$file

