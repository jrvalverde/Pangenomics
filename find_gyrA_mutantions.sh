#!/bin/bash


#echo $1
#echo $2
#echo $@
#exit

GENOMES_FAS=fasta	# directory with the genome fasta files

# Blast against gyrA (the nucleotide sequence of the gyrA gene)
cd $GENOMES_FAS

bash ../script/blast.gyrA.sh |& tee zygLOG	# do the blast against GyrA
mv blast.gyrA ..

cd ..

# analyze blast results
cd ../blast.gyrA

# find genomes with more than one match
#   grep -c will give the count of hits with "NC_002695.2" in each genome
#   grep -v ':1$' will remove the genomes with exactly one hit
grep -c "NC_002695.2" *.b7 | grep -v ":1$" > zyg+1  

# we inspected the output and there were a few but only one will give us
# trouble
# in all case the gene was split across several contigs
# in one case, the part of the gene corresponding to hotspot 83(protein)
# was missing ( NC_002695.2" )
# in two cases the gene was split in three parts ( SRR3098899_contigs.b7
# and SRR3098953_contigs.b7 ) which in both cases, after inspection, 
# seems to be due to fragmented reconstrucion of the genome 
cat zyg+1 | while IFS=':' read name count ; do
    echo "------------- $name $count ------------------"
    echo "query.acc.ver, subject.acc.ver, %identity, alignment-length, mismatches, gap-opens, q.start, q.end, s.start, s.end, evalue, bitscore"
    grep "NC_002695.2" $name
done > zyg+1.hits


# extract genome-name and start/end of match *in* GyrA
grep -H  "NC_002695.2" * |tr ':' '\t' | cut -f 1,11,12 > zygGyrStEnd 

# check which genomes have matches not starting/ending at position 1 of the 
# gene
grep -v -w '1' zygGyrStEnd > zygNo1

read -p "Check by hand zygNo1 and do any corrections before continuing" yn
# we checked by hand and saw that all cases corresponded to secondary
# fragments from the splitted genes

# save the id of the trouble genome so we can process it by hand later
echo "ERR871394_contigs.b7" > zygINCOMPLETE   

# when we grep more than one file, the output contains the name of the
# file a ':' and the matching line(s)

# create directories for gene and protein sequences
mkdir -p ../gyrA
mkdir -p ../GyrA

# extract the gene sequence to ../gyrA, and translate it saving the protein
# sequence to ../GyrA

grep "NC_002695.2" *.b7 \
| tr ':' ' ' \
| while read name contig gene gcoord identity score \
             mismatches gaps qstart qend gstart gend \
             evalue bscore ; do 
      echo $name $contig $qstart $qend $gstart $gend 
      
      genome=`basename $name .b7`
      fasta=../fasta/$genome.fasta
      echo $fasta $contig $qstart $qend $gstart $gend 
      
      if [ $gstart -eq 1 ] ; then
          samtools faidx $fasta "$contig:$qstart-$qend" \
          > ../gyrA/$genome.gyrA.fna

          transeq -sequence ../gyrA/$genome.gyrA.fna \
                  -frame 1 \
                  -outseq ../GyrA/$genome.gyrA.faa

      elif [ $gend -eq 1 ] ; then
          samtools faidx $fasta "$contig:$qstart-$qend" \
          > ../gyrA/$genome.gyrA.rev.fna
          
          revseq -sequence ../gyrA/$genome.gyrA.rev.fna \
                 -outseq   ../gyrA/$genome.gyrA.fna

          transeq -sequence ../gyrA/$genome.gyrA.fna \
                  -frame 1 \
                  -outseq ../GyrA/$genome.gyrA.faa
      else
          echo "match ignored"
      fi
  done

cd ..
# we are done processing the blast output
# the gyrA directory containing the gene sequences is not further used, 
# but we keep it so we can manually check that the genes are correctly
# extracted, and in case we later want to analyze the genomic gene sequences.


# move to the protein directory and check for mutant patterns
#	IMPORTANT NOTE: we have previously made sure manyally that these 
# patterns cover all the possible cases. We did by ensuring that no non-wild
# type genome had a different pattern (i.e. that there are no other mutations
# besides the ones we are looking for in the matched fragment in any genome).
#
cd GyrA

grep GKYHPHGD.AVY *.faa | grep -v GKYHPHGDSAVY > zyg_S83L_mutants.txt


grep AVY.TIVRMAQ *.faa | grep -v AVYDTIVRMAQ > zyg_D87N_mutants.txt

cat *.nam | sort | uniq | wc -l

comm -2 -3 zyg_S83L_mutants.nam zyg_D87N_mutants.nam > zyg_S83L_only_mutants.nam

comm -1 -3 zyg_S83L_mutants.nam zyg_D87N_mutants.nam > zyg_D87N_only_mutants.nam

comm -1 -2 zyg_S83L_mutants.nam zyg_D87N_mutants.nam > zyg_S83L_D87N_mutants.nam

cat *.nam | sort | uniq > zyg_S83L_D87N_any_mutants.txt


cd ..

# we are done
