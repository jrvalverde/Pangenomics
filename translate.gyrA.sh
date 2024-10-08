
mkdir -p ../gyrA
mkdir -p ../GyrA

i=0	# we'll start with no match
prevgenome=""

#grep "NC_002695.2" *.b7 .
grep "NC_002695.2" zygMATCHES \
| tr ':' ' ' \
| while read name contig gene gcoord identity score \
             mismatches gaps qstart qend tstart tend \
             evalue bscore ; do 
             
      genome=`basename $name .b7`
      fasta=../fasta/$genome.fasta
      echo ""
      echo -n "$fasta $i $contig $qstart $qend $tstart $tend "

      if [ "$genome" != "$prevgenome" ] ; then
          i=0		# reset counter
      fi
      prevgenome="$genome"	# for next iteration

      i=$((i + 1))	# we have a new match      
      # if already done, skip it
      if [ -s ../GyrA/$genome.gyrA.$i.faa ] ; then continue ; fi
     
      if [ $tstart -lt $tend ] ; then
          samtools faidx $fasta "$contig:$qstart-$qend" \
          > ../gyrA/$genome.gyrA.$i.fna

          transeq -sequence ../gyrA/$genome.gyrA.$i.fna \
                  -frame 1 \
                  -outseq ../GyrA/$genome.gyrA.$i.faa

      elif [ $tend -lt $tstart ] ; then
          samtools faidx $fasta "$contig:$qstart-$qend" \
          > ../gyrA/$genome.gyrA.$i.rev.fna
          
          revseq -sequence ../gyrA/$genome.gyrA.$i.rev.fna \
                 -outseq   ../gyrA/$genome.gyrA.$i.fna

          transeq -sequence ../gyrA/$genome.gyrA.$i.fna \
                  -frame 1 \
                  -outseq ../GyrA/$genome.gyrA.$i.faa
      else
          echo -n "match ignored!"
      fi
  done
