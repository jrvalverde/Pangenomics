# Select sequences not wild-type
grep GKYHPHGD.AVY *.gyrA.*.faa | grep -v GKYHPHGDSAVY > zyg_S83_mutants.txt

grep AVY.TIVRMAQ *.gyrA.*.faa | grep -v AVYDTIVRMAQ > zyg_D87_mutants.txt

# select exactly S83L or D87N mutants
grep GKYHPHGDLAVY *.gyrA.*.faa > zyg_S83L_mutants.txt

grep AVYNTIVRMAQ *.gyrA.*.faa | grep -v AVYDTIVRMAQ > zyg_D87N_mutants.txt

# get full list
ls -1 *.gyrA.*.faa > zyg_all.nam

for i in *mutants.txt ; do sed -e 's/:.*//g' $i > `basename $i txt`nam ; done

# make combinations
comm  -12 zyg_S83L_mutants.nam zyg_D87N_mutants.nam > zyg_S83L+D87N_mutants.nam

cat zyg_S83L_mutants.nam zyg_D87N_mutants.nam | sort | uniq > zyg_S83L_D87N_any.nam

comm -23 zyg_D87N_mutants.nam zyg_S83L_mutants.nam > zyg_D87N_only.nam

comm -13 zyg_D87N_mutants.nam zyg_S83L_mutants.nam > zyg_S83L_only.nam

comm -23 zyg_all.nam zyg_S83L_D87N_any.nam > zyg_sensible.nam
comm -13 zyg_all.nam zyg_S83L_D87N_any.nam > zyg_resistant.nam

#echo -e "ID\tfilename\tS83mut\tD87mut\tS83L\tD87N\tresistant" > zyg_sensible.tsv
cat zyg_sensible.nam | \
	sed -e 's/^\(.*\)_/\1\t\1_/g' \
	    -e 's/$/\tN\tN\tN\tN\tN/g' \
         >> zyg_sensible.tsv

#echo -e "ID\tfilename\tS83mut\tD87mut\tS83L\tD87N\tresistant" > zyg_S83L_only.tsv
cat zyg_S83L_only.nam | \
	sed -e 's/^\(.*\)_/\1\t\1_/g' \
	    -e 's/$/\tY\tN\tY\tN\tY/g' \
         >> zyg_S83L_only.tsv

#echo -e "ID\tfilename\tS83mut\tD87mut\tS83L\tD87N\tresistant" > zyg_D87N_only.tsv
cat zyg_D87N_only.nam | \
	sed -e 's/^\(.*\)_/\1\t\1_/g' \
	    -e 's/$/\tN\tY\tN\tY\tY/g' \
         >> zyg_D87N_only.tsv

#echo -e "ID\tfilename\tS83mut\tD87mut\tS83L\tD87N\tresistant" > zyg_S83L+D87N_mutants.tsv
cat zyg_S83L+D87N_mutants.nam | \
	sed -e 's/^\(.*\)_/\1\t\1_/g' \
	    -e 's/$/\tY\tY\tY\tY\tY/g' \
         >> zyg_S83L+D87N_mutants.tsv


echo -e "ID\tfilename\tS83mut\tD87mut\tS83L\tD87N\tS83L.D87N\tresistant" > zyg_metadata.tsv
cat zyg_all.nam | \
    while read name ; do
        echo $name
        genome=`echo $name | sed -e 's/_.*//g'`
        RES='N'
        if grep -q $name zyg_D87_mutants.nam ; then D87='Y' ; else D87='N' ; fi
        if grep -q $name zyg_S83_mutants.nam ; then S83='Y' ; else S83='N' ; fi
        if grep -q $name zyg_D87N_mutants.nam ; then D87N='Y' ; RES='Y' ; else D87N='N' ; fi
        if grep -q $name zyg_S83L_mutants.nam ; then S83L='Y' ; RES='Y' ; else S83L='N' ; fi
        if grep -q $name zyg_S83L+D87N_mutants.nam ; then S83L_D87N='Y' ; RES='Y' ; else S83L_D87N='N' ; fi
        echo -e "$genome\t$name\t$S83\t$D87\t$S83L\t$D87N\t$S83L_D87N\t$RES" \
            >> zyg_metadata.tsv
done
