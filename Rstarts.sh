# edit zygRstats.R to set data.dir to 0050 and run it
$EDIT zygRstats.R
Rscript zygRstats.R

# edit zygRstats.R to set data.dir to 0250 and run it
$EDIT zygRstats.R
Rscript zygRstats.R

cd 0050.Rstats
grep -e '\(Glass\|SUMMARY\)' R.Rstats.log > effect.sum.txt
cd ..

cd 0250.Rstats
grep -e '\(Glass\|SUMMARY\)' R.Rstats.log > effect.sum.txt
cd ..
