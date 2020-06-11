#To be run in workspace directory

for sim in ./*
do
echo $sim >> ../charges.txt
grep qtot $sim/init.top | tail -n 1 >> ../charges.txt
done

echo NONZERO CHARGES: >> ../charges.txt

grep -v 0.000000 ../charges.txt | grep -B 1 qtot >> ../charges.txt
