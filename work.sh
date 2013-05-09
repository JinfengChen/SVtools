ls `pwd`/../input/*.vcf | grep -v "cmp" | sed "s/@//" > vcflist
qsub -q js venn.sh

echo "largeSV"
ls `pwd`/../input/largeSV/*.vcf | grep -v "cmp" | sed "s/@//" > vcflist


