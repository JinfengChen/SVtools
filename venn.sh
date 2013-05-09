#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=10:00:00

#export PERL5LIB = /opt/vcftools/0.1.8.1/lib/perl5/site_perl/
PATH=/rhome/cjinfeng/software/tools/tabix/tabix-0.2.6:$PATH

cd $PBS_O_WORKDIR

date

perl VCFvenn.pl --vcf vcflist

date

echo "Done"
