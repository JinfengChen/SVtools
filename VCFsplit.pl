#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"vcf:s","help");


my $help=<<USAGE;
Split VCF file into SNP,indel(<50), SV(>50bp). Also could be used to filter duplicated or other types of fails in filters (only output PASS).
perl $0 --vcf
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

filter($opt{vcf});


##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HEG4
##Chr1    21547   UNION_BC_k31_var_28204  A       T       .       PASS    KMER=31;SVLEN=0;SVTYPE=SNP      GT:COV:GT_CONF  1/1:0,36:24.95
sub filter
{
my ($file)=@_;
my $head=$1 if ($file=~/(.*)\.vcf$/);
my $snp="$head.SNP.vcf";
my $sv ="$head.SV.vcf";
my $indel="$head.INDEL.vcf";
my %hash;
open OUT1,">$snp" or die "$!";
open OUT2,">$sv" or die "$!";
open OUT3,">$indel" or die "$!"; 
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    if ($_=~/^#/){
       print OUT1 "$_\n";
       print OUT2 "$_\n";
       print OUT3 "$_\n";
    }
    my @unit=split("\t",$_);
    my $line=$_;
    if ($unit[6] eq "PASS" and $unit[7]=~/SVTYPE\=SNP/){  ### $unit[7]=~/SVTYPE\=SNP/ including SNP and SNP_FROM_COMPLEX 
       print OUT1 "$line\n";
    }elsif ($unit[6] eq "PASS" and $unit[7]!~/SVTYPE\=PH_SNPS/){ ### $unit[7]!~/SVTYPE\=PH_SNPS/ excluding PHased SNP, this part was in decom.vcf
       if (length $unit[3] < 50 and length $unit[4] < 50){  
          print OUT3 "$line\n";
       }else{
          print OUT2 "$line\n";
       }
    }
}
close IN;
close OUT1;
close OUT2;
close OTU3;
return \%hash;
}
 
