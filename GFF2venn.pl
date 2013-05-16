#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","project:s","weight","help");


my $help=<<USAGE;
perl $0 --gff a.gff,b.gff,c.gff,d.gff --weight
Use R package Vennerable here. This package can draw weighted venn diagram for 2 or 3 sets.
The resulting pdf file could be modified in AI, so can make a perfect figure.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||="all";
my @gff=split(",",$opt{gff});
my $gff=join(" ",@gff);
my $num=@gff;
`cat $gff | sort -k1,1 -k4n -k5n | awk '{print \$1"\t"\$4"\t"\$5"\t"\$2}' > $opt{project}.gff`;
`mergeBed -n -nms -i $opt{project}.gff > $opt{project}.gff.merge`;
my $feature=`cut -f 4 $opt{project}.gff | sort -u`;
my @feat=split("\n",$feature);
for(my $i=0;$i<@feat;$i++){
   $feat[$i]="\"$feat[$i]\"";
}
my $feat=join(",",@feat);
my $weight=readmerge("$opt{project}.gff.merge",$feature);
my $w= $opt{weight} ? "TRUE" : "FALSE";
my $plot= @feat < 4 ? "plot(Vcombo,doWeights=$w,type=\"circles\")\n" : "plot(Vcombo,type=\"ellipses\")\n";


my $cmd =<<R;
require(Vennerable)
pdf("$opt{project}.gff.merge.pdf")
Vcombo <- Venn(SetNames=c($feat), Weight= c($weight))
$plot
dev.off()
R

open OUT, ">$opt{project}.gff.merge.R" or die "$!";
     print OUT "$cmd\n";
close OUT;
`cat $opt{project}.gff.merge.R | R --vanilla --slave`;



####
sub readmerge
{
my ($file,$feature)=@_;
my @feat=split("\n",$feature);
my $seq="0" x @feat;
#print "0:$seq\n";
my %hash;
for(my $i=0;$i<@feat;$i++){
   my $s=$seq;
   substr($s,$i,1)=1; 
   $hash{$s}=0;
   #print "1:$s\n";
   for(my $j=$i+1;$j<@feat;$j++){
         my $s1=$s;
         substr($s1,$j,1)=1;
         $hash{$s1}=0;
         #print "2:$s1\n";
         for (my $k=$j+1;$k<@feat;$k++){
             my $s2=$s1;
             substr($s2,$k,1)=1;
             $hash{$s2}=0;
             #print "3:$s2\n"; 
             for (my $h=$k+1;$h<@feat;$h++){
                 my $s3=$s2;
                 substr($s3,$h,1)=1;
                 $hash{$s3}=0; 
                 #print "4:$s3\n";
             }
         }
   }
}
my %rank;
for(my $i=0;$i<@feat;$i++){
   $rank{$feat[$i]}=$i;
}
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    print "$_\n";
    my @f=split(";",$unit[3]);
    my $s=$seq;
    for(my $i=0;$i<@f;$i++){
       my $r=$rank{$f[$i]};
       substr($s,$r,1)=1;
    }
    $hash{$s}++;
}
close IN;
my @temp;
for(keys %hash){
   push @temp,"\`$_\`=$hash{$_}";
}
my $cmd=join(",",@temp);
print "$cmd\n";
return $cmd;
}
 
