#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"vcf:s","project:s","help");


my $help=<<USAGE;
Compare a list of vcf file and draw venn diagram
Note that the annotation in venn need to input by user in R code, according the order of output in venn.sh.o1094253, and run R by cat .R | R --slave
perl $0 --vcf vcflist
--vcf: list containing vcf file to be compared
echo vcflist
1.vcf
2.vcf
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}


$opt{project} ||="HEG4.LargeSV";
my $vcftool="/rhome/cjinfeng/software/tools/SVcaller/vcftools_0.1.10/bin";
my $tabix="/rhome/cjinfeng/software/tools/tabix/tabix-0.2.6";
my @vcf=readlist($opt{vcf});

#open OUT, ">vcfcompare.$opt{project}.summary0" or die "$!";
my @vcfgz;
for(my $i=0;$i<@vcf;$i++){
   my $title=$1 if ($vcf[$i]=~/(.*)\.vcf$/);
   print "$vcf[$i]\n$title\n";
   `cp $vcf[$i] $title.cmp.vcf` unless (-e "$title.cmp.vcf");
   `$tabix/bgzip $title.cmp.vcf` unless (-e "$title.cmp.vcf.gz");
   `$tabix/tabix -p vcf $title.cmp.vcf.gz` unless (-e "$title.cmp.vcf.gz.tbi");
   push @vcfgz, "$title.cmp.vcf.gz";
}
my $vcfcom=join(" ",@vcfgz);
print "$vcfcom\n";
`$vcftool/vcf-compare $vcfcom > vcfvenn.$opt{project}.raw` unless (-e "vcfvenn.$opt{project}.raw");
venndiagram("vcfvenn.$opt{project}.raw",\@vcfgz);
#`$vcftool/vcf-compare $vcf[$i].cmp.vcf.gz $refvcf.cmp.vcf.gz > $vcf[$i].cmp.vcf.gz.compare`;
#`$vcftool/vcf-compare $refvcf.cmp.vcf.gz $vcf[$i].cmp.vcf.gz | grep -v "#" > $vcf[$i].cmp.vcf.gz.compare`;
#my ($se,$sp)=readcom("$vcf[$i].cmp.vcf.gz.compare");
#print OUT "$depth\t$se\t$sp\n";
#`rm $vcf[$i].cmp.vcf.gz $vcf[$i].cmp.vcf.gz.tbi`; 
#`rm $refvcf.cmp.vcf.gz $refvcf.cmp.vcf.gz.tbi`;
#close OUT;

#`cat vcfcompare.$opt{project}.summary0 | sort -n > vcfcompare.$opt{project}.summary`;
#`rm vcfcompare.$opt{project}.summary0`;

#& drawcom("vcfcompare.$opt{project}.summary");


#############################
sub venndiagram
{
my ($file,$vcfgz)=@_;

####use word A/B/C to replace file name
my %name;
my %name2;
my @word=("1","2","3","4","5","6","7","8");
for(my $i=0;$i<@$vcfgz;$i++){
   $name{$1}=$word[$i] if ( $vcfgz[$i]=~ /\/([\w+|\.]*?\.vcf\.gz)/);
   print "$1\t$word[$i]\n";
   $name2{$word[$i]}=$1;
}

####make data for venn diagram
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    if ($_=~/^VN/){
        my @unit=split("\t",$_);
        print "$unit[1]\t";
        my @seq;
        for(my $i=2;$i<@unit;$i++){
           if ($unit[$i]=~ /\/([\w+|\.]*?\.vcf\.gz)/ ){
              print "$1\t$name{$1}\t";
              push @seq,$name{$1};
           }
        }
        print "\n";
        my $sequence=join("",sort @seq);
        $hash{$sequence}=$unit[1] unless (exists $hash{$sequence}); ### only need all number of each area, not specific of each area
    }
}
close IN;
###change the intersec to format of venn diagram in R. vcf-compare produce specific number in each area (like n12 have only shared by 1 and 2
###but R need all that shared between 1 and 2.
my @intersec;
my %hash2; ## store new intersec data
foreach my $s (sort {$a <=> $b} keys %hash){
    print "$s\t$hash{$s}\n";
    foreach my $k (keys %hash){
        if (matchnum($k,$s)){
           $hash2{$s}+=$hash{$k};
           print "$s\t$hash2{$s}\t$k:$hash{$k}\n";
        }else{
           print "$s:$k not match\n";
        }
         
    }
}
foreach my $s (sort {$a <=> $b} keys %hash2){
    print "$s\t$hash2{$s}\n";
    push @intersec, $hash2{$s};
}
#### draw venn diagram
drawvenn(\@intersec,\%name2);
}

###match two number, like 123 match 12 and 13 and 23, 1,2,3 etc
sub matchnum
{
my ($ref,$qry)=@_;
my @array=split("",$qry);
my $count;
for(my $i=0;$i<@array;$i++){
   if ($ref=~/$array[$i]/){
      $count++;
   }
}
my $flag=0;
if ($count == @array){
   $flag=1;
}
return $flag;
}


sub drawvenn
{
my ($intersec,$name)=@_;
my $n=keys %$name;
my $num;
if ($n == 3){ ### when using venn.triple, the order is n12, n23, n13 not n12, n13, n23.
   my $tmp=$intersec->[4];
   $intersec->[4]=$intersec->[5];
   $intersec->[5]=$tmp;
   $num=join(",",@$intersec);
}else{
   $num=join(",",@$intersec);
}
print "$num\n";


my $cmd =<<R;
require(VennDiagram)
pdf("vcfvenn.$opt{project}.pdf")
colortable <- c("orange", "red", "green", "blue","yellow")
color <- colortable[1:$n]
nametable <- c("first", "second", "third", "fouth","fivth")
name <- nametable[1:$n]

if ($n == 4) {
venn.plot <- draw.quad.venn($num,category=name,fill=color,lty="blank",alpha=rep(0.3,$n),margin=0.2)
}
if ($n == 3){
venn.plot <- draw.triple.venn($num,category=name,fill=color,lty="blank",alpha=rep(0.3,$n),margin=0.2)
}
if ($n == 2) {
venn.plot <- draw.pairwise.venn($num,category=name,fill=color,lty="blank",alpha=rep(0.3,$n),margin=0.2)
}
grid.draw(venn.plot)
dev.off()
R
open OUT, ">vcfvenn.$opt{project}.R" or die "$!";
     print OUT "$cmd\n";
close OUT;
`cat vcfvenn.$opt{project}.R | R --vanilla --slave`;
}




#VN	18381	../../input/HEG4_dbSNP.vcf.cmp.vcf.gz (18.5%)	./SNPcall.1X.gatk.calls.cmp.vcf.gz (43.4%edcom
sub readcom
{
my ($file)=@_;
open IN, "$file" or die "$!";
    my $line1=<IN>;
    my $line2=<IN>;
    my $line3=<IN>;
    my $temp=length $line1 > length $line2 ? $line1 : $line2;
    $max = length $temp > length $line3 ? $temp : $line3;
    #my @unit=split("\t",$first);
    my $sensity;
    my $specify;
    if ($max=~/.*\((.*?)\%\).*\((.*?)\%\)/){
       $sensity=$1;
       $specify=$2;
    }
close IN;
return ($sensity,$specify);
}

sub readlist
{
my ($file)=@_;
my @vcf;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    push @vcf, $_;
}
close IN;
return @vcf;
}
 
