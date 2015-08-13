#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.mod> <out.matrix> <out.phy>\n" unless @ARGV==3;
my ($infile,$outMatrix,$outPhy)=@ARGV;

my $order;
open(OUT,">$outMatrix") || die;
print OUT "NMER_RAW\n";
open(IN,$infile) || die $infile;
while(<IN>) {
  if(/ORDER:\s*(\d+)/) { 
    $order=$1; 
    print OUT "$order\n" 
  }
  elsif(/BACKGROUND:\s*(\S.*)/) { 
    my $n=4**($order+1);
    print OUT "$n\n$1\n"
  }
  elsif(/RATE_MAT:/) {
    my $numLines=4**($order+1);
    print OUT "$numLines $numLines\n";
    for(my $i=0 ; $i<$numLines ; ++$i) {
      $_=<IN>;
      print OUT "$_";
    }
  }
}
close(OUT);

system("~/evolution/PhyLib/newick-to-phy.pl $infile > $outPhy");




