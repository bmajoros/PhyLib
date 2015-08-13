#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.mod> <out.matrix> <out.phy>\n" unless @ARGV==3;
my ($infile,$outMatrix,$outPhy)=@ARGV;

open(OUT,">$outMatrix") || die $outMatrix;
print OUT "COND\n0\n0\nDNA\nREV\n0\n1\nREV\n";
open(IN,$infile) || die $infile;
my ($alpha,$beta,$kappa,$chi,$omega,$tau);
my $tree;
while(<IN>) {
  if(/RATE_MAT:/) {
    $_=<IN>; my @fields=split/\s+/,$_; shift @fields;
    $beta=$fields[1]; $alpha=$fields[2]; $chi=$fields[3];

    $_=<IN>; @fields=split/\s+/,$_; shift @fields;
    $kappa=$fields[2]; $omega=$fields[3];

    $_=<IN>; @fields=split/\s+/,$_; shift @fields;
    $tau=$fields[3];

    print OUT "$alpha\t$beta\t$kappa\t$chi\t$omega\t$tau\n";
    print OUT "7\n1 1 1 1 0 0 0\n";
  }
  #elsif(/TREE:\s*(\S+);/) {$tree=$1}
}
close(IN);
close(OUT);

system("~/evolution/PhyLib/newick-to-phy.pl $infile > $outPhy");


