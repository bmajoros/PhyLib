#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <mpirun-*.out>\n" unless @ARGV==1;
my ($infile)=@ARGV;

while(1) {
  system("grep GRAD $infile");
  sleep(30);
}

