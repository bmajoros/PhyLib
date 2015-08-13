#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <in.newick>\n" unless @ARGV==1;
my ($infile)=@ARGV;

my $nextID=1;

open(IN,$infile) || die $infile;
while(<IN>) {
  if(/tree.*=\s*(\(.*\S)/) { convert($1) }
  elsif(/TREE:\s*(\(.*\S)/) { convert($1) }
}
close(IN);

sub convert
  {
    my ($newick)=@_;
    $newick=~s/\s+//g;
    my $tree=parse($newick);
    emitPhy($tree);
  }

sub generateName
  {
    my $id=$nextID;
    ++$nextID;
    return "A$id";
  }

sub emitPhy
  {
    my ($node)=@_;
    if($node->{type} eq "internal") {
      print "INTERNAL\n";
      my $name=$node->{name};
      if(!$name) { $name=generateName() }
      print "$name\n";
      my $dist=$node->{left}->{distance}+0;
      print "$dist\n";
      emitPhy($node->{left});
      $dist=$node->{right}->{distance}+0;
      print "$dist\n";
      emitPhy($node->{right});
    }
    else {
      print "LEAF\n";
      my $name=$node->{name};
      print "$name\n";
    }
  }


# ================= PARSING ROUTINES ================

sub makeInternalNode
  {
    my ($left,$right,$name,$distance)=@_;
    my $rec=
      {
       type=>"internal",
       left=>$left,
       right=>$right,
       name=>$name,
       distance=>$distance
      };
    return $rec;
  }

sub makeLeaf
  {
    my ($name,$distance)=@_;
    my $rec=
      {
       type=>"leaf",
       name=>$name,
       distance=>$distance
      };
    return $rec;
  }

sub parse
  {
    my ($newick)=@_;
    my $pos=0;
    return pp_tree($newick,\$pos);
  }

sub peek
  {
    my ($buf,$pos)=@_;
    return substr($buf,$$pos,1);
  }

sub match
  {
    my ($buf,$pos,$c)=@_;
    my $nextChar=substr($buf,$$pos,1);
    unless($nextChar eq $c)
      { die "Syntax error in tree, near: \"$nextChar\" (pos=$$pos, buf=\"$buf\")\n" }
    ++$$pos;
  }

sub pp_tree
  {
    my ($buf,$pos)=@_;
    if(peek($buf,$pos) eq "(") { return pp_pair($buf,$pos) }
    else { return pp_primary($buf,$pos) }
  }


sub pp_pair
  {
    my ($buf,$pos)=@_;
    match($buf,$pos,"(");
    my $left=pp_tree($buf,$pos);
    match($buf,$pos,",");
    my $right=pp_tree($buf,$pos);
    match($buf,$pos,")");
    my $name=pp_optName($buf,$pos);
    my $dist=pp_optDist($buf,$pos);
    return makeInternalNode($left,$right,$name,$dist);
  }

sub pp_name
  {
    my ($buf,$pos)=@_;
    my $c=peek($buf,$pos);
    unless($c=~/[a-zA-Z_]/)
      {die "Syntax error in tree, near \"$c\" (pos=$$pos, buf=\"$buf\")\n"}
    my $L=length($buf);
    my $name;
    for( ; $$pos<$L ; ++$$pos) {
      my $c=peek($buf,$pos);
      if($c=~/[a-zA-Z_\d]/) { $name.=$c }
      else { last }
    }
    return $name;
  }

sub pp_primary
  {
    my ($buf,$pos)=@_;
    my $name=pp_name($buf,$pos);
    my $dist=pp_optDist($buf,$pos);
    return makeLeaf($name,$dist);
  }


sub pp_optName
  {
    my ($buf,$pos)=@_;
    if(peek($buf,$pos)=~/[a-zA-Z_]/) { return pp_name($buf,$pos) }
  }


sub pp_optDist
  {
    my ($buf,$pos)=@_;
    if(peek($buf,$pos) eq ":") {
      match($buf,$pos,":");
      my $lexeme;
      my $L=length($buf);
      for( ; $$pos<$L ; ++$$pos) {
	my $c=peek($buf,$pos);
	if($c=~/[\d\.\-]/) { $lexeme.=$c }
	else { last }
      }
      return 0+$lexeme;
    }
  }





