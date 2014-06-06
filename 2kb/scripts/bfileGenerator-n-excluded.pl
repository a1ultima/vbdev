#!/usr/bin/perl -w

#
# usage: ./this-prog.pl -order 2 input.fasta > output.bg2
#

#
# warning, doesn't seem to handle masked sequences at all??
#

use Getopt::Long;

my $order = 3;
my $skipids = '^CR';

GetOptions("order=i"=>\$order,
           "skipids=s"=>\$skipids,
          );

my %seqs;
my $id;
while (<>) {
  if (/^>(\S+)/) {
    $id = $1;
    undef $id if ($id =~ /$skipids/);
  } elsif ($id) {
    s/\W//;
    $seqs{$id} .= $_;
  }
}

foreach my $k (1 .. $order) {
  my %count;
  my $n = 0;
  foreach $id (keys %seqs) {
    my $seq = $seqs{$id};
    my $len = length($seq);
    for (my $i=0; $i<$len-($k-1); $i++) {
      my $ktuple = substr($seq, $i, $k);
      next if ($ktuple =~ /[^ACGT]/);
      $count{$ktuple}++;
      $n++;
    }
  }
  my $expected = 4**$k;
  if (keys %count < $expected) {
    die "not all $k-tuplets seen in input seqs\n";
  }
  print "# order = $k\n";
  foreach my $subseq (sort keys %count) {
    printf "%-10s %f\n", lc($subseq), $count{$subseq}/$n;
  }
}