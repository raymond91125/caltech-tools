#!/usr/bin/perl

# generate mapping of wbgenes to their possible names, with same names to different genes in reverse priority in @tables.  2015 15 16

use strict;
use warnings;
use DBI;
my $dbh = DBI->connect ( "dbi:Pg:dbname=testdb;host=131.215.52.76", "", "") or die "Cannot connect to database!\n";     # for remote access
# my $dbh = DBI->connect ( "dbi:Pg:dbname=wobrdb", "", "") or die "Cannot connect to database!\n";
my $result;

  my %geneNameToId; my %geneIdToName;
#   my @tables = qw( gin_locus );
  my @tables = qw( gin_wbgene gin_seqname gin_synonyms gin_locus );
#   my @tables = qw( gin_seqname gin_synonyms gin_locus );
  foreach my $table (@tables) {
    my $result = $dbh->prepare( "SELECT * FROM $table;" );
    $result->execute();
    while (my @row = $result->fetchrow()) {
      my $id                 = "WBGene" . $row[0];
      my $name               = $row[1];
      my ($lcname)           = lc($name);
      $geneIdToName{$id}     = $name;
#       $geneNameToId{$lcname} = $id; 
      $geneNameToId{$name}   = $id; } }
#   return (\%geneNameToId, \%geneIdToName);

my $outfile = '/home/azurebrd/cron/gin_names/gin_names.txt';
open (OUT, ">$outfile") or die "Cannot create $outfile : $!";
foreach my $name (sort { $geneNameToId{$a} cmp $geneNameToId{$b} } keys %geneNameToId) {
  my $id = $geneNameToId{$name};
  my $primary = '';
  if ($geneIdToName{$id} eq $name) { $primary = 'primary'; }
  print OUT qq($id\t$name\t$primary\n);
} # foreach my $name (sort keys %geneNameToId)
close (OUT) or die "Cannot close $outfile : $!";
