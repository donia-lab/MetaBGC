#!/usr/bin/perl -w
use strict;

my ($seq, $list, $matches, %list_hash);
my $switch = "0";

if ( @ARGV == 3 ) { $seq = $ARGV[0]; $list = $ARGV[1]; $matches = $ARGV[2] }
else { die "Usage:$0 fasta_file list_file matches_file\n"; }

open(LIST,"<$list") or die "Could not open $list!\n";

while ( <LIST> ) {
    if ( $_ =~ /^(\S+)/ ) { $list_hash{$1} = $_; }
    elsif ( $_ =~ /^(.+)\n/ ) { $list_hash{$1} = $_; }
    else { print STDERR "Problem: Can't read $_!\n"; }
}
close LIST;

open(SEQ,"<$seq") or die "Could not open $seq!\n";
open(OUT,">$matches") or die "Could not open $matches\n"; 

while ( <SEQ> ) {
    
    if ( $_ =~ /^>(\S+)/ ) {
	if ( exists($list_hash{$1}) ) { $switch = "1"; }
	else { $switch = "0"; }}
	
    if ( $switch eq "1" ) { print OUT $_; }
}
	
close SEQ;
close OUT;

print STDERR "Result file has been written to $matches\n"; 
