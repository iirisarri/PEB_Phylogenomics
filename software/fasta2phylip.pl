#!/usr/bin/perl -w

use strict;
$|= 1;

my $fasta_file = $ARGV[0];
my %seqs;
my $len = 0;
my $max_name_len = 0;
my $name;
open(I,$fasta_file) || die "Unable to open fasta file $fasta_file\n";
while(<I>) {
	if(/^>/) {
		s/^>//;	chomp;
		#$name = $_;
		$name = (split)[0];
		$seqs{$name} = "";
	} else {
		s/[\n\t \[\]0-9]//g;
		$seqs{$name} .= $_;
	}
}
close(I);

foreach my $name ( keys %seqs ) {
	my $tmplen = length($seqs{$name});
	if($len == 0) {
		$len = $tmplen;
	} elsif($len != $tmplen) {
		print STDERR "Length for $name is $tmplen, expected $len\n";
		exit;
	}
	my $tmp_name_len = length($name);
	if($tmp_name_len > $max_name_len) {
		$max_name_len = $tmp_name_len;
	}
}

my $num_taxa = scalar(keys(%seqs));

print " $num_taxa $len\n";
foreach my $name ( keys %seqs ) {
	my $orig = $name;
	if(length($name) > 69) {
		$name = substr($name,0,69);
		$max_name_len = 70;
	}
	print $name.' 'x($max_name_len-length($name))."  $seqs{$orig}\n";
}


