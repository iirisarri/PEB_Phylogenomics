#!/usr/bin/perl
 
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

# concat_fasta.pl
# Iker Irisarri. October 2017
# Modification of concat_fasta.pl that prints out the corresponding partition file in RAxML style
# This script will concatenate the files in the order given in the command line
# requirements: all taxa in individual fasta files should have matching headers
# data type is hard-coded!!
# DNA will use "DNA," and "NNNNN" for missing data
# PROTS will use "LG," and "?????" for missing data

my $usage = "concat_fasta.pl infile1.fa infile2. etc. [Data type is hardcoded!]> outfile.fa partitionfile.part";
my @infiles = @ARGV or die $usage; 

# declare hashes
my %new;
my $order = 0; # traces order of genes
my %concat;
my %partitions;
my $total_concat_length = 0;
#my $datatype = "DNA";	# DNA
my $datatype = "LG";	# PROTS

# outfile to print out partition information
open (OUT, ">", "partitionfile.part");

# loop through the list of infiles
foreach my $infile (@infiles) {

	my $new_seq_length = 0;
	%new = ();
	$order++;
	
	# read file in with seqio
	my $seqio_obj = Bio::SeqIO->new(
		'-file' => "<$infile", 
		'-format' => "fasta", 
		);

	# store sequences of new file into hash %new
	while (my $inseq = $seqio_obj->next_seq) {

		$new{$inseq->primary_id}=$inseq->seq;
		$new_seq_length = $inseq->length;
		
		# save aln lengths into %partitions
		if ( !exists $partitions{$order} ) {
		
			$partitions{$order} = [$infile, $new_seq_length];
		}
	}


	# add sequences to %concat
	foreach my $k ( keys %new ) {

		# if not present in %concat
		if ( !exists $concat{$k} ) {
		
			# prior to concatenation, fill previous positions with "Ns" for $total_concat_length
			# for gene #1 this $total_concat_length == 0
#			my $seq1 = 'N' x $total_concat_length . $new{$k};
			my $seq1 = '?' x $total_concat_length . $new{$k};

			$concat{$k} = $seq1;
			
		}
		else {
	
			# concatenate
			my $seq2 = $concat{$k} . $new{$k};
			# reassign
			$concat{$k} = $seq2;
		}
	}
	# fill with Ns any taxa present in %concat but not in %new
	foreach my $j ( keys %concat ) {
	
		if ( !exists $new{$j} ) {
		
#			my $seq3 = $concat{$j} . 'N' x $new_seq_length;
			my $seq3 = $concat{$j} . '?' x $new_seq_length;
			$concat{$j} = $seq3;
		}
	}
	# add new gene length to total
	$total_concat_length += $new_seq_length;
}

# print out concatenated sequence
foreach my $taxa ( sort keys %concat ) {

	print ">$taxa\n";
	print $concat{$taxa}, "\n";
}

# print out partition file
my $prev_end = "";
foreach my $part ( sort keys %partitions ) {

	my $aln = @{ $partitions{$part} }[0];
	my $len = @{ $partitions{$part} }[1];
	
	# for first alignment
	if ( $part == 1 ) {
		print OUT "$datatype, $aln = 1-$len\n";
		
		#update $prev_end
		$prev_end = $len;
	}
	# for all other alignments
	else {

		my $start = $prev_end + 1;
		my $end = $prev_end + $len;
		print OUT "$datatype, $aln = $start-$end\n";
		
		# update $prev_end
		$prev_end = $end;
	}
}
close(OUT);

print STDERR "\ndone!\n\n";