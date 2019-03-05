#!/usr/bin/perl

use strict;
use warnings; 

use Bio::DB::Fasta;
use Bio::SeqIO;
use Data::Dumper;

my $usage = "extract_from_fasta_by_name.pl fasta_file query_file > STDOUT\n";
my $fasta = $ARGV[0] or die $usage;
my $query = $ARGV[1] or die $usage;

# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
                	        '-format' => "fasta");
                		        
# declare hash
my %queries;
my %queries_found;

#open query file, chomp line and save it into the array with push
#the hash will contain gene names (variable $line) as keys and 1 as a value (that's random)
open (IN , $query) or die "Can't open $query, $!\n";

while (my $line = <IN>){
	chomp $line;
	#print $line, "\n";
	$queries{$line} = 1;	
	}

#check structure of hash	
#print Dumper \%queries, "\n";

while (my $seq_obj = $seqio_obj->next_seq){

    my $seqname = $seq_obj->primary_id;
    my $description = $seq_obj->description;

    if ( exists ( $queries{$seqname} ) ) {

	# if the sequence is found, store it in new hash
	$queries_found{$seqname} = 1;

#     	print ">",  $seq_obj->description, "\n";
        print ">",  $seq_obj->primary_id, " ",  $seq_obj->description, "\n";
       	print $seq_obj->seq, "\n";
    }

}

# print out sequences that were not found
foreach my $q ( keys %queries ) {

    if ( !exists $queries_found{$q} ) {

	print STDERR "sequence $q not found!\n";

    }

}
