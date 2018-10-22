#!/usr/bin/perl

#############################################
#A Perl script to create a .rmap and a .baitmap
#file from a HiCUP digest file and a list of oligo positions
#Oligo format (tab-delimited):
#Chromsome    Start    End    Feature    
#############################################

use strict;
use warnings;
use POSIX;

use Data::Dumper;

#Read in digest file
my ($digest_file, $oligo_file) = @ARGV;

die "Specify i) digest file and ii) oligos file\n" unless( scalar @ARGV == 2);
 
if ( $digest_file =~ /.*\.gz$/ ) {
    open( DIGEST, "zcat $digest_file |" ) or die "Cannot open '$digest_file' : $!";
} else {
    open( DIGEST, $digest_file ) or die "Cannot open file: '$digest_file' : $!";
}

my %digest_fragments;    # %{csome \t ten_kb_region}{"first_base\t$last_base"} = ''
scalar <DIGEST>;    #Ignore headers
scalar <DIGEST>;
while (<DIGEST>) {

	my $chromosome_name            = ( split /\t/ )[0];
	my $first_base                 = ( split /\t/ )[1];
	my $last_base                  = ( split /\t/ )[2];
	#my $fragment_number            = ( split /\t/ )[3];
	my $ten_kb_region              = ceil( $first_base / 10000 );
	my $fragment_end_ten_kb_region = ceil( $last_base / 10000 );

	do {
	    $digest_fragments{"$chromosome_name\t$ten_kb_region"}{"$first_base\t$last_base"} = '';
	    $ten_kb_region++;
	} while ( $ten_kb_region <= $fragment_end_ten_kb_region );
    
}
close DIGEST or die "Could not close filehandle on '$digest_file' : $! ";

#print Dumper \%digest_fragments;





#Read in oligos file and identify baited restriction fragment
if ( $oligo_file =~ /.*\.gz$/ ) {
    open( OLIGO, "zcat $oligo_file |" ) or die "Cannot open '$oligo_file' : $!";
} else {
    open( OLIGO, $oligo_file ) or die "Cannot open file: '$oligo_file' : $!";
}

while(<OLIGO>){
	my $line = $_;
	chomp $line;

	my ($csome, $start, $end, undef, $feature) = split(/\t/, $line);
	$csome =~ s/^chr//;
	add_feature($csome, $start, $end, $feature);
}
close OLIGO or die "Could not close filehandle on '$oligo_file' : $!";

#print Dumper \%digest_fragments;
#exit;

#Read in digest file again and print out results to file and reported the 
#number of baited fragments
if ( $digest_file =~ /.*\.gz$/ ) {
    open( DIGEST2, "zcat $digest_file |" ) or die "Cannot open '$digest_file' : $!";
} else {
    open( DIGEST2, $digest_file ) or die "Cannot open file: '$digest_file' : $!";
}

my $rmap_file = "$digest_file.rmap";
my $baitmap_file = "$digest_file.baitmap";
open(RMAP, '>', $rmap_file) or die "Couldn't write to '$rmap_file' : $!";
open( BAITMAP, '>', $baitmap_file) or die "Couldn't read '$baitmap_file' : $!";

scalar <DIGEST2>;    #Ignore headers
scalar <DIGEST2>;
my $index = 1;


#print Dumper \%digest_fragments;

while (<DIGEST2>) {

	my $chromosome_name            = ( split /\t/ )[0];
	my $first_base                 = ( split /\t/ )[1];
	my $last_base                  = ( split /\t/ )[2];
	#my $fragment_number            = ( split /\t/ )[3];
	my $ten_kb_region              = ceil( $first_base / 10000 );

	print RMAP "$chromosome_name\t$first_base\t$last_base\t$index\n";

	if($digest_fragments{"$chromosome_name\t$ten_kb_region"}->{"$first_base\t$last_base"} ne ''){		
		my $feature = $digest_fragments{"$chromosome_name\t$ten_kb_region"}->{"$first_base\t$last_base"};
		print BAITMAP "$chromosome_name\t$first_base\t$last_base\t$index\t$feature\n";
	}
	$index++;
}

close DIGEST2 or die "Could not close filehandle on '$digest_file' : $! ";
close RMAP or die "Could not close filehandle on '$rmap_file' : $! ";
close BAITMAP or die "Could not close filehandle on '$baitmap_file' : $! ";

print "Done";

exit (0);


###################################################
#Subrotine
###################################################

#Subroutine: add_feature
#Adds feature id to %digest_fragments (declared outside of subroutine)
sub add_feature{

	my ($csome, $start, $end, $feature) = @_;
	#print join("\t", $csome, $start, $end, $feature);
	#print "\n";
	my $pos = int( ($end - $start) / 2 ) + $start;
	my $pos_ten_kb_region = ceil( $pos / 10_000 );

	#Find the fragment
	unless(exists $digest_fragments{"$csome\t$pos_ten_kb_region"}) {
		die "Could not find '$csome\t$pos_ten_kb_region' in digest_fragments hash - the digest file and oligo probes file appear to be incompatible.\n";   
	}


	my ($start_frag, $end_frag);
	foreach ( keys %{ $digest_fragments{"$csome\t$pos_ten_kb_region"} } ) {
	    my $lookup_start_end_site = $_;                               
	    ($start_frag, $end_frag) = split(/\t/, $lookup_start_end_site);
	    
	    if ( ( $start_frag <= $pos ) and ( $end_frag >= $pos ) ) {
	    	last; 
	    }
	}



	#Now populate the hash
	my $start_ten_kb_region = ceil( $start_frag / 10_000 );
	my $end_ten_kb_region = ceil( $end_frag / 10_000 );


	unless(exists $digest_fragments{"$csome\t$start_ten_kb_region"}) {
		die "Could not find '$csome\t$start_ten_kb_region' in digest_fragments hash - the digest file and oligo probes file appear to be incompatible.\n";   
	}

	unless(exists $digest_fragments{"$csome\t$end_ten_kb_region"}) {
		die "Could not find '$csome\t$end_ten_kb_region' in digest_fragments hash - the digest file and oligo probes file appear to be incompatible.\n";   
	}

	for(my $ten_kb_region = $start_ten_kb_region; $ten_kb_region <= $end_ten_kb_region; $ten_kb_region++){    #Ensure all regions overlapping the frament have been included
	    
		#print "$ten_kb_region\n";
	    foreach ( keys %{ $digest_fragments{"$csome\t$ten_kb_region"} } ) {
	    	my $lookup_start_end_site = $_;                               
	    	my ($lookup_start, $lookup_end) = split(/\t/, $lookup_start_end_site);
	    	

	    	#print "$lookup_start_end_site\n";

	    	if ( ( $lookup_start <= $pos ) and ( $lookup_end >= $pos ) ) {
	    		$digest_fragments{"$csome\t$ten_kb_region"}->{$lookup_start_end_site} = $feature;
	    		#last;
	    	}
	    }
	}
}









