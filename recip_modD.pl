
#
###################################################################################################
#                                                                                                 #
# This program tests for full orthology by reciprocal blast     UPDATE 9/22/15                    #
#                                                                                                 #
###################################################################################################

#!/usr/bin/perl
use strict ;

my $eCutoff		= 1e-25	; # only accept hits below this e-value cutoff - should be 1 for the first run!!!!

my $base 		= 'Djak' ;
my $execPath	= '/usr/bin/blastall' ;
my $outFile2	= $base . "Reciprocals$eCutoff.tab" ;

my @species  = ('Mant', 'Sine', 'Tanz', 'TAlp', 'Gbif', 'Chan', 'Gurn', 'Marm', 'Niva', 'Occi', 'Prav', 'Rain', 'Shas', 'Sier', 'Whit', 'Yezo', 'Yuas') ;	# excluding Djak
my %fullName = (	# these are the names of the BLAST databases
		'Djak' => 'djakonovi.fasta.clean.ready' ,
		'Mant' => 'Mantophasma.fasta.clean.ready' ,
		'Sine' => 'sinensis.fasta.clean.ready' ,
		'Tanz' => 'Tanzaniophasma_gb.fasta.clean.ready' ,
		'TAlp' => 'Trinity_Alps_Original.fasta.clean.ready' ,
		'GBif' => 'Trinity_Bifratrilecta_Mix.fasta.clean.ready' ,
		'Chan' => 'Trinity_Chandleri_Original.fasta.clean.ready' ,
		'Gurn' => 'Trinity_Gurneyi_Original.fasta.clean.ready' ,
		'Marm' => 'Trinity_Marmoreus_Original.fasta.clean.ready' ,
		'Niva' => 'Trinity_Nivalis_Original.fasta.clean.ready' ,
		'Prav' => 'Trinity_Pravdini_Original.fasta.clean.ready' ,
		'Rain' => 'Trinity_Rainier_Original.fasta.clean.ready' ,
		'Shas' => 'Trinity_Shasta_Original.fasta.clean.ready' ,
		'Sier' => 'Trinity_Sierrabuttes_Original.fasta.clean.ready' ,
		'Whit' => 'Trinity_Whitechuck_Original.fasta.clean.ready' ,
		'Yezo' => 'Trinity_Yezoensis_Original.fasta.clean.ready' ,
		'Yuas' => 'yuasai_gb.fasta.clean.ready' ,
		'Occi' => 'Trinity_Occidentalis_Original.fasta.clean.ready' ) ;

my $blastExec 	=  'blastall -p blastn -a 4 -v 1 -b 1 -I T -e $eCutoff';

my @species1 = @species ; # or: my @species1 = @species[0..2]

###################################################################################################
# take Locusta hits in each species from files, and then blast everybody's sequences against everybody else:
my %hits = () ;
my $base = 'Djak';
foreach my $s1 (@species1) {
    my $file = "($base)Djak$s1.tab";
    unless (open FILE, '>' .$file) {
	die "\nUnable to create $file\n";
    }
	# make lookup table query->hit for Locusta:
        %{$hits{Djak}{$s1}} = readBlast("($base)Djak$s1.tab") ;
	# blast against all:
	my $file = "($base)($s1)Acc.fa";
        unless (open FILE, '>' .$file) {
	    die "\nUnable to create $file\n";
	}
my $query 	  = "($base)($s1)Acc.fa" ;
	foreach my $s2 ('Djak', @species) {
		next if ($s1 eq $s2) ;
		my $outFile  = "$base$s1$s2.bls" ;
		unless (open FILE, '>' .$outFile) {
		    die "\nUnable to create $outFile\n";
		}
		my $outFileS = "$base$s1$s2.tab" ;
		unless (open FILE, '>' .$outFileS) {
		    die "\nUnable to create $outFileS\n";
		}
		print "blasting $fullName{$s1} -> $fullName{$s2} ...\n" ;
		system("$blastExec -i $query -o $outFile -d $fullName{$s2}") unless ((-f $outFile) || (-f $outFileS)) ; #skip if done before...
		processBLAST ($outFile, $outFileS) unless (-f $outFileS) ;
		# make lookup table query->hit (check for reciprocal  best hits with this later):
		%{$hits{$s1}{$s2}} = readBlast($outFileS) ;
	}
}

###################################################################################################
# check for reciprocal  best hits (all pairwise comparisons):
open (OUT, ">$outFile2") || die ;
print OUT "L_migratoria" ;
push @species, 'Djak' ;
foreach my $s (@species) {
	print OUT "\t$fullName{$s}" ;
}
foreach my $s1 (0..$#species) {
	foreach my $s2 ($s1+1 .. $#species) {
		print OUT "\t$species[$s1]$species[$s2]" ;
	}
}
print OUT "\n" ;

foreach my $t (keys %{$hits{Djak}{Djak}}) {
	next unless $hits{Djak}{Djak}{$t} ;
	$hits{Djak}{Djak}{$t} = $t ;	# to find back the Locusta gene itself
	print OUT $t ;
	foreach my $i1 (0..$#species) {
		my $h1 = $hits{Djak}{$species[$i1]}{$t} || 0 ;
		print OUT "\t$h1" ;
	}

	# check for pairwise reciprocity:
	foreach my $i1 (0..$#species) {
		my $s1 = $species[$i1] ;
		my $h1 = $hits{Djak}{$s1}{$t} ;
		foreach my $i2 ($i1+1 .. $#species) {
			my $s2 = $species[$i2] ;
			my $h2 = $hits{Djak}{$s2}{$t} ;
			my $s1_s2 = $h1 && $hits{$s1}{$s2}{$h1} ;
			my $s2_s1 = $h2 && $hits{$s2}{$s1}{$h2} ;
			my $rec = ($h1 && $h2 && ($s1_s2 eq $h2) && ($s2_s1 eq $h1)) || 0 ;
			print OUT "\t$rec" ;
		}
	}
	print OUT "\n" ;
}
close OUT ;

###################################################################################################

sub readBlast {
	my $file = $_[0] ;
	my %hash ;
	$hash{'-'} = 0 ;
	open (IN, $file) || die "couldn't open $file" ;
	<IN> ;
	while (<IN>) {
		chomp ;
		my @l = split ("\t", $_) ;
		next if ($eCutoff && ($l[4] > $eCutoff)) ;	## new line to exclude distant hits
		$hash{$l[0]} = ($l[2] ne '-') ? $l[2] : 0 ;
	}
	close IN ;
	return %hash ;
}

###################################################################################################
#                                                                                                 #
# Extract info from blast output files                                                            #
#                                                                                                 #
###################################################################################################

sub processBLAST {
	my $inFile = shift ;
	my $outFile = shift ;

	open (IN, $inFile) || die "couldn't open $inFile" ;
	open (OUT, ">$outFile") || die "couldn't open $outFile" ;

	print OUT "query	frame	hit	hitFull	Expect	identities	total	startQuery	endQuery\n" ;

	my $n = my $noHit = 0 ;

	while (<IN>) {
		last if /^Query=/ ;
	}
	my $x = $_ ;
	while (<IN>) {
		$x .= $_ ;
		next unless /^Query=/ ;
		$n++ ;
		print OUT process($x, $noHit) ;
		$x = $_ ;
	}
	close IN ;
	$n++ ;
	print OUT process($x, $noHit) ;
	close OUT ;
	return ($n, $noHit) ;
}

###################################################################################################
sub process {
	my $x = shift ;
	my $noHit = shift ;

	# extract info for this query:
	(my $query)	= ($x =~ /^Query=\s*(?:gi\|)?([^\s\|]+)/) ;
	(my $hit)		= ($x =~ /\n>([^\n]+)\n/) ;
	(my $expect)	= ($x =~ /\n Score =\s*[\d\.]+ bits.*Expect(?:\(2\))? =\s*([\d\.e\+\-]+)\s/) ;
	my ($identities, $total) = ($x =~ /\n Identities =\s*(\d+)\/(\d+) /) ;
	(my $frame) 	= ($x =~ /\n Frame =\s*([\+\-\d]+)\s*\n/) ;
	(my $start) 	= ($x =~ /(?:Frame =\s*[\+\-\d]+\s*|\)\n)\nQuery: (\d+) /) ;
	(my $end) 		= ($x =~ /\nQuery:\s+\d+\s+\S+\s+(\d+)\s[^:]*:[^:]*(Subset of the database|Score)/s) ;
	my $acc = '' ;
	if ($hit) {
		($acc) 	= ($hit =~ /^(?:gi\|)?([\w\-\.]+)/) ;
		$acc ||= '-' ;
		($expect = 1 . $expect) if ($expect && $expect =~ /^e/) ;
	} else {
		$noHit++ ;
		$acc = $hit = '-' ;
	}
	$frame	||= '-' ;
	$expect	||= '-' ;
	$total	||= '-' ;
	$start	||= '-' ;
	$end	||= '-' ;
	$identities	||= '-' ;
	return "$query	$frame	$acc	$hit	$expect	$identities	$total	$start	$end\n" ;
}
