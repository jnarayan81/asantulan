#!/usr/bin/perl

#use strict;
use 5.010;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Data::Dumper;
use Pod::Usage;
use Parallel::ForkManager;
use Bio::SeqIO;

#Author: Jitendra Narayan and Ashish K Shah
#Usage: perl asantulan.pl <file1 <file2> > <outfile>
#perl ld_all_ex_m.pl -f TRAK2_cpm.txt -t TRAK2_g.txt -c README.md -l 1

my (
	$cpm,
	$trakG,
	$config,
	$debug,
	$help,
	$man,
	$length,
	$tfile_corrected
);

my $version=0.1;
GetOptions(
    'cpm|f=s' => \$cpm,
    'trakG|t=s' => \$trakG,
    'cfile|c=s' => \$config,
    'length|l=i' => \$length,
    'help|h' => \$help
) or die &help($version);
&help($version) if $help;
#pod2usage("$0: \nI am afraid, no files given.")  if ((@ARGV == 0) && (-t STDIN));

if (!$cpm or !$trakG or !$config) { help($version) }
#if (!$thread) { $thread = `grep -c -P '^processor\\s+:' /proc/cpuinfo` }
#Keep all $parameters detail
my $parameters = readConfig ($config);
my $param = join (' ', @$parameters);

#Deal with file
my $fh= read_fh($cpm);
my @array2;
while(<$fh>) {
	chomp($_);
	push @array2, $_;
		}

my $c=0;
my $dnaf= read_fh($trakG);
while($c <= $#array2) {
	$final[$c]="\n".$array2[$c];
	@dna=<$dnaf>;
 close $dnaf;

 @nums;

 foreach $line(@dna){
	 print "\n$array2[$c]"."";
	 $c++ ;
	 $AA = 0; $TA = 0; $GA = 0; $CA = 0;
	 $AT = 0; $TT = 0; $GT = 0; $CT = 0;
	 $AG = 0; $TG = 0; $GG = 0; $CG = 0;
	 $AC = 0; $TC = 0; $GC = 0; $CC = 0; $NN = 0;
	 $e = 0;
	 $total = 0;
	 $freq_A = 0;
	 $freq_T = 0;
	 $freq_G = 0;
	 $freq_C = 0;
	 $t2 = 0;

	 $total= 0;

	 $GT_TG = 0;
	 $GC_CG = 0;
	 $AT_TA = 0;
	 $CA_AC = 0;
	 $GA_AG = 0;
	 $TC_CT = 0;

	 push @line_o, $line;

	 while($line =~ /AA/ig){$AA++}
	 while($line =~ /AC/ig){$AC++}
	 while($line =~ /AG/ig){$AG++}
	 while($line =~ /AT/ig){$AT++}

	 while($line =~ /TA/ig){$TA++}
	 while($line =~ /TC/ig){$TC++}
	 while($line =~ /TG/ig){$TG++}
	 while($line =~ /TT/ig){$TT++}

	 while($line =~ /GA/ig){$GA++}
	 while($line =~ /GC/ig){$GC++}
	 while($line =~ /GG/ig){$GG++}
	 while($line =~ /GT/ig){$GT++}

	 while($line =~ /CA/ig){$CA++}
	 while($line =~ /CC/ig){$CC++}
	 while($line =~ /CG/ig){$CG++}
	 while($line =~ /CT/ig){$CT++}
	 while($line =~ /NN/ig){$NN++}
	 while($line =~ /[^ATGC]/ig){$e++}

	 $geno_AA = 0;
	 $geno_AT_TA = 0;
	 $geno_AG_GA = 0;
	 $geno_AC_CA = 0;
	 $geno_TT = 0;
	 $geno_TG_GT = 0;
	 $geno_TC_CT = 0;
	 $geno_GG = 0;
	 $geno_GC_CG = 0;
	 $geno_CC = 0;
	 $freq_A = 0;
	 $freq_T = 0;
	 $freq_G = 0;
	 $freq_C = 0;
	 $t2 = 0;
	 $total=$AA+$AT+$AG+$AC+$CA+$CT+$CG+$CC+$GA+$GT+$GG+$GC+$TA+$TT+$TG+$TC ;

	 if ( $total > 0 ){
		 $geno_AA = 0;
		 $geno_AT_TA = 0;
		 $geno_AG_GA = 0;
		 $geno_AC_CA = 0;
		 $geno_TT= 0;
		 $geno_TG_GT = 0;
		 $geno_TC_CT = 0;
		 $geno_GG = 0;
		 $geno_GC_CG = 0;
		 $geno_CC= 0;

		 $GT_TG = $GT +$TG ;
		 $GC_CG = $GC+$CG ;
		 $AT_TA = $AT +$TA ;
		 $CA_AC = $CA +$AC ;
		 $GA_AG = $GA+$AG ;
		 $TC_CT =  $TC +$CT ;

		 $geno_AA = $AA/$total ;
		 $geno_AT_TA = ($AT+$TA)/$total ;
		 $geno_AG_GA = ($AG+$GA)/$total ;
		 $geno_AC_CA = ($AC+$CA)/$total ;
		 $geno_TT = $TT/$total ;
		 $geno_TG_GT = ($TG+$GT)/$total ;
		 $geno_TC_CT = ($TC+$CT)/$total;
		 $geno_GG = $GG/$total ;
		 $geno_GC_CG = ($GC+$CG)/$total;
		 $geno_CC = $CC/$total ;
		 $t2 = $total+$total;

		 $total_A = $AA+$AA+($AT+$AG+$AC+$CA+$GA+$TA);
		 $freq_A = $total_A/$t2;
		 $total_T = $TT+$TT+($TA+$TG+$TC+$AT+$GT+$CT);
		 $freq_T = $total_T/$t2;

		 $total_G = $GG+$GG+($GA+$GT+$GC+$AG+$TG+$CG) ;
		 $freq_G = $total_G/$t2;

		 $total_C = $CC+$CC+($CA+$CT+$CG+$AC+$TC+$GC) ;
		 $freq_C = $total_C/$t2 ;

		 $ef_AA = 0;
		 $ef_TT  = 0;
		 $ef_GG  = 0;
		 $ef_CC  = 0;

		 $ef_AT_TA = 0;
		 $ef_CA_AC = 0;
		 $ef_GA_AG = 0;
		 $ef_TC_CT = 0;
		 $ef_GC_CG = 0;
		 $ef_GT_TG = 0;

		 $ef_AA = $freq_A * $freq_A  ;
		 $ef_TT = $freq_T * $freq_T ;
		 $ef_GG = $freq_G * $freq_G ;
		 $ef_CC = $freq_C * $freq_C ;

		 $ef_AT_TA = ( $freq_A * $freq_T) + ($freq_T * $freq_A) ;
		 $ef_CA_AC = ($freq_C * $freq_A ) + ($freq_A * $freq_C) ;
		 $ef_GA_AG = ($freq_G * $freq_A) + ($freq_A * $freq_G) ;
		 $ef_TC_CT = ($freq_T * $freq_C ) + ($freq_C * $freq_T);
		 $ef_GC_CG = ($freq_G * $freq_C) + ($freq_C * $freq_G );
		 $ef_GT_TG = ( $freq_G * $freq_T) + ($freq_T * $freq_G );

		 $ld_AA = 0;
		 $ld_TT  = 0;
		 $ld_GG  = 0;
		 $ld_CC  = 0;

		 $ld_AT_TA = 0;
		 $ld_CA_AC = 0;
		 $ld_GA_AG = 0;
		 $ld_TC_CT = 0;
		 $ld_GC_CG = 0;
		 $ld_GT_TG = 0;

		 $ld_AA = $ef_AA - $geno_AA ;
		 $ld_TT = $ef_TT - $geno_TT ;
		 $ld_GG = $ef_GG - $geno_GG ;
		 $ld_CC = $ef_CC - $geno_CC ;

		 $ld_AT_TA = ( $ef_AT_TA - $geno_AT_TA ) ;
		 $ld_CA_AC = ( $ef_CA_AC - $geno_CA_AC ) ;
		 $ld_GA_AG = ( $ef_GA_AG - $geno_GA_AG ) ;
		 $ld_TC_CT = ( $ef_TC_CT - $geno_TC_CT ) ;
		 $ld_GC_CG = ( $ef_GC_CG - $geno_GC_CG ) ;
		 $ld_GT_TG = ( $ef_GT_TG - $geno_GT_TG ) ;

		 $ld_total = 0;
		 $ld_total = (
		 ( abs $ld_AA ) +
		 ( abs $ld_TT  ) +
		 ( abs $ld_GG  ) +
		 ( abs $ld_CC  ) +
		 ( abs $ld_AT_TA ) +
		 ( abs $ld_CA_AC )+
		 ( abs $ld_GA_AG ) +
		 ( abs $ld_TC_CT ) +
		 ( abs $ld_GC_CG ) +
		 ( abs $ld_GT_TG ) ) ;

		 print "\t$ld_total ";
	 }
	 else {
		 print "\t 0 " ;
	 }
 }
 exit ;
}

#Open and Read a file
sub read_fh {
	my $filename = shift @_;
	my $filehandle;
	if ($filename =~ /gz$/) {
		open $filehandle, "gunzip -dc $filename |" or die $!;
	}
	else {
		open $filehandle, "<$filename" or die $!;
	}
	return $filehandle;
}

#Read config files
sub readConfig {
my ($file) = @_;
my $fh= read_fh($file);
my @lines;
while (<$fh>) {
	chomp;
	next if /^#/;
	next if /^$/;
	$_ =~ s/^\s+|\s+$//g;
	push @lines, $_;
}
close $fh or die "Cannot close $file: $!";
return \@lines;
}

#Help section
sub help {
	my $ver = $_[0];
  print "\n asantulan $ver\n\n";

  print "Usage: $0 --cpm <> --trakG <> --cfile <> \n\n";
  print	"Options:\n";
  print "	--cpm|-q	infile file\n";
  print "	--trakG|-t	target file\n";
  print "	--cfile|-c	config file\n";
  print "	--length|-l	length below this is ignored\n";
  print " --help|-h	brief help message\n";

exit;
}
