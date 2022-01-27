#!/usr/bin/perl

my $inputWIG = $ARGV[0];
my $mappedReads = $ARGV[1];
my @lineArr;

my $outputWIG = $inputWIG;
$outputWIG =~ s/\.gz$//g;
$outputWIG =~ s/\.wig$//g;
$outputWIG .= ".RPM.wig";
$mappedReads = $mappedReads/1000000;
print $mappedReads . "\n";

my $cntr=0;
open(FD,"gunzip -c $inputWIG |");
open(FOUT,"> $outputWIG");



while(<FD>){
	chomp($_);
# 	if($_ =~ m/^(\d+\.{0,1}\d*)$/){
	if($_ !~ m/track/ && $_ !~ m/variable/) {
		@lineArr = split("\t", $_);
		my $height = $lineArr[1];
# 		print "h1".  $height . "\n";
		$height = $height/$mappedReads;
# 		$height =~ m/^(\d+\.{0,1}\d{0,4})/;
# 		$height = $1;
# 		print $height . "\n";
		print FOUT $lineArr[0] . "\t" . $height."\n";
	}else{
# 		print $_ . "\n";
		print FOUT $_."\n";
	}
	$cntr++;
#	if($cntr>10){last;}
}
close(FD);
close(FOUT);
my $res = `gzip $outputWIG`;
