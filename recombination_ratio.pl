use strict;

if (@ARGV ne 4) {die("Insufficient parameter!")}

open (TAB,"<$ARGV[0]");
open (DEPTH,"<$ARGV[1]");
open (IN,"<$ARGV[2]");
open (OUTT,">$ARGV[3]");
my $head;
my %fasta;
while (<IN>) {
	chomp;
	if (/>(\S+)/) {$head = $1}
	else {$fasta{$head} .= $_}
}

my %reads_number = ();
while (<DEPTH>) {
	chomp;
	my @tmp = split /\s+/,$_;
	if ($tmp[4] > 0) {
		if (($tmp[2] < 299) and (($tmp[2]+$tmp[4]+299) > (length($fasta{$tmp[1]})))) {
			$reads_number{$tmp[1]}++;
		}
	}
}

while(<TAB>){
	chomp;
	my @tmp = split /\s+/,$_;
	print OUTT $tmp[0],"\t",$tmp[1],"\t",$tmp[2],"\t",$tmp[3],"\t",
	$reads_number{$tmp[4]},"\t",$reads_number{$tmp[5]},"\t",$reads_number{$tmp[6]},"\t",$reads_number{$tmp[7]},"\n";
}

#Repeat_11	149	contig3|59736-59884	contig1|68832-68684	Repeat_11-2	Repeat_11-3	Repeat_11-1	Repeat_11-4