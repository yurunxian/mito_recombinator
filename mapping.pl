use strict;
use List::MoreUtils ':all';

open (IN,"<$ARGV[0]") or die();
open (TAB,"<$ARGV[1]") or die();
my $head;
my %fasta;
while (<IN>) {
	chomp;
	if (/>(\S+)/) {$head = $1}
	else {$fasta{$head} .= $_}
}
my @table = <TAB>;
my @name = ();
foreach my $i (@table){
	chomp($i);
	my @tmp = split /\s+/,$i;
	push @name,$tmp[0];
}
@name = uniq(@name);


foreach my $i (@name){
	print $i,"\n";
	open (OUT,">tmp.txt");
	my @seq = grep {$_ =~ /$i\-/} keys %fasta;
	foreach my $m (@seq){
		print OUT ">$m\n",$fasta{$m},"\n";
	}
	system("bowtie2-build tmp.txt tmp");
	system("bowtie2 -x tmp -1 $ARGV[2] -2 $ARGV[3] --end-to-end --no-discordant --no-mixed --very-fast -p 20 --no-unal -S tmp.sam");
	system("cat tmp.sam >> $ARGV[4]");
	system("rm tmp*");
}