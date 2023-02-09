use strict;
use List::MoreUtils ':all';
use List::Util qw/max min/;

if (@ARGV ne 4) {die("Insufficient parameter!")}

open (TAB,"<$ARGV[0]") or die();
open (IN,"<$ARGV[1]") or die();
open (OUT,">$ARGV[2]") or die();
open (OUTT,">$ARGV[3]") or die();
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

my (%ref1,%ref2,%alt1,%alt2);
my @order = ();
my %fasta_all = %fasta;
foreach my $name (@name){
	my @repeat = grep {$_ =~ /^$name\t/} @table;
	my %fasta_N = %fasta;
	my $copy_number = scalar @repeat;
	my $count = 0;
	my @leading = ();
	my @tailing = ();
	my $core = "";
	@order = ();
	my $first = 0;
	if ($copy_number eq 1) {next}
	for (my $m = 0; $m < $copy_number; $m++) {
		my @tmp = split /\s+/,$repeat[$m];
		if ($tmp[1] >= 350) {next}
		if (substr($fasta_N{$tmp[5]},min($tmp[2],$tmp[3])-1,$tmp[1]) =~ /0/) {next}
		else{
			substr($fasta_N{$tmp[5]},min($tmp[2],$tmp[3])-1,$tmp[1]) = "0" x $tmp[1];
			my $string2 = substr($fasta_all{$tmp[5]},min($tmp[2],$tmp[3])-1,$tmp[1]);
			my $ratio = 0;
			if (my @N = $string2 =~ /0/gi) {$ratio = scalar(@N) / length($string2)}
			if (($ratio > 0)and($first eq 0)) {
				$first++;
				substr($fasta_all{$tmp[5]},min($tmp[2],$tmp[3])-1,$tmp[1]) = "0" x $tmp[1];
				push @order,$tmp[5]."|".$tmp[2]."-".$tmp[3];
				my ($leading,$tailing) = string($tmp[2],$tmp[3],$tmp[5]);
				push @leading,$leading;
				push @tailing,$tailing;
				if ($tmp[2] < $tmp[3]) {
					$core = substr($fasta{$tmp[5]},$tmp[2]-1,$tmp[3]-$tmp[2]+1);
				}else{
					$core = revcom(substr($fasta{$tmp[5]},$tmp[3]-1,$tmp[2]-$tmp[3]+1));
				}
			}
			elsif ($ratio eq 0){
				substr($fasta_all{$tmp[5]},min($tmp[2],$tmp[3])-1,$tmp[1]) = "0" x $tmp[1];
				push @order,$tmp[5]."|".$tmp[2]."-".$tmp[3];
				my ($leading,$tailing) = string($tmp[2],$tmp[3],$tmp[5]);
				push @leading,$leading;
				push @tailing,$tailing;
				if ($tmp[2] < $tmp[3]) {
					$core = substr($fasta{$tmp[5]},$tmp[2]-1,$tmp[3]-$tmp[2]+1);
				}else{
					$core = revcom(substr($fasta{$tmp[5]},$tmp[3]-1,$tmp[2]-$tmp[3]+1));
				}
			}
		}
	}
	if (scalar @leading >= 2) {
		conformation(\@leading,\@tailing,$core,$name);
	}
}

foreach my $i (keys %fasta_all){
	my $ratio = 0;
	if (my @N = $fasta_all{$i} =~ /0/gi) {$ratio = scalar(@N) / length($fasta_all{$i})}
	print $i,"\t",$ratio,"\n";
}

sub conformation {
	my ($tmp1,$tmp2,$core,$name) = @_;
	my @leading = @$tmp1;
	my @tailing = @$tmp2;
	my $count = 0;
	my $length = scalar @leading;
	my @combine = ();
	my %output = ();
	for (my $m = 0; $m < $length; $m++) {
		for (my $n = $m+1; $n < $length; $n++) {
			my $combine_name = $order[$m]."\t".$order[$n];
			push @combine,$combine_name;
			$ref1{$combine_name} = $leading[$m].$core.$tailing[$m];
			$output{$ref1{$combine_name}} = 1;
			$ref2{$combine_name} = $leading[$n].$core.$tailing[$n];
			$output{$ref2{$combine_name}} = 1;
			$alt1{$combine_name} = $leading[$m].$core.$tailing[$n];
			$output{$alt1{$combine_name}} = 1;
			$alt2{$combine_name} = $leading[$n].$core.$tailing[$m];
			$output{$alt2{$combine_name}} = 1;
		}
	}
	my %ref_name = ();
	my %ref_fasta = ();
	foreach my $i (sort keys %output){
		$count++;
		$ref_name{$i} = $name."-".$count;
		$ref_fasta{$name."-".$count} = $i;
	}
	my @redundency = ();
	my @nonredundency_ref = ();
	foreach my $i (@combine){
		my @uniq_sort = ();
		push @uniq_sort,$ref_name{$ref1{$i}};
		push @uniq_sort,$ref_name{$ref2{$i}};
		push @uniq_sort,$ref_name{$alt1{$i}};
		push @uniq_sort,$ref_name{$alt2{$i}};
		@uniq_sort = uniq(@uniq_sort);
		if (scalar @uniq_sort eq 4) {
			my $tmp = join " ",@uniq_sort;
			if (!grep {$_ eq $tmp} @redundency) {
				my $pos = $i;
				$pos =~ s/\|/\t/gi;
				$pos =~ s/\-/\t/gi;
				my @pos = split /\s+/,$pos;
				if ($pos[0] eq $pos[3]) {
					if ((min($pos[4],$pos[5])-max($pos[1],$pos[2]) > 50)or(min($pos[1],$pos[2])-max($pos[4],$pos[5]) > 50)) {
						push @redundency,$tmp;
						my ($ref1,$ref2) = ($ref_fasta{$ref_name{$ref1{$i}}},$ref_fasta{$ref_name{$ref2{$i}}});
						my $tolerate1 = substr($ref1,300-int((350 - length($core))/2),350);
						my $tolerate2 = substr($ref2,300-int((350 - length($core))/2),350);
						if (($ref2 =~ /$tolerate1/)or($ref1 =~ /$tolerate2/)) {
							substr($fasta_all{$pos[0]},min($pos[1],$pos[2])-1,max($pos[1],$pos[2])-min($pos[1],$pos[2])+1) = "1" x (max($pos[1],$pos[2])-min($pos[1],$pos[2])+1);
							substr($fasta_all{$pos[3]},min($pos[4],$pos[5])-1,max($pos[4],$pos[5])-min($pos[4],$pos[5])+1) = "1" x (max($pos[4],$pos[5])-min($pos[4],$pos[5])+1);
							next;
						}
						else{
						print OUTT $name,"\t",length $core,"\t",$i,"\t",$ref_name{$ref1{$i}},"\t",$ref_name{$ref2{$i}},"\t",$ref_name{$alt1{$i}},"\t",$ref_name{$alt2{$i}},"\n";
						push @nonredundency_ref,$ref_name{$ref1{$i}};
						push @nonredundency_ref,$ref_name{$ref2{$i}};
						push @nonredundency_ref,$ref_name{$alt1{$i}};
						push @nonredundency_ref,$ref_name{$alt2{$i}};
						}
					}
					else{
						substr($fasta_all{$pos[0]},min($pos[1],$pos[2])-1,max($pos[1],$pos[2])-min($pos[1],$pos[2])+1) = "1" x (max($pos[1],$pos[2])-min($pos[1],$pos[2])+1);
						substr($fasta_all{$pos[3]},min($pos[4],$pos[5])-1,max($pos[4],$pos[5])-min($pos[4],$pos[5])+1) = "1" x (max($pos[4],$pos[5])-min($pos[4],$pos[5])+1);
					}
				}
				else{
					print OUTT $name,"\t",length $core,"\t",$i,"\t",$ref_name{$ref1{$i}},"\t",$ref_name{$ref2{$i}},"\t",$ref_name{$alt1{$i}},"\t",$ref_name{$alt2{$i}},"\n";
					push @nonredundency_ref,$ref_name{$ref1{$i}};
					push @nonredundency_ref,$ref_name{$ref2{$i}};
					push @nonredundency_ref,$ref_name{$alt1{$i}};
					push @nonredundency_ref,$ref_name{$alt2{$i}};
				}
			}
		}
	}
	@nonredundency_ref = uniq(@nonredundency_ref);
	foreach my $i (@nonredundency_ref){
		print OUT ">",$i,"\n",$ref_fasta{$i},"\n";
	}
}

sub string {
	my ($start,$end,$contig) = @_;
	my ($leading,$tailing);
	my $length = length $fasta{$contig};
	if ($start < $end) {
		if ($start <= 300) {
			$leading = substr($fasta{$contig},$length-(300-$start)-2,300-$start+1).substr($fasta{$contig},0,$start-1);
		}else{
			$leading = substr($fasta{$contig},$start-301,300);
		}
		if ($end + 300 > $length) {
			$tailing = substr($fasta{$contig},$end,$length-$end).substr($fasta{$contig},0,300-($length-$end)+1);
		}else{
			$tailing = substr($fasta{$contig},$end,300);
		}
	}
	if ($start > $end) {
		if ($end <= 300) {
			$tailing = revcom(substr($fasta{$contig},$length-(300-$end)-2,300-$end+1)).revcom(substr($fasta{$contig},0,$end-1));
		}else{
			$tailing = revcom(substr($fasta{$contig},$end-301,300));
		}
		if ($start + 300 > $length) {
			$leading = revcom(substr($fasta{$contig},$start,$length-$start)).revcom(substr($fasta{$contig},0,300-($length-$start)+1));
		}else{
			$leading = revcom(substr($fasta{$contig},$start,300));
		}
	}
	return $leading,$tailing;
}

sub revcom {
	my $string = shift;
	my $revcom=reverse($string);
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	return $revcom;
}