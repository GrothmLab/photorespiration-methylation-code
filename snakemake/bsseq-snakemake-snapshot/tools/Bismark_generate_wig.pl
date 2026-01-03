#!/usr/bin/perl

# file: Bismark_generate_wig.pl
# author: Zhenhui Zhong & Wanlu Liu, Jacobsen lab, UCLA 
# date: 11-July-2017

########################################################################################
#                                                                                       #
#   This is a perl script to generate wig files for CX_report.txt of Bismark(v0.18.0).  #
#   Site with CT count smaller than 4 will be filtered.                                 #                                
#   Input files: CX_report.txt                                                          #
#   Output files: 1. xxx_CG.wig 2. xxx_CHG.wig 3. xxx_CHH.wig                           #
#   Usage: perl Bismark_generate_wig.pl xxx.CX_report_sorted.txt(.gz) outdir            #
########################################################################################


my ($file)=$ARGV[0]=~m/^.*\/(\S+).CX_report_sorted.txt.*$/;
#print "$file\n";
if ($ARGV[0]=~m/\S+.CX_report_sorted.txt.gz$/) {
	open(IN, "gunzip -c $ARGV[0] |") or die "Usage: perl Bismark_generate_wig.pl xxx.CX_report.txt.gz outdir";
}
else {
	open IN,"$ARGV[0]" or die "Usage: perl Bismark_generate_wig.pl xxx.CX_report.txt outdir";
}
my ($outdir)=$ARGV[1]=~/(\S+)/;
open OUT1,">$outdir\/$file\_CG.wig" or die "Usage: perl Bismark_generate_wig.pl xxx.CX_report.txt outdir";
open OUT2,">$outdir\/$file\_CHG.wig" or die "Usage: perl Bismark_generate_wig.pl xxx.CX_report.txt outdir";
open OUT3,">$outdir\/$file\_CHH.wig" or die "Usage: perl Bismark_generate_wig.pl xxx.CX_report.txt outdir";
while(<IN>){
	chomp;
	my @a=split/\t/;
	$total_count = $a[3]+$a[4];
	next if($total_count <=3);
	$level = $a[3]/$total_count;
	$type=$a[5];
	$ratio=$a[2] eq "+"?$level:-1*$level;
	$sort_c{$a[0]."\t".$type}=$a[0];
	$sort_p{$a[0]."\t".$type}=$a[1];
	if(!$hash{$a[0]."\t".$type}){
		$hash{$a[0]."\t".$type}.="variableStep  chrom=$a[0]\n";
	}
	$hash{$a[0]."\t".$type}.="$a[1]	$ratio\n" if (abs($ratio) > 0);
}
for (sort {$sort_c{$a} cmp $sort_c{$b} || $sort_p{$a} <=> $sort_p{$b} } keys %hash){
	if($_=~/CG/){$out="OUT1";}
	if($_=~/CHG/){$out="OUT2";}
	if($_=~/CHH/){$out="OUT3";}
	print $out $hash{$_};
}
close IN;
close OUT1;
close OUT2;
close OUT3;
