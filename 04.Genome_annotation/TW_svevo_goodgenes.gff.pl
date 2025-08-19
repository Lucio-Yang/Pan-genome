#!perl -w
use strict;
my %goodGeneSet;
open INPUT, "Kronos_geneSet.list";
while( my $line=<INPUT> ){
	chomp $line;
	$goodGeneSet{$line}=1;
}
close INPUT;


open OUTPUT, ">Kronos_goodgenes.gff";
open INPUT, "/home/yangg/my_data/project/GeneAnnot/06.ab_initio_evidence/01_abinitio_prep/TW_svevo/transcriptome_homology_merge_evidence.gff";
while(my $line=<INPUT>){
	if ( $line=~/\smRNA\s.*ID=(\S+?);Parent=(\S+?);/ ){
		if(exists $goodGeneSet{$1}){
			print OUTPUT "$line";
		}
	}elsif ( $line=~/Parent=(\S+)$/ ){
		if(exists $goodGeneSet{$1}){
			print OUTPUT "$line";
		}
	}
}
close INPUT;
close OUTPUT;
