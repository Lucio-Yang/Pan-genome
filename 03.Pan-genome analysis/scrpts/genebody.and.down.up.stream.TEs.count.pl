use strict;
use List::Util qw(sum);
use List::Util qw(min);
use POSIX;

my $fn_gene="sample.table.new";
my $fn_out=">sample.table.TEs.all.count.txt";
my $fn_ave=">sample.table.TEs.sum.200.20.count.txt";

my $fn_gff="sample.new.gff";
my $fn_len="sample.lens";
my $fn_tes="sample.target.bed";
my $ll=3; my $ini=0;

my $ud=2000; my $window=200; my $step=0.1;

## len file read
my %hs_len; open(LEN,$fn_len) or die $!;
while(<LEN>){
	chomp(); my @a=split(); $hs_len{$a[0]}=$a[1];
} close(LEN); print "$fn_len file done\n";

## gff file read
my %hs_gff; open(GFF,$fn_gff) or die $!;
while(<GFF>){
	chomp(); my @a=split();
	my $up=get_max($a[2]-$ud,1);                         #########上下游2000
	my $dw=get_min($a[3]+$ud,$hs_len{$a[0]});				#########上下游2000
	$hs_gff{$a[1]}={'c'=>$a[0],'s'=>$a[2],'e'=>$a[3],'u'=>$up,'d'=>$dw};
} close(GFF); print "$fn_gff file done\n";

open(OUT,$fn_out) or die $!; open(AVE,$fn_ave) or die $!;
my %hs_gene; my %hs_stat;
for(my $m=$ini;$m<$ll;$m++){
	## bias gene file read
	undef %hs_gene; open(GENE,$fn_gene) or die $!;
	while(<GENE>){
		#if($_!~/Om/){next;}     ########################
		chomp(); my @a=split();
		$hs_gene{$a[$m]}=[];
		for(my $i=$hs_gff{$a[$m]}{'s'}-$ud;$i<=$hs_gff{$a[$m]}{'e'}+$ud;$i++){  
			if($i<1 or $i>$hs_gff{$a[$m]}{'d'}){push(@{$hs_gene{$a[$m]}},-1);}
			else{push(@{$hs_gene{$a[$m]}},0);}
		}
		# print "$a[$m]\t$hs_gff{$a[$m]}{'s'}\t$hs_gff{$a[$m]}{'e'}\t$ud\n@{$hs_gene{$a[$m]}}\n"; last;
	} close(GENE); print "bias gene file $m done\n";

	## TEs file read and set 1 when overlap
	open(TE,$fn_tes) or die $!;
	while(<TE>){
		chomp(); my @a=split();
		my ($chr,$start,$end,$type)=($a[0],$a[1],$a[2],$a[4]);
		# if($type ne "LTR/Gypsy"){next;}
		foreach my $mk (sort{$a cmp $b} keys %hs_gene){
			if($chr eq $hs_gff{$mk}{'c'}){
				if(($hs_gff{$mk}{'u'}-$end)*($hs_gff{$mk}{'d'}-$start)<=0){
					for(my $i=$start;$i<=$end;$i++){
						my $x=$i-$hs_gff{$mk}{'s'}+$ud;
						my $y=$i-$hs_gff{$mk}{'e'}-$ud;
						if($y>0){last;}
						elsif($x>=0 and $y<=0){$hs_gene{$mk}[$x]=1;}
					}
				}
			}
		}
	}

	## calculater windows density and print into file
	my $ln=0; undef %hs_stat;
	foreach my $mk (sort{$a cmp $b} keys %hs_gene){
		my $s1=$hs_gff{$mk}{'s'}-$ud;
		my $s2=$hs_gff{$mk}{'u'};
		my $s3=$hs_gff{$mk}{'s'};
		my $e3=$hs_gff{$mk}{'e'};
		my $e2=$hs_gff{$mk}{'d'};
		my $e1=$hs_gff{$mk}{'e'}+$ud;
		print OUT "$mk\t$s1\t$s2\t$s3\t$e3\t$e2\t$e1";
		## 计算上游TEs在个窗口的密度
		sl_wd($mk,0,$ud-1,$window,$window*$step,0);
		## 计算基因区TEs在各窗口的密度
		my $st_g=ceil(($e3-$s3+1)/($ud/($window*$step)));
		if($st_g*(($ud-$window)/$window/$step)>=$e3-$s3+1){$st_g-=1;}
		my $wd_g=($e3-$s3+1)-$st_g*($ud-$window)/($window*$step);
		sl_wd($mk,$ud,$e3-$s1,$wd_g,$st_g,5000);
		## 计算下游TEs在各窗口的密度
		sl_wd($mk,$e3-$s1+1,$e1-$s1,$window,$window*$step,10000);
		print OUT "\n";
		$ln++; print "$ln gene done\n";
	}

	foreach my $mk(sort{$a <=> $b}keys %hs_stat){
		my $ave=$hs_stat{$mk}{'s'}/$hs_stat{$mk}{'n'};
		print AVE $ave."\t";
	} print AVE "\n";
} close(OUT); close(AVE);

sub sl_wd{
	my ($id,$is,$ie,$wd,$st,$fg)=@_; my @a=@{$hs_gene{$id}};
	print "$id,$is,$ie,$wd,$st\n";
	for(my $i=$is;$i+$wd-1<$ie+$st;$i+=$st){  ###后面的边界值注意超过了gene的区域
		my $ine=$i+$wd-1; if($ine>$ie){$ine=$ie;}
		# print "$i,$ine\t@a[$i..$ine]\n";
		my $sum=sum(@a[$i..$ine]); ### 判读有-1了怎么处理
		if(min(@a[$i..$ine])==-1){$sum=-1;print "Warning. there is a -1 in $i\n";}
		my $wl=$ine-$i+1;
		my $den=$sum/$wl;
		print OUT "\t".$den;
		my $order=1+($i-$is)/$st;
		if($den>=0){
			$hs_stat{$fg+$order}{'n'}++;
			$hs_stat{$fg+$order}{'s'}+=$den;}
	}
}
sub get_min{
	my $min = 1e15;
	foreach my $item (@_){
		if($item < 0){next;}
		if($item < $min){$min = $item;}
	}
	return $min;
}
sub get_max{
	my $max = -1;
	foreach my $item (@_){
		if($item < 0){next;}
		if($item > $max){$max = $item;}
	}
	return $max;
}
