<a name="eKZn2"></a>
# 写在前面
> 最近在写毕业设计，要绘制降解组的图片，要求miRNA切割位点+保守结构域一起绘制（如下，从左变成右），问了师姐只能手搓出来，实在是懒，不想搓，想起python应该可以搞定，就是可能需要一些微调，遂捡起来看看。这其中还是发生了一些有趣的事情。折腾，人活着不就是得折腾。

<a name="JT12l"></a>
## 前期准备
<a name="ElxRh"></a>
#### 修改ClevaeLand4源代码
[GitHub - MikeAxtell/CleaveLand4: CleaveLand4: Analysis of degradome data to find sliced miRNA and siRNA targets](https://github.com/MikeAxtell/CleaveLand4)
> ClevaeLand4 只能生成 pdf，本来是想用 pdf 转换成 python 可以修改的格式，但是网上找了一圈都没有很好的办法。想到它绘图的源数据为什么没有存下来，于是翻了一下它的源代码，发现它绘制完图之后就删掉了绘图数据 (好恶心！) ，于是乎，篡改，拿到原数据

```perl
## 一般是在1142行，特可以直接搜索函数名,这是没改过的样子
sub make_t_plot {
    my($dir_name,$deg_data,$tx_len,$gline,$cat_digit,$pval) = @_;

    # open tmp file for data
    my $args1 = "Tplot_tmp.txt";
    (open(TMP, ">$args1")) || die "ABORT: Failed to open a temp file for Tplot data in sub-routine make_t_plot\n";
    
    # get required arguments
    my @gfields = split ("\t", $gline);
    my $args2 = "$dir_name" . "\/" . "$gfields[0]" . "_" . "$gfields[1]" . "_" . "$gfields[4]" . "_TPlot.pdf";
    my $args3 = "T\=$gfields[1]" . "_" . "Q\=$gfields[0]" . "_" . "S\=$gfields[4]";
    my $args4 = "category\=$cat_digit" . "_" . "p\=$pval";
    
    # print out the data to temp file
    print TMP "Position\tAll\tSite\n";
    my $i;
    my $junk;
    my @deg_lines = split ("\n", $deg_data);
    for($i = 1; $i <= $tx_len; ++$i) {
	print TMP "$i\t";
	if(@deg_lines) {
	    if($deg_lines[0] =~ /^$i\t(\d+)/) {
		print TMP "$1\t";
		if($i == $gfields[4]) {
		    print TMP "$1\n";
		} else {
		    print TMP "NA\n";
		}
		$junk = shift @deg_lines;
	    } else {
		print TMP "0\tNA\n";
	    }
	} else {
	    print TMP "0\tNA\n";
	}
    }
    close TMP;
    
    # Call R
    system "R --vanilla --quiet --slave --args $args1 $args2 $args3 $args4 < CleaveLand4_Tplotter.R 2> /dev/null";
    
    # clean up TMP
    system "rm -f $args1";
    
    # return the file path
    return $args2;
}
```
更改后
```perl
sub make_t_plot {
    my($dir_name,$deg_data,$tx_len,$gline,$cat_digit,$pval) = @_;

    # open tmp file for data
    #my $args1 = "Tplot_tmp.txt";
    my $args1 = "$gfields[0]" . "_" . "$gfields[1]" . "_" . "$gfields[4]" . "_TPlot_tmp.txt";
    (open(TMP, ">$args1")) || die "ABORT: Failed to open a temp file for Tplot data in sub-routine make_t_plot\n";
    
    # get required arguments
    my @gfields = split ("\t", $gline);
    #my $args2 = "$dir_name" . "\/" . "$gfields[0]" . "_" . "$gfields[1]" . "_" . "$gfields[4]" . "_TPlot.pdf";
    my $args2 = "$dir_name" . "\/" . "$gfields[0]" . "_" . "$gfields[1]" . "_" . "$gfields[4]" . "_TPlot.png";
    my $args3 = "T\=$gfields[1]" . "_" . "Q\=$gfields[0]" . "_" . "S\=$gfields[4]";
    my $args4 = "category\=$cat_digit" . "_" . "p\=$pval";
    
    # print out the data to temp file
    print TMP "Position\tAll\tSite\n";
    my $i;
    my $junk;
    my @deg_lines = split ("\n", $deg_data);
    for($i = 1; $i <= $tx_len; ++$i) {
	print TMP "$i\t";
	if(@deg_lines) {
	    if($deg_lines[0] =~ /^$i\t(\d+)/) {
		print TMP "$1\t";
		if($i == $gfields[4]) {
		    print TMP "$1\n";
		} else {
		    print TMP "NA\n";
		}
		$junk = shift @deg_lines;
	    } else {
		print TMP "0\tNA\n";
	    }
	} else {
	    print TMP "0\tNA\n";
	}
    }
    close TMP;
    
    # Call R
    system "R --vanilla --quiet --slave --args $args1 $args2 $args3 $args4 < CleaveLand4_Tplotter.R 2> /dev/null";
    
    # clean up TMP
    # system "rm -f $args1";
    
    # return the file path
    return $args2;
}
```
如此我们就获取了绘图的原始数据
<a name="WaCvt"></a>
#### 获取CDS序列上的结构与坐标
嗯，这个我是去CD-search上一个一个复制下来的(因为Betch CD--search 只支持蛋白序列)。当然也可以直接用蛋白搜索出的结构域坐标进行变换，前后乘上3就行。<br />格式如下<br />第一列是基因or转录本名称，与降解组分析中的名称要相同
```
LITCHI015235,2000,6000
LITCHI015235,8000,10000
LITCHI005435,2000,4000
LITCHI005435,4000,6000
```
到这里，前期的准备就结束了，马上开始画图！
```
