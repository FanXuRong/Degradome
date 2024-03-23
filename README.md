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
```python
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.patches as patches
import re
import os


class Dragome:
    @staticmethod
    # 读取 ClevaeLand 输出的信息文件以及整理的含保守结构域的 csv 文件
    def Read_result_file (result_file,Conserved_file):
        
        Outfile = result_file +".csv"
        print (f'Inputfile is {result_file} {Conserved_file}\nOutfile is {Outfile}')
        ## 读取含保守结构域的 csv 文件
        Conserved_df = pd.read_csv(Conserved_file,header=None,names=['Conserved_Gene','Conserved_start','Conserved_end'])
        
        ## 读取 ClevaeLand 输出的信息文件
        with open(result_file, 'r') as f:
            miRNA_pos = []  # 靶位点比对信息的列表
            Category = 0  # 存储 Category 行
            T_plot_message_file = ""  # 存储 T-Plot file行
            score = 0  # 存储分数
            Message = []  # 存储所有信息
            for line in f:
                if line.startswith("5'"):
                    ## 检测到 miRNA 切割信息,开始处理上一段的切割信息,并整理其他信息到列表中
                    if miRNA_pos:
                        miRNA_target = ''
                        Tr_seq = re.findall(r'(.*?)Transcript: (.*?):\d+-\d+', miRNA_pos[0])[0]
                        miRNA_pos[0] = Tr_seq[0]+Tr_seq[1]+'\n'
                        miRNA_seq = re.findall(r'(.*?)Query: (.*?)$', miRNA_pos[2])[0]
                        miRNA_pos[2] = miRNA_seq[0] + miRNA_seq[1]
                        miRNA_pos[2] = miRNA_seq[0]+miRNA_seq[1]
                        miRNA_target = miRNA_pos[0] + miRNA_pos[1] + miRNA_pos[2]
                        # print (miRNA_target)
                        Message.append([gene,miRNA,T_plot_message_file,score,Category,p_value,miRNA_target,T_plot_fig_message_file])
                    miRNA_pos = []  # 将前三行列表重置为空列表
                    miRNA_pos.append(line)  # 将当前行添加到前三行列表中

                elif miRNA_pos and len(miRNA_pos) < 3:  # 如果前三行列表非空且长度小于 3 ，则继续添加行至前三行列表中，获取完整的 miRNA 靶向序列
                    miRNA_pos.append(line)
                    
                    ## 提取 score
                elif (re.findall(r'Allen et al. score: (\d+(\.\d+)?)',line)):
                    score = float(re.findall(r'Allen et al. score: (\d+(\.\d+)?)',line)[0][0])
                    # print (score)
                    
                    ## 提取 Category
                elif (re.findall(r'Degardome Category: (\d?)',line)):
                    Category = int((re.findall(r'Degardome Category: (.\d?)',line)[0]))
                    # print (Category)
                    
                    ## 提取 p 值
                elif (re.findall(r'Degradome p-value: (.*\d.\d?)',line)):
                    p_value = float(re.findall(r'Degradome p-value: (.*\d.\d?)',line)[0])
                    p_value = round(p_value, 6)
                    #print (p_value)
                    
                    ## 提取绘图的数据文件,参考 CleaveLand4.pl 中 make_t_plot 函数进行修改
                elif (line.startswith("T-Plot file:")):
                    a = re.findall(r'(.*?).png',line)[0]
                    T_plot_message_file = os.path.basename(a) + '_tmp.txt'
                    T_plot_fig_message_file = os.path.basename(a) + '.svg'
                    # print (T_plot_message_file)
                    
                    ## 提取基因以及miRNA信息
                elif (line.startswith("SiteID:")):
                    _,gene,miRNA = line.strip().split(':')
                    miRNA = 'miR' + miRNA
                    gene = gene.strip()
            # 处理最后一个结构的前三行数据
            if miRNA_pos:
                Tr_seq = re.findall(r'(.*?)Transcript: (.*?):\d+-\d+', miRNA_pos[0])[0]
                miRNA_pos[0] = Tr_seq[0]+Tr_seq[1]+'\n'
                miRNA_seq = re.findall(r'(.*?)Query: (.*?)$', miRNA_pos[2])[0]
                miRNA_pos[2] = miRNA_seq[0] + miRNA_seq[1]
                miRNA_pos[2] = miRNA_seq[0]+miRNA_seq[1]
                miRNA_target = miRNA_pos[0] + miRNA_pos[1] + miRNA_pos[2]
            # print (miRNA_target)
            Message.append([gene,miRNA,T_plot_message_file,score,Category,p_value,miRNA_target,T_plot_fig_message_file])
        ## 生成信息矩阵
        Message_df = pd.DataFrame(Message)
        Message_df.columns=["gene","miRNA","T_plot_message_file","score","Category","p_value","miRNA_target","T_plot_fig_message_file"]
        
        ## 匹配保守结构域
        Conserved_gene = set(Conserved_df['Conserved_Gene'])
        for gene in Conserved_gene:
            Select_df = Conserved_df[Conserved_df['Conserved_Gene']==gene]
            for i in  range(len(Select_df)):
                Conserved_start = f'Conserved{i+1}_start'
                Conserved_end = f'Conserved{i+1}_end'
                Message_df.loc[:,Conserved_start] = Select_df.iloc[i]['Conserved_start']
                Message_df.loc[:,Conserved_end] = Select_df.iloc[i]['Conserved_end']
        Message_df.to_csv(Outfile,index=False)
        return Message_df
        '''
        这段代码定义了一个类 Dragome，其中包含一个静态方法 Read_result_file。这个方法用于读取 ClevaeLand 输出的信息文件和整理的含保守结构域的 csv 文件，并返回一个信息矩阵。

        首先，将输出文件名设置为 result_file + ".csv"，然后读取含保守结构域的 csv 文件到 Conserved_df 数据帧中。

        接下来，打开 ClevaeLand 输出的信息文件，并依次读取每一行。根据行的内容进行不同的处理：

        如果行以 "5'" 开头，表示检测到 miRNA 切割信息，开始处理上一段的切割信息，并将其他信息整理到列表 Message 中。
        如果行不为空且 miRNA_pos 列表的长度小于 3，将行添加到 miRNA_pos 列表中，以获取完整的 miRNA 靶向序列。
        如果行中包含 "Allen et al. score"，提取分数并存储到变量 score 中。
        如果行中包含 "Degardome Category"，提取 Category 并存储到变量 Category 中。
        如果行中包含 "Degradome p-value"，提取 p 值并存储到变量 p_value 中。
        如果行以 "T-Plot file:" 开头，提取绘图的数据文件，并修改文件名后缀为 ".txt" 和 ".svg"。
        如果行以 "SiteID:" 开头，提取基因和 miRNA 的信息。
        处理完最后一个结构的前三行数据后，将所有信息存储到 Message 列表中，并将其转换为名为 Message_df 的数据帧。
        
        接下来，通过匹配保守结构域，将保守结构域的起始和结束位置添加到 Message_df 数据帧中的相应列中。
        
        最后，将 Message_df 数据帧保存为 csv 文件，并返回该数据帧。
        '''
    @staticmethod
    ## 绘图
    def Draw(Series):
        # 读取整理的数据
        gene,miRNA,Draw_data_file,score,Category,p_value,miRNA_target,Save_fig,Conserved1_start,Conserved1_end,Conserved2_start,Conserved2_end = Series[0:12]
        data = pd.read_csv(Draw_data_file,sep='\t')
        ###################################################################################################
        ## 绘制 T-plot
        ax1 = plt.subplot2grid((14,12),(0,0),rowspan=8,colspan=6)
        ax1.plot(data['Position'], data['All'], color='black')
        # 获取切割位点 Site 所在的 Series
        df_notnan = data[~data['Site'].isna()]
        # 获取位点
        Site = df_notnan.iloc[0,0]
        # 获取位点切割分数
        Site_height = df_notnan.iloc[0,2]
        # 获取 y 轴长度
        y_len = plt.ylim()[1] - plt.ylim()[0]
        # 获取分数相对于 y 轴长度的比值，用于做图
        Site_len = float((Site_height-(plt.ylim()[0]))/(y_len))
        ax1.axvline(x=Site, color=(251/255,87/255,84/255), ymin=0,ymax=Site_len,linestyle='-')
        plt.xlabel('Transcript Position')
        plt.ylabel("Degradome 5' end Frequency")
        plt.title('Degradome Data')
        x_len = plt.xlim()[1] - plt.xlim()[0]
        ax1_annotation = 'Category=' + str(Category) + "\n" + 'P=' + str(p_value)
        plt.text(0.6,0.8,ax1_annotation,transform=ax1.transAxes,fontsize=12,font={'family':'Consolas'})
        print (x_len)
        ####################################################################################################
        ## 绘制保守结构域
        ax2 = plt.subplot2grid((14,12),(10,0),rowspan=1,colspan=6,sharex=ax1)
        ## 按照 CDS 长度绘制最下层的图
        rect1 = plt.Rectangle((0, 0), x_len, 0.5, facecolor=(240/255,90/255,217/255), alpha=0.5)
        ## 第一个结构域
        rect2 = plt.Rectangle((Conserved1_start, 0), (Conserved1_end-Conserved1_start), 0.5, facecolor=(255/255,217/255,48/255), alpha=0.5)
        ## 第二个结构域
        rect3 = plt.Rectangle((Conserved2_start, 0), (Conserved2_end-Conserved2_start), 0.5, facecolor=(255/255,217/255,48/255), alpha=0.5)
        ax2.add_patch(rect1)
        ax2.add_patch(rect2)
        ax2.add_patch(rect3)
        ## 隐藏坐标轴
        ax2.axis('off')
        ###################################################################################################
        ## 绘制miRNA的切割信息
        ax3 = plt.subplot2grid((14,12),(13,0),rowspan=1,colspan=6)
        ax3.text(0.1,-0.5,miRNA_target,fontproperties='Consolas')
        ax3.axis('off')
        ###################################################################################################
        ## 绘制连接线
        xy = (Site, 0) 
        con = patches.ConnectionPatch( 
            xyA = xy, coordsA = ax2.transData, 
            xyB = xy, coordsB = ax1.transData, 
            linestyle="--", linewidth=1.5,edgecolor="red")
        ax2.add_artist(con) 
        text_xy = (0.5,3)
        con_text = patches.ConnectionPatch( 
            xyA = xy, coordsA = ax2.transData, 
            xyB = text_xy, coordsB = ax3.transData, 
            linestyle="--", linewidth=1.5,edgecolor="red")
        ax3.add_artist(con_text)
        final_xy = (0.5,1.5)
        con_final = patches.ConnectionPatch(xyA = text_xy, xyB = final_xy, 
                              coordsA=ax3.transData, coordsB=ax3.transData, arrowstyle="->", shrinkB=5,edgecolor="red")
        ax3.add_artist(con_final)
        # ax3.plot([xyA[0]], [xyA[1]], "o") 
        subtitle = f'{gene} -----> {miRNA}'
        plt.suptitle(subtitle,x=(8/12/2),y=1)
        plt.tight_layout()
        plt.savefig(Save_fig)
        print (f'Out Fig is {Save_fig}')
        '''
        这段代码定义了一个静态方法 Draw，用于绘制图表。传入的参数是一个 Series（序列），其中包含了需要绘制图表所需的数据。

        首先，从 Series 中提取出需要的数据，包括基因、miRNA、绘图的数据文件名、分数、Category、p 值、miRNA 剪切信息、保存的图表文件名以及保守结构域的起始和结束位置。
        
        然后，使用 pd.read_csv 方法读取绘图数据文件到数据帧 data 中。
        
        接下来，创建 T-plot 的主图表 ax1，绘制曲线图，添加垂直线段以标识切割位点 Site，设置 x 轴和 y 轴的标签和标题，计算 x 轴长度和分数相对于 y 轴长度的比值，添加注释文字。
        
        继续，创建保守结构域的图表 ax2，绘制矩形图以表示不同的结构域，隐藏坐标轴。
        
        然后，创建 miRNA 的剪切信息图表 ax3，使用 ax3.text 方法添加文字。
        
        最后，使用 patches.ConnectionPatch 创建连接线的对象，将连接线和箭头分别添加到相应的图表中。
        
        最后，设置图表的标题，调整布局，保存图表，并打印输出保存的图表文件名。
        '''
```
