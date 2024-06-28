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
                    T_plot_fig_message_file = os.path.basename(a) + '_draw.pdf'
                    # print (T_plot_message_file)
                    miRNA,gene,Site = line.strip().split('/')[1].split('_')[0:3]
                    ## 提取基因以及miRNA信息
                # elif (line.startswith("SiteID:")):
                #     _,gene,miRNA = line.strip().split(':')
                #     miRNA = 'miR' + miRNA
                #     gene = gene.strip()
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
        Message_df.set_index("gene")
        ## 匹配保守结构域
        # Conserved_gene = set(Conserved_df['Conserved_Gene'])
        # for gene in Conserved_gene:
        #     Select_df = Conserved_df[Conserved_df['Conserved_Gene']==gene]
        #     Series = Message_df.loc[Message_df['gene']==gene]
        #     if Series.empty:
        #         continue
        #     for i in range(len(Select_df)):
        #         Conserved_start = f'Conserved{i+1}_start'
        #         Conserved_end = f'Conserved{i+1}_end'
        #         print (Series)
        #         Series.loc[:,Conserved_start] = Select_df.iloc[i]['Conserved_start']
        #         Message_df.loc[gene:,Conserved_end] = Select_df.iloc[i]['Conserved_end']
        #         Message_df.loc[Message_df['gene']==gene].Conserved_end = Select_df.iloc[i]['Conserved_end']
        #         Message_df.loc[Message_df['gene']==gene, Conserved_start] = Select_df.iloc[i]['Conserved_start']
        # Message_df.to_csv(Outfile,index=False)
        Conserved_gene = set(Conserved_df['Conserved_Gene'])
        for gene in Conserved_gene:
            Select_df = Conserved_df[Conserved_df['Conserved_Gene']==gene]
            if gene in Message_df['gene'].values:
                Series = Message_df.loc[Message_df['gene']==gene]
                if Series.empty:
                    continue
                for i in range(len(Select_df)):
                    Conserved_start = f'Conserved{i+1}_start'
                    Conserved_end = f'Conserved{i+1}_end'
                    Message_df.loc[Message_df['gene']==gene,Conserved_start] = Select_df.iloc[i]['Conserved_start']
                    Message_df.loc[Message_df['gene']==gene,Conserved_end] = Select_df.iloc[i]['Conserved_end']
        Message_df.to_csv(Outfile, index=False)
        return Message_df
    @staticmethod
    ## 绘图
    def Draw(Series):
        # 读取整理的数据
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
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
        plt.text(0.6,0.8,ax1_annotation,transform=ax1.transAxes,fontsize=10)
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
        ax3.text(0.1,-0.5,miRNA_target,fontproperties=font_manager.FontProperties(fname='/home/xurong_fan/anaconda3/lib/python3.9/site-packages/matplotlib/mpl-data/fonts/ttf/consola.ttf'))
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

## 读取文件
a = Dragome.Read_result_file('flower.txt','AP2_cds_cd_search.txt')

## 筛选
Select_df = a[ (a['score']<=5) & (a['Category']<=2)]
Select_df = Select_df.reset_index(drop=True)

Select_gene = set(Select_df['gene'])
Target_df = pd.DataFrame
Target_gene =[]
with open ('Lch.id','r') as f:
    for line in f:
        gene = line.strip()
        for index, row in Select_df.iterrows():
            if row['gene']  == gene:
                Target_gene.append(index)

## 绘图
for index in Target_gene:
    Series = Select_df.iloc[index,:]
    Dragome.Draw(Series)