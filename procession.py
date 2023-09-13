#!/usr/bin/env python
import pandas as pd
import argparse 


parser = argparse.ArgumentParser(
    description=f"""
    目的1:对NerveInfect项目的中间文件进行处理,主要是添加“检出病原体”中遗漏的部分、“已知结果rpk”.
        需要注意的是,输入文件是procession.txt,在此之前要检查“已知结果来源”、“已知结果” ===> 这些需要人工excel核对之后生成procession.txt。
        同时也要检查procession.txt中的'已知结果'和'检出结果'中病原体是否对应一致的,否则统计出错：
        例如 [黑曲霉、黄曲霉、土曲霉、布鲁菌、乙型流感病毒BV ==> 黑曲霉复合群、黄曲霉复合群、土曲霉复合群、布鲁菌属、乙型流感病毒]等
        
    目的2:完成对处理的"汇总"统计,得到'培养物检测限’。目前针对的是'培养物’的样本。
        建议处理完之后需要核对

    特别要对procession.txt病原体的名称进行核对,例如,不同批次,不同亚型,检出病原命名不规范等
        """,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f1",type=str,
                    help= " procession file")
args = parser.parse_args()


########################################################################################脚本1
def procession(procession_file):
# 读取人工处理的procession.txt文件并将数据存储到DataFrame中
    with open(procession_file, 'r') as file:
        lines = file.readlines()

    columns = lines[0].strip().split('\t') # 读取列名

    
    data = []   # 读取数据
    for line in lines[1:]:
        parts = line.strip().split('\t')
        data.append(parts)

    data_df = pd.DataFrame(data, columns=columns)
    result = []

    # 遍历每个样本类型
    for sample_type in data_df['样本编号'].unique():
        sample_data = data_df[data_df['样本编号'] == sample_type]
        known_results = set(sample_data['已知结果'].apply(lambda x: x.strip().split("|")[0]).unique())
        detected_pathogens = set(sample_data['检出病原体'].apply(lambda x:x.strip().split("|")[0]).unique())

        # 找到未检出的已知结果
        missing_results = known_results - detected_pathogens  #两个集合之间的差运算，结果是known_results存在，而detected_pathogens不存在的结果

        if missing_results:
            for missing_result in missing_results:
                missing_row = sample_data.iloc[0].copy()  # 复制第一行数据
                missing_row['检出病原体'] = str(missing_result).split("|")[0] + '|漏'
                missing_row['read'] = 0
                missing_row['rpk'] = 0
                missing_row['过滤标志'] = 'na'
                missing_row['扩增子覆盖度'] = 'na'
                missing_row['扩增子'] = "na"
                missing_row['扩增子read'] = 0
                missing_row['扩增子rpk']=0
                missing_row['扩增子rate'] = 0 
                missing_row['最佳比对率']="na"
            
                result.append(missing_row)
        result_df1 = pd.DataFrame(result,columns=columns) 

# 将新增的遗漏情况合并到原始数据
    result_df = pd.concat([data_df,result_df1],ignore_index=True)
    return result_df


def fix(result_df):
    '对上一步的已知结果进一步修订,包括("内参")和 (know_bact 和 detect_bact不匹配)的情况,赋值"/"'
    for sample_type in result_df['样本编号'].unique():
        sample_data = result_df[result_df['样本编号'] == sample_type]
        for index, row in sample_data.iterrows():
            know_bact = row['已知结果'].split("|")[0].strip()
            known_result = row['已知结果']
            know_rpk = known_result.split("|")[1].strip() if '|' in known_result else "/"  #这里还需要判断“|”不存在的那种情况。
            detect_bact = row['检出病原体'].split("|")[0].strip()

            if row['过滤标志'] == "内参":
                result_df.at[index, '已知结果rpk'] = "/"   #由于迭代的值没有输出出来，因此不适用：row['已知结果rpk'] == "/"。相反，是对原始的DF进行修改，因此应该使用at或iat方法来设置特定行和列的值，而不是直接使用=。
            else:
                if know_bact == detect_bact:
                    result_df.at[index, '已知结果rpk'] = know_rpk
                else:
                    result_df.at[index, '已知结果rpk'] = "/"
            ##不需要使用return函数，因为只是在原始的 result_df 上进行修复操作
####################################################################################








###################################################################################脚本2
def process_data_frame(xlxs_file):
    '目的:搭建统计检验限的表格框架'
    '大概逻辑:先读取“汇总”表格;筛选“培养物(或其它)”的表格部分 ==> Uniq 已知病原 作为index;==> Uniq 非“\”的rpk作为column'
    data_df = pd.read_excel(io = xlxs_file)
    data_df['已知病原'] = data_df["已知结果"].apply(lambda x : x.strip().split("|")[0] if (x and "|" in x) else x)  #新添一列

    kown_path_type = []
    kown_rpk_type = []
    for index,row in data_df.iterrows():
        if str(row['已知结果来源']) == "培养物":   ##这里可能经常要修改,因为有新冠的“核酸标准物质”这一栏
            path_type = str(row['已知病原'])
            kown_path_type.append(path_type)

            if "/" not in str(row["已知结果rpk"]):
                rpk_type = str(row["已知结果rpk"]).strip()
                kown_rpk_type.append(rpk_type)

    kown_path_type = list(set(kown_path_type))
    kown_rpk_type = list(set(kown_rpk_type))

    result_pd = pd.DataFrame(index=kown_path_type,columns=kown_rpk_type)    #创建“检测限”的表格框架
    result_pd = result_pd[sorted(result_pd.columns, key=lambda x: int(x))]  #列名排序
    return result_pd,data_df




def fill_data_frame (result_pd,data_df):
    '目的:对上一步搭建的检验限的表格框架进行填充'
    '大概逻辑:先对data_df滤掉 非培养物 rpk为/ 的列表部分;==> 按“已知病原”的每一种病原 [result_type] 去截取表格,筛选部分列并去重'
    '==> 对上一步截取的表格,按“已知结果rpk”中的每一类rpk [rpk_type] 去截取表格'
    '==> 以i;j;rpk记录检测限 ==> 按照[result_type, rpk_type]将其写入到result_pd'
    data_df = data_df[data_df['已知结果来源'] == "培养物"]  ##这里可能经常要修改,因为有新冠的“核酸标准物质”这一栏
    data_df = data_df[data_df['已知结果rpk'] != "/"]

    for result_type in data_df['已知病原'].unique():  ##按“已知病原”去截取表格段落
        result_data = data_df[data_df['已知病原'] == result_type]
        selected_columns = ['已知结果来源', '样本编号', '已知病原','已知结果', '检出病原体', 'read', 'rpk', '过滤标志', '已知结果rpk']
        df_subset = result_data[selected_columns].drop_duplicates()   #去重   

        
        for rpk_type in df_subset['已知结果rpk'].unique():  ##按“已知结果rpk”去截取df_subset的表格段落
            i=0;j=0;rpk = []    #i:应当出现的个数(样本数);j:实际出现的个数
            result_rpk_data = df_subset[df_subset['已知结果rpk'] == rpk_type]

            for index,row in result_rpk_data.iterrows():
                detected_pathogens = str(row["检出病原体"]).strip()
                row_index =str(row['已知病原']).strip()
                if  row_index == detected_pathogens.split("|")[0]:
                    i = i + 1
                    if "报" in detected_pathogens.split("|")[1]:
                        j = j + 1
                        rpk.append(str(row['rpk']))
                    else:
                        j = j+0
                        rpk.append(str(row['rpk']))
                else:
                    pass

            output_str = f"{j}/{i}({'-'.join(rpk)})"  ##列表转换
            result_pd.loc[str(result_type),str(rpk_type)] = output_str
            result_pd = result_pd.fillna('')  ##将NaN

            # print (output_str)
    return result_pd
####################################################################################




def main():
    result_df = procession(args.f1)
    fix(result_df)
    result_df.to_excel('result.xlsx', index=False)
    result_pd,data_df = process_data_frame('result.xlsx')
    lod_pd = fill_data_frame(result_pd=result_pd,data_df=data_df)
    lod_pd.to_excel('LOD.xlsx', index=True)

if __name__ == "__main__":
    main()

