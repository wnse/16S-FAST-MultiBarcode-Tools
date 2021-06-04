import pandas as pd
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import time
import os

def get_overlap_umi(df_merge_f, umi='umi_10', col_barcode='R1_Barcode', col_count='Read name'):
    tmp = df_merge_f.reset_index().groupby([col_barcode, umi])[col_count].count().reset_index()
    out_dict={}
    for i in combinations(tmp[col_barcode].unique(), 2):
        barcode1 = i[0]
        barcode2 = i[1]
        key = "{} vs {}".format(barcode1, barcode2)
        value = (tmp[tmp[col_barcode]==barcode2][umi].isin(tmp[tmp[col_barcode]==barcode1][umi])).sum()
        out_dict[key] = value
    return df_merge_f[col_barcode].value_counts(), tmp[col_barcode].value_counts(), pd.DataFrame.from_dict(out_dict, orient='index')

def get_umi1_dict(df):
    umi_dict = (df['umi_10']+"_"+df['umi_11']).value_counts().to_dict()
    umi1_dict = {}
    umi2_dict = {}
    for i,v in umi_dict.items():
        u1 = i.split('_')[0]
        u2 = i.split('_')[1]
        umi1_dict[u1] = umi1_dict.get(u1, {})
        umi1_dict[u1][u2] = v
        umi2_dict[u2] = umi2_dict.get(u2, {})
        umi2_dict[u2][u1] = v
    return umi_dict, umi1_dict, umi2_dict

def define_UMI_ID(UMI_dict, UMI1_dict, UMI2_dict, cutoff=0.5, min_counts=1):
    ''' 根据连接文库序列的umi，定义umiID
    参数：
        UMI_dict: 连接文库所有umi配对关系及数量（type=dict）
        UMI1_dict: 连接文库所有umi1对应的umi2及数量
        UMI2_dict: 连接文库所有umi2对应的umi1及数量
        cutoff: 可以作为umiID的配对关系所需要大于第一配对关系数量的最小比例
        min_counts: 配对关系的最小数量
    返回：
        out_list: [umiID, umi1, umi2, counts]
    '''
    out_list = []
    UMI_list = sorted(UMI_dict.items(), key=lambda UMI_dict: UMI_dict[1], reverse=True)
    for idx_count in UMI_list:
        if len(UMI1_dict) == 0:
            break
        (u1, u2) = idx_count[0].split('_')
        if u1 in UMI1_dict.keys():
            if len(UMI1_dict[u1]) == 0:
                del UMI1_dict[u1]
                continue
            if u2 not in UMI2_dict.keys():
                del UMI1_dict[u1][u2]
                continue
            if UMI1_dict[u1][u2] < min_counts:
                logging.info(' exit for {} counts {} < {}' \
                             .format(idx_count[0], UMI1_dict[u1][u2], min_counts))
                break
            out_list.append([idx_count[0], u1, u2, UMI1_dict[u1][u2]])
            list1 = list(UMI2_dict[u2].keys())
            list2 = list(UMI1_dict[u1].keys())
            max_counts = UMI1_dict[u1][u2]
            del UMI1_dict[u1]
            del UMI2_dict[u2]
            for deu in list1:
                if deu in UMI1_dict.keys():
                    if max(UMI1_dict[deu], key=UMI1_dict[deu].get) == u2:
                        if UMI1_dict[deu][u2] >= cutoff * max_counts:
                            out_list.append([idx_count[0], deu, u2, UMI1_dict[deu][u2]])
                        del UMI1_dict[deu]
            for deu in list2:
                if deu in UMI2_dict.keys():
                    if max(UMI2_dict[deu], key=UMI2_dict[deu].get) == u1:
                        if UMI2_dict[deu][u1] >= cutoff * max_counts:
                            out_list.append([idx_count[0], u1, deu, UMI2_dict[deu][u1]])
                        del UMI2_dict[deu]
    return out_list

def get_umi_barcode(df, u1, u2, col1='umi_10', col2='umi_11', col_barcode='R1_Barcode'):
#     u1 = out_list[0][1]
#     u2 = out_list[0][2]
    return df[(df[col1]==u1)&(df[col2]==u2)][col_barcode].unique()

def get_umi_barcode_total(df, umiID_list):
    umiID_barcode=[]
    df_barcode = df.groupby(['umi_10', 'umi_11', 'R1_Barcode'])['Read name'].count()
    t0 = time.perf_counter()
    for n, i in enumerate(umiID_list):
        u1 = i[1]
        u2 = i[2]
        tmp = df_barcode.loc[(u1,u2)].index.to_list()
#         tmp = get_umi_barcode(df, u1, u2)
        if len(tmp)>1:
            print(i, tmp)
        else:
            umiID_barcode.append(tmp[0])
        print('{:.2f}%\t{:.2f}s'.format(n/len(umiID_list)*100, time.perf_counter()-t0), end='\r')
    return umiID_barcode

def main_info(file):
#     file = 'umi_info/2.umi_info.txt'
    df = pd.read_csv(file)
    print('{}\t{}'.format(file, df.shape[0]))
    tmp1, tmp2, tmp3 = get_overlap_umi(df)
    tmp1.index = pd.MultiIndex.from_product([['Reads'],tmp1.index])
    tmp2.index = pd.MultiIndex.from_product([['Left UMIs'],tmp2.index])
    tmp3.index = pd.MultiIndex.from_product([['Left UMIs Overlap'],tmp3.index])
    df_o = pd.concat([tmp1, tmp2, tmp3])
    print(df_o)
    tmp1, tmp2, tmp3 = get_overlap_umi(df, umi='umi_11')
    tmp2.index = pd.MultiIndex.from_product([['Right UMIs'],tmp2.index])
    tmp3.index = pd.MultiIndex.from_product([['Right UMIs Overlap'],tmp3.index])
    print(tmp2)
    print(tmp3)
    df_o = pd.concat([df_o, tmp2, tmp3])
    
    umi_dict, umi1_dict, umi2_dict = get_umi1_dict(df)
    out_list = define_UMI_ID(umi_dict, umi1_dict, umi2_dict)
    umiID_barcode = get_umi_barcode_total(df, out_list)
    
    tmp1 = pd.DataFrame(umiID_barcode)[0].value_counts()
    tmp1.index = pd.MultiIndex.from_product([['UMI_IDs'],tmp1.index])
    tmp2 = pd.DataFrame(out_list, umiID_barcode).reset_index().groupby('index')[3].sum()
    tmp2.index = pd.MultiIndex.from_product([['UMI_IDs Reads'],tmp2.index])
    print(tmp1)
    print(tmp2)
    df_o = pd.concat([df_o, tmp1, tmp2])
    
    return df_o

d = 'umi_info/'
df_result = pd.DataFrame()
for i in os.listdir(d):
    file = os.path.join(d,i)
    df_o = main_info(file)
    if df_result.empty:
        df_result = df_o.rename(columns={0:i})
    else:
        df_result = pd.merge(df_result, df_o.rename(columns={0:i}), left_index=True, right_index=True, how='outer')
df_result.to_excel('get_L_umi_info.xlsx')
