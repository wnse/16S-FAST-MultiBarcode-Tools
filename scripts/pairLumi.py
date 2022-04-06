# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import pandas as pd
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import time
import os
import logging
import time
import numpy as np

def get_umi_dict(df):
    dup_idx = df[df.duplicated(subset=['umi_21','umi_11'], keep=False)].index
    for i in dup_idx:
        logging.info('duplicate umi paires:\t{}'.format(df.loc[i].to_list()))
    
    df_dropdup = df.drop(dup_idx)
    logging.info('umi paires:\t{}'.format(df_dropdup.shape[0]))
    
    df_dropdup['umi_pair'] = df_dropdup['umi_21']+"_"+df_dropdup['umi_11']
    umi_dict = df_dropdup.groupby('umi_pair')['Read name'].sum().to_dict()
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
                logging.info('exit for {} counts {} < {}' \
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

def dropRepeatUid(df, col1='umiID', col2='Barcode'):
    df_count = df.groupby(col1)[col2].apply(lambda x: x.value_counts().to_list())
    repeat_uid = df_count[df_count.apply(len)>1].index
    
    repeat_idx = df[df[col1].isin(repeat_uid)].index
    logging.info('remain {} different barcode:\t{}'.format(col1, df_count.shape[0]-len(repeat_uid)))
    logging.info('remain {} count:\t{}'.format(col1, df.shape[0]-len(repeat_idx)))
    logging.info('repeat {} different barcode:\t{}'.format(col1, len(repeat_uid)))
    logging.info('repeat {} count:\t{}'.format(col1, len(repeat_idx)))
    logging.info('\n{}'.format(df.loc[repeat_idx].groupby([col1,col2]).count()))
    return df.drop(repeat_idx), df.loc[repeat_idx]


def pairLumi(df):
#     df = pd.read_csv(file, index_col=0)
    umi_dict, umi1_dict, umi2_dict = get_umi_dict(df)
    out_list = define_UMI_ID(umi_dict, umi1_dict, umi2_dict)

    df_o = pd.DataFrame(out_list, columns=['umiID', 'umi_21', 'umi_11', 'Lreads'])
    df_o = pd.merge(df_o.set_index(['umi_21','umi_11']), df.set_index(['umi_21','umi_11']), 
                    left_index=True, right_index=True, how='left'
                   )
    logging.info('uIDs in total barcode:\t{}'.format(len(df_o['umiID'].unique())))
    logging.info('uID paires in total barcode:\t{}'.format(df_o.shape[0]))
    df_o_dropRepeat, df_o_repeat = dropRepeatUid(df_o)
    return df_o_dropRepeat.reset_index(), df_o_repeat.reset_index()

def plot_Lreads(df_list, outfile_list):
    for i, df in enumerate(df_list):
        outfile = outfile_list[i] + '_uID_umiPaire_Lreads.png'
        tmp_data = df.set_index('umiID')['Lreads']
        plot_multiUMIreads(tmp_data, outfile=outfile)
        
        outfile = outfile_list[i] + '_uID_barcode_Lreads.png'
        plot_umiID_Lreads(df, outfile=outfile)
        
        outfile = outfile_list[i] + '_uID_max_Lreads.png'
        tmp_data = df.drop_duplicates(subset=['umiID'], keep='first').sort_values(by='Lreads', ascending=False)
        plot_uID_each_Lreads(tmp_data, outfile=outfile)
    
    
# logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
# out_dir = '../test/out_dir_2/'
# File_Tag='test'
# # L_umi_info_file = os.path.join(out_dir, File_Tag+'_L_umi.csv')
# # L_umi_sta = os.path.join(out_dir, File_Tag+'_L_umi_sta.csv')
# # df_uID, df_uID_dropRepeat = pairLumi(L_umi_info_file)

# L_uID_info_file = os.path.join(out_dir, File_Tag+'_L_uID_dropRepeat.csv')
# # df_uID.to_csv(L_uID_info_file)
# L_uID_info_Repeat_file = os.path.join(out_dir, File_Tag+'_L_uID_Repeat.csv')
# # df_uID_dropRepeat.to_csv(L_uID_info_file)
# df_uID = pd.read_csv(L_uID_info_file, index_col=0)
# df_uID_repeat = pd.read_csv(L_uID_info_Repeat_file, index_col=0)
# plot_Lreads([df_uID, df_uID_repeat], [os.path.join(out_dir, File_Tag+'_'+i) for i in ['dropRepeat', 'Repeat']])

# %%
# def plot_uID_each_Lreads(df, outfile='uID_max_Lreads.png'):
#     fig, ax = plt.subplots()
#     tmp_data = df.reset_index().drop('index', axis=1)
#     sns.lineplot(y='Lreads', x='index',data=tmp_data.reset_index(), ax=ax, label='Total')
#     idx = df['Barcode'].value_counts().index[:10]
#     for i in idx:
#         tmp_data = df[df['Barcode']==i].reset_index().drop('index', axis=1)
#         sns.lineplot(y='Lreads', x='index',data=tmp_data.reset_index(), ax=ax, label=i)
#     plt.xlabel('range of uID')
#     plt.tight_layout()
#     plt.show()
# #     plt.tight_layout()
# #     fig.savefig(outfile)

# %%
def plot_uID_each_Lreads(df, outfile='uID_max_Lreads.png'):
    fig, ax = plt.subplots()
    tmp_data = df.reset_index().drop('index', axis=1)
    plt.plot(list(range(0, tmp_data.shape[0])), tmp_data['Lreads'], label='Total')
    idx = df['Barcode'].value_counts().index[:10]
    tmp_no = 0
    y_max = tmp_data['Lreads'].max()
    for i in idx:
        tmp_data = df[df['Barcode']==i].reset_index().drop('index', axis=1)
        plt.plot(range(tmp_no, tmp_no+tmp_data.shape[0]), tmp_data['Lreads'].to_list(), label=i)
        tmp_no += tmp_data.shape[0] 
        plt.vlines(tmp_no, 0, y_max, linestyles='dashed', colors='grey', alpha=0.3)
    plt.legend()
    plt.xlabel('range of uID')
    plt.tight_layout()
#     plt.show()
    fig.savefig(outfile)


# %% tags=[]
def plot_multiUMIreads(data, max_col=10, outfile='uID_umiPaire_reads.png'):
    df_tmp = data.groupby(level=0).apply(lambda x:x.to_list())
    df_tmp = df_tmp.to_list()
    max_len = max([len(i) for i in df_tmp])
    if max_len>max_col:
        max_len = max_col
    np_o = np.zeros([len(df_tmp),max_len])
    for i1,v1 in enumerate(df_tmp):
        for i2,v2 in enumerate(v1):
            if i2 < max_len:
                np_o[i1,i2] = v2
            else:
                break
    np_o = np.sort(np_o)[np.lexsort(np.sort(np_o).T, -1)][:,::-1]
    fig, ax = plt.subplots()
    sns.heatmap(np_o, cmap="YlGnBu", ax=ax)
    plt.tight_layout()
    fig.savefig(outfile)
    

# df_tmp = df_uID[df_uID.duplicated(subset='umiID', keep=False)].set_index('umiID')['Lreads']
# plot_multiUMIreads(df_tmp)
# plt.show()

# df_tmp = df_uID_dropRepeat[df_uID_dropRepeat.duplicated(subset='umiID', keep=False)].set_index('umiID')['Lreads']
# plot_multiUMIreads(df_tmp)
# plt.show()

# %% tags=[]
def plot_umiID_Lreads(df_umiID, min_umiID=100, outfile='uID_barcode_reads.png'):
    df_uid_count = df_umiID.drop_duplicates(subset=['Barcode','umiID']).groupby('Barcode')['umiID'].count()
    df_uid_Lreads = df_umiID.groupby('Barcode')['Lreads'].sum()
    df_tmp = pd.merge(df_uid_count, df_uid_Lreads, left_index=True, right_index=True)
    df_tmp = df_tmp.sort_values(by=['umiID','Lreads'], ascending=False)
    df_tmp['Lreads/umiID'] = df_tmp['Lreads']/df_tmp['umiID']
    df_tmp = df_tmp[df_tmp['umiID']>min_umiID]

    fig, ax1 = plt.subplots(1,1)
    plt.ylabel('#')
    
    ax1.plot(df_tmp['umiID'], c='b', label='umiID')
    ax1.plot(df_tmp['Lreads'], c='g', label='Lreads')
    plt.legend()
    plt.xticks(rotation=90)
    
    ax2 = ax1.twinx()
    plt.ylabel('Lreads/uID')
    ax2.plot(df_tmp['Lreads/umiID'], c='r', label='Lreads/umiID')
    
    plt.legend(loc='center right')
    plt.tight_layout()
    fig.savefig(outfile)

# plot_umiID_Lreads(df_uID)
# plot_umiID_Lreads(df_uID_dropRepeat)
