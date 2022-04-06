# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.9.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import matplotlib.pyplot as plt
import seaborn as sns
from pairLumi import dropRepeatUid
import logging

# %%
import pandas as pd
import argparse
import logging
import os
import numpy as np

def corresponding_seq2umi(adapter_info_file, primer_info_file, umi_set):
    '''处理拼接文库数据R2的接头，确定UMI，并序列ID与连接文库UMI进行关联
    参数：
        adapter_info_file: adapter序列文件
        primer_info_file: primer序列文件
        umi_dict: 连接文库确定的umi与umiID关系（type=dict）
    返回：
        aUMI: 序列ID与umi的对应关系（type=dict）
        sta_list: 统计信息（type=list)
    '''
#     sta_list = []
    total_seq = 0
    seq_type = {}
    umi_seq = set()
    umi_seq_paired = []
    umi_seq_unpaired = []
    n = 0
    for i in open(adapter_info_file):
        total_seq += 1
        tmp = i.split('\t', )
        seq_id = tmp[0].split(' ')[0]
        if tmp[1] == str(-1):
            continue
        else:
            n += 1
            seq_type[seq_id] = tmp[7]

    logging.info('A reads:\t{:d}'.format(total_seq))
    logging.info('A reads with adapters:\t{}'.format(n))

    n = 0
    m = 0
    for i in open(primer_info_file):
        tmp = i.split('\t')
        seq_id = tmp[0].split(' ')[0]
        if seq_id in seq_type.keys():
            if tmp[1] == str(-1):
                continue
            else:
                n += 1
                seq = tmp[6]
                if len(seq)<=20 and len(seq)>=10:
                    if seq_type[seq_id] == tmp[7][0]:
                        m += 1
                        if seq in umi_set:
                            umi_seq.add(seq)
                            umi_seq_paired.append([seq_id, seq, tmp[7]])
                        else:
                            umi_seq_unpaired.append([seq_id, seq, tmp[7]])                                                     
    logging.info('A reads with linkers:\t{}'.format(n))
    logging.info('A reads with umi:\t{}'.format(m))
    logging.info('A reads with umi In L:\t{}'.format(len(umi_seq_paired)))
    logging.info('A reads with umi NotIn L:\t{}'.format(len(umi_seq_unpaired)))
    df_umi_paired = pd.DataFrame(umi_seq_paired, columns=['AreadID', 'umi', 'Barcode'])
    df_umi_unpaired = pd.DataFrame(umi_seq_unpaired, columns=['AreadID', 'umi', 'Barcode'])
    return df_umi_paired, df_umi_unpaired


# %%
# def dropRepeatUid(df, col1='umiID', col2='Barcode'):
#     df_count = df.groupby(col1)[col2].apply(lambda x: x.value_counts().to_list())
#     repeat_uid = df_count[df_count.apply(len)>1].index
    
#     repeat_idx = df[df[col1].isin(repeat_uid)].index
#     logging.info('repeat {} different barcode: {}'.format(col1, len(repeat_uid)))
#     logging.info('repeat {} count: {}'.format(col1, len(repeat_idx)))
#     logging.info('\n{}'.format(df.loc[repeat_idx].groupby([col1,col2]).count()))
#     return df.drop(repeat_idx), df.loc[repeat_idx]

# %%
def staAumi(A2_cut_adapter_info, A2_umi_file_info, umi_set):
    df_umi_paired, df_umi_unpaired = corresponding_seq2umi(A2_cut_adapter_info, A2_umi_file_info, umi_set)
    logging.info('umi In L:\t')
    df_umi_paired_dropRepeat, df_umi_paired_repeat = dropRepeatUid(df_umi_paired, col1='umi', col2='Barcode')
    logging.info('umi Not In L:\t')
    df_umi_unpaired_dropRepeat, df_umi_unpaired_repeat = dropRepeatUid(df_umi_unpaired, col1='umi', col2='Barcode')
    return df_umi_paired_dropRepeat, df_umi_unpaired_dropRepeat

# %%
# logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level='INFO')
# out_dir = '../test/out_dir'
# File_Tag='test'

# A2_cut_adapter = os.path.join(out_dir,File_Tag+'_cut_adapter_P.fq.gz')
# A2_cut_adapter_info = A2_cut_adapter + '.cutadapt.info.file'
# A2_umi_file = os.path.join(out_dir,File_Tag+'_cut_linker_P.fq')
# A2_umi_file_info = A2_umi_file + '.cutadapt.info.file'


# L_uID_info_file = os.path.join(out_dir, File_Tag+'_L_uID_dropRepeat.csv')
# df_uID = pd.read_csv(L_uID_info_file, index_col=0)

# umi_set = set(df_uID['umi_21'].to_list() + df_uID['umi_11'].to_list())

# umi_seq_paired, umi_seq_unpaired = corresponding_seq2umi(A2_cut_adapter_info, A2_umi_file_info, umi_set)
# df_umi_paired = pd.DataFrame(umi_seq_paired, columns=['AreadID', 'umi', 'Barcode'])
# df_umi_unpaired = pd.DataFrame(umi_seq_unpaired, columns=['AreadID', 'umi', 'Barcode'])

# df_umi_paired_dropRepeat, df_umi_paired_repeat= dropRepeatUid(df_umi_paired, col1='umi', col2='Barcode')
# A_umi_info_file = os.path.join(out_dir, File_Tag+'_A_umi_dropRepeat.csv')
# df_umi_paired_dropRepeat.to_csv(A_umi_info_file)

# df_umi_unpaired_dropRepeat, df_umi_unpaired_repeat= dropRepeatUid(df_umi_unpaired, col1='umi', col2='Barcode')

# plot_each_count(df_umi_paired_dropRepeat, title='Total Paired')
# plot_each_count(df_umi_paired_dropRepeat[df_umi_paired_dropRepeat['Barcode'].str.contains('F')], title='Total Paired F')
# plot_each_count(df_umi_paired_dropRepeat[df_umi_paired_dropRepeat['Barcode'].str.contains('R')], title='Total Paired R')

# plot_each_count(df_umi_unpaired_dropRepeat, title='Total UnPaired')
# plot_each_count(df_umi_unpaired_dropRepeat[df_umi_unpaired_dropRepeat['Barcode'].str.contains('F')], title='Total UnPaired F')
# plot_each_count(df_umi_unpaired_dropRepeat[df_umi_unpaired_dropRepeat['Barcode'].str.contains('R')], title='Total UnPaired R')
